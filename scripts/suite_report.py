#!/usr/bin/env python3
"""Aggregate the per-item reports of one suite run into a single verdict.

Every item in the suite -- a case script or a python checker -- already writes
its own JSON with the same shape: a list of named checks, each PASS/WARN/FAIL/
ERROR/SKIP/INFO. That is the right granularity while debugging one item and the
wrong one for answering "can we ship this build", which is the question a
release run actually has to answer.

This turns the pile into one table and one exit code. Two details matter:

  an item that crashed is not an item that failed
      A case that could not run -- missing data, no license proxy, no JDK 8 --
      says nothing about the build. Reporting it as a failure trains people to
      ignore red, and reporting it as a pass hides that the scenario was never
      tested. It gets its own status, and the summary counts it separately.

  a missing report is itself a finding
      If an item exits non-zero before writing JSON it is usually a crash in
      the harness rather than in the product. Rather than silently dropping the
      row, the item is listed with whatever the exit code was, so the run id
      still accounts for every item that was supposed to execute.

Exit codes: 0 everything passed, 1 something failed, 2 something could not run.
"""

import argparse
import json
import os
import sys

# Same rank order as lib/report.py, so a suite verdict can never be better than
# the worst check inside the items it summarises.
RANK = {"PASS": 0, "INFO": 0, "SKIP": 1, "WARN": 2, "FAIL": 3, "ERROR": 4}

COLOR = {"PASS": "\033[32m", "WARN": "\033[33m", "SKIP": "\033[36m",
         "FAIL": "\033[31m", "ERROR": "\033[35m"}
RESET = "\033[0m"


def status_from_exit(code):
    return {0: "PASS", 1: "FAIL"}.get(code, "ERROR")


def load_item(spec):
    """Parse a 'name|kind|exit|wall|json' record and fold in its report."""
    name, kind, code, wall, path = spec.split("|", 4)
    item = {
        "name": name,
        "kind": kind,
        "exit_code": int(code),
        "wall_s": float(wall),
        "status": status_from_exit(int(code)),
        "report": path,
        "counts": {},
        "problems": [],
    }

    if not path or not os.path.exists(path):
        item["note"] = "no report written"
        return item

    try:
        with open(path) as handle:
            data = json.load(handle)
    except (OSError, ValueError) as exc:
        item["note"] = "unreadable report: %s" % exc
        return item

    checks = data.get("checks", [])
    for check in checks:
        state = check.get("status", "INFO")
        item["counts"][state] = item["counts"].get(state, 0) + 1
        if state in ("FAIL", "ERROR"):
            item["problems"].append("%s: %s" % (check.get("name", "?"),
                                                check.get("message", "")))

    # Trust the report over the exit code: a case that aborts early still calls
    # case_end, and the checks it did record are the useful part.
    worst = max((RANK.get(c.get("status", "INFO"), 0) for c in checks), default=0)
    for state, rank in RANK.items():
        if rank == worst and state in ("FAIL", "ERROR", "WARN"):
            item["status"] = state
            break
    if worst <= 1 and item["exit_code"] == 0:
        item["status"] = "SKIP" if worst == 1 else "PASS"

    # An item may declare its own verdict. Take it only when it is worse than
    # what the checks show: a harness that forgets to account for a status --
    # a skipped case still reporting PASS, say -- must not be able to talk the
    # summary into calling an untested scenario green.
    declared = data.get("status")
    if declared and RANK.get(declared, 0) > RANK.get(item["status"], 0):
        item["status"] = declared

    item["duration_s"] = data.get("duration_s")
    item["steps"] = data.get("steps", [])
    return item


def summarise(items):
    counts = {}
    for item in items:
        counts[item["status"]] = counts.get(item["status"], 0) + 1
    worst = max((RANK.get(i["status"], 0) for i in items), default=0)
    verdict = "PASS"
    for state, rank in RANK.items():
        if rank == worst and state in ("FAIL", "ERROR", "WARN", "SKIP"):
            verdict = state
            break
    return verdict, counts


def print_table(items, verdict, counts, run_id, color):
    def paint(state):
        if not color:
            return state
        return "%s%s%s" % (COLOR.get(state, ""), state, RESET)

    width = max([len(i["name"]) for i in items] + [4])
    print()
    print("=" * (width + 46))
    print("suite run %s" % run_id)
    print("=" * (width + 46))
    print("  %-*s  %-7s  %9s  %s" % (width, "item", "status", "wall", "detail"))
    for item in items:
        detail = ""
        if item["problems"]:
            detail = item["problems"][0]
            if len(item["problems"]) > 1:
                detail += "  (+%d more)" % (len(item["problems"]) - 1)
        elif item.get("note"):
            detail = item["note"]
        elif item["counts"]:
            detail = " ".join("%s=%d" % (k.lower(), v)
                              for k, v in sorted(item["counts"].items()))
        print("  %-*s  %-7s  %8.1fs  %s"
              % (width, item["name"], paint(item["status"]), item["wall_s"],
                 detail[:80]))

    print("-" * (width + 46))
    tally = " ".join("%s=%d" % (k.lower(), v) for k, v in sorted(counts.items()))
    print("  verdict: %s   (%s)" % (paint(verdict), tally))

    # Repeat the failures in full underneath the table: the one-line detail
    # column is a pointer, not a diagnosis.
    for item in items:
        if item["status"] in ("FAIL", "ERROR") and item["problems"]:
            print()
            print("  %s:" % item["name"])
            for problem in item["problems"]:
                print("      %s" % problem)
            print("      report: %s" % item["report"])
    print()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--item", action="append", default=[], metavar="SPEC",
                        help="'name|kind|exit|wall|json' record (repeatable)")
    parser.add_argument("--run-id", default="", help="identifier for this run")
    parser.add_argument("--tier", default="", help="tier that was requested")
    parser.add_argument("--out", metavar="PATH", help="write suite JSON here")
    parser.add_argument("--no-color", action="store_true")
    args = parser.parse_args(argv)

    if not args.item:
        parser.error("nothing to aggregate: pass at least one --item")

    items = [load_item(spec) for spec in args.item]
    verdict, counts = summarise(items)

    color = not args.no_color and sys.stdout.isatty()
    print_table(items, verdict, counts, args.run_id, color)

    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as handle:
            json.dump({"tool": "suite", "schema_version": 1,
                       "run_id": args.run_id, "tier": args.tier,
                       "verdict": verdict, "counts": counts,
                       "items": items}, handle, indent=2)

    if any(i["status"] == "ERROR" for i in items):
        return 2
    if any(i["status"] == "FAIL" for i in items):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
