"""Shared check/report plumbing for the validation scripts.

Every script emits the same JSON shape so the release runner can aggregate
results across cases without knowing what each script does:

    {"tool": ..., "status": ..., "checks": [{"name", "status", "message", "detail"}]}

Exit codes are uniform too: 0 pass, 1 fail, 2 error.
"""

from __future__ import annotations

import datetime
import json
import os
import sys
import time

PASS = "PASS"
WARN = "WARN"
FAIL = "FAIL"
ERROR = "ERROR"
SKIP = "SKIP"
INFO = "INFO"

_SEVERITY = {INFO: 0, SKIP: 0, PASS: 1, WARN: 2, FAIL: 3, ERROR: 4}

EXIT_PASS = 0
EXIT_FAIL = 1
EXIT_ERROR = 2

_COLORS = {PASS: "\033[32m", WARN: "\033[33m", FAIL: "\033[31m", ERROR: "\033[35m"}


def _use_color(stream):
    return stream.isatty() and not os.environ.get("NO_COLOR")


def _paint(status, stream):
    if not _use_color(stream) or status not in _COLORS:
        return status
    return "%s%s\033[0m" % (_COLORS[status], status)


class Report(object):
    def __init__(self, tool, **context):
        self.tool = tool
        self.context = context
        self.checks = []
        self.started_at = time.time()

    def add(self, name, status, message, detail=None):
        self.checks.append(
            {"name": name, "status": status, "message": message, "detail": detail or {}}
        )
        return status

    def check(self, name, ok, ok_message, fail_message, detail=None, fail_status=FAIL):
        """Add a pass/fail check in one call."""
        return self.add(
            name,
            PASS if ok else fail_status,
            ok_message if ok else fail_message,
            detail,
        )

    @property
    def status(self):
        worst = PASS if self.checks else INFO
        for check in self.checks:
            if _SEVERITY[check["status"]] > _SEVERITY[worst]:
                worst = check["status"]
        return worst

    def has(self, status):
        return any(check["status"] == status for check in self.checks)

    def to_dict(self):
        return {
            "tool": self.tool,
            "schema_version": 1,
            "status": self.status,
            "started_at": datetime.datetime.fromtimestamp(self.started_at).isoformat(timespec="seconds"),
            "duration_s": round(time.time() - self.started_at, 3),
            "context": self.context,
            "checks": self.checks,
        }

    def write_json(self, path):
        if not path:
            return
        directory = os.path.dirname(os.path.abspath(path))
        if directory:
            os.makedirs(directory, exist_ok=True)
        with open(path, "w") as handle:
            json.dump(self.to_dict(), handle, indent=2, sort_keys=False, default=str)
            handle.write("\n")

    def print_human(self, stream=None, verbose=False):
        stream = stream or sys.stdout
        width = max([len(c["name"]) for c in self.checks] + [4])
        stream.write("\n%s\n" % self.tool)
        for key, value in self.context.items():
            stream.write("  %-10s %s\n" % (key + ":", value))
        stream.write("\n")
        for check in self.checks:
            stream.write(
                "  [%s] %-*s  %s\n"
                % (_paint(check["status"], stream), width, check["name"], check["message"])
            )
            detail = check["detail"]
            if detail and (verbose or check["status"] in (FAIL, ERROR, WARN)):
                for line in _format_detail(detail):
                    stream.write("        %s\n" % line)
        stream.write(
            "\n  => %s (%.1fs)\n\n"
            % (_paint(self.status, stream), time.time() - self.started_at)
        )

    def exit_code(self, warn_is_failure=False):
        status = self.status
        if status == ERROR:
            return EXIT_ERROR
        if status == FAIL:
            return EXIT_FAIL
        if status == WARN and warn_is_failure:
            return EXIT_FAIL
        return EXIT_PASS

    def finish(self, json_path=None, verbose=False, warn_is_failure=False):
        self.write_json(json_path)
        self.print_human(verbose=verbose)
        return self.exit_code(warn_is_failure=warn_is_failure)


def _format_detail(detail, limit=20):
    lines = []
    if isinstance(detail, dict):
        items = list(detail.items())
        for key, value in items[:limit]:
            if isinstance(value, (list, tuple)):
                rendered = ", ".join(str(v) for v in value[:6])
                if len(value) > 6:
                    rendered += ", ... (%d total)" % len(value)
            elif isinstance(value, dict):
                rendered = json.dumps(value, default=str)
                if len(rendered) > 200:
                    rendered = rendered[:200] + " ..."
            else:
                rendered = str(value)
            lines.append("%s: %s" % (key, rendered))
        if len(items) > limit:
            lines.append("... (%d more fields)" % (len(items) - limit))
    else:
        lines.append(str(detail))
    return lines


def rel_diff(a, b):
    """Symmetric relative difference; 0 when both sides are 0."""
    if a == b:
        return 0.0
    denom = max(abs(a), abs(b))
    if denom == 0:
        return 0.0
    return abs(a - b) / denom


def load_config(defaults, path=None):
    """Overlay a JSON config file on top of ``defaults``."""
    config = dict(defaults)
    if path:
        with open(path) as handle:
            config.update(json.load(handle))
    return config


def add_common_args(parser):
    parser.add_argument("--json", metavar="PATH", help="write the machine-readable report here")
    parser.add_argument("--config", metavar="PATH", help="JSON file overriding the default thresholds")
    parser.add_argument("-v", "--verbose", action="store_true", help="print detail for passing checks too")
    parser.add_argument(
        "--warn-is-failure",
        action="store_true",
        help="exit non-zero on WARN (use for release gating)",
    )
    return parser
