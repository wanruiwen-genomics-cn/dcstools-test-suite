#!/usr/bin/env python3
"""Render the top-level test entry page: reports/index.html.

The per-run html_report.py renders one suite run. This script aggregates the
*latest* result of every case across all case reports under reports/<case>/,
so opening reports/index.html gives the current state of the whole suite at a
glance -- the same overview + case-cards + click-to-expand form as a single
run report, plus a "run history" section linking to each suite run.

It scans case.json files directly (not just suite.json), so a standalone
`cases/<x>.sh --compare-only` rerun is reflected in the index without needing a
full suite run.

    python3 scripts/index_report.py [--reports-dir reports] [-o reports/index.html]
"""
from __future__ import annotations
import argparse, glob, json, os, sys
from collections import Counter

# Reuse the per-run renderer so the index looks identical to a run report.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import html_report as hr

RANK = {"PASS": 0, "INFO": 0, "SKIP": 1, "WARN": 2, "FAIL": 3, "ERROR": 4}


def _load(p):
    if not p or not os.path.exists(p):
        return None
    try:
        with open(p) as f:
            return json.load(f)
    except Exception:
        return None


def discover_case_reports(reports_dir):
    """All case.json files under reports/<case>/<run_id>/, as
    (case_name, run_id, path, data). Newest run_id sorts last."""
    out = []
    for sp in glob.glob(os.path.join(reports_dir, "*", "*", "case.json")):
        # skip reports/suite/<run_id>/case.json (suite dir is not a case)
        parts = os.path.relpath(sp, reports_dir).split(os.sep)
        if len(parts) != 3 or parts[0] == "suite":
            continue
        case, run_id = parts[0], parts[1]
        data = _load(sp)
        if not data:
            continue
        out.append((case, run_id, sp, data))
    out.sort(key=lambda r: r[1])  # by run_id ascending; newest last
    return out


def build_item(case, run_id, sp, data, reports_dir):
    """Build a suite-style item dict from a case.json for html_report.render."""
    base = os.path.dirname(sp)
    checks = data.get("checks", [])
    counts = Counter(c.get("status", "INFO") for c in checks)
    problems = []
    for c in checks:
        st = c.get("status", "INFO")
        if st in ("FAIL", "ERROR", "WARN"):
            problems.append("%s: %s" % (c.get("name", ""), c.get("message", "")))
    item = {
        "name": case,
        "kind": "case",
        "run_id": run_id,
        "status": data.get("status", "INFO"),
        "wall_s": data.get("duration_s", 0),
        "counts": dict(counts),
        "problems": problems,
        "report": sp,
        "_detail": data,
        "_comparisons": [],
    }
    for f in sorted(glob.glob(os.path.join(base, "*.json"))):
        if os.path.basename(f) == "case.json":
            continue
        d = _load(f)
        if d and "checks" in d:
            item["_comparisons"].append((os.path.basename(f), d))
    return item


def latest_per_case(reports_dir):
    """Newest case.json per case name."""
    reps = discover_case_reports(reports_dir)
    latest = {}  # case -> item
    for case, run_id, sp, data in reps:  # ascending run_id, so last wins
        latest[case] = build_item(case, run_id, sp, data, reports_dir)
    items = list(latest.values())
    verdict = "PASS"
    counts = {"PASS": 0, "FAIL": 0, "ERROR": 0, "SKIP": 0, "WARN": 0, "INFO": 0}
    for it in items:
        st = it.get("status", "INFO")
        counts[st] = counts.get(st, 0) + 1
        if RANK.get(st, 0) > RANK.get(verdict, 0):
            verdict = st
    suite = {
        "tool": "suite",
        "schema_version": 1,
        "run_id": "latest",
        "tier": "aggregate (latest per case)",
        "verdict": verdict,
        "counts": counts,
        "items": items,
    }
    return suite


def discover_suite_runs(reports_dir):
    """[(run_id, suite_dict, suite_json_path)] newest-first, for run history."""
    runs = []
    for sp in glob.glob(os.path.join(reports_dir, "suite", "*", "suite.json")):
        s = _load(sp)
        if not s:
            continue
        rid = s.get("run_id") or os.path.basename(os.path.dirname(sp))
        runs.append((rid, s, sp))
    runs.sort(key=lambda r: r[0], reverse=True)
    return runs


# Index grouping. Species is a row; the index mode is a column. Cases not listed
# here fall into "其他" so a newly added case still shows up.
GROUPS = [
    {
        "id": "smoke",
        "title": "冒烟与工具契约",
        "blurb": "几分钟内跑完，确认工具本身能跑通",
        "layout": "cards",
        "cases": [
            "wgs_human_small", "seqarc_roundtrip",
            "aligner_shared_mem", "jointcall",
        ],
    },
    {
        "id": "wgs",
        "title": "胚系 WGS · 二倍体",
        "blurb": "同一行是同一套数据，三列只换对齐索引",
        "layout": "matrix",
        "columns": [("default", "默认索引"), ("ert", "ERT"), ("minibwa", "minibwa")],
        "rows": [
            ("人 30x", {
                "default": "wgs_human_30x",
                "ert": "wgs_human_30x_ert",
                "minibwa": "wgs_human_30x_minibwa",
            }),
            ("小鼠", {
                "default": "wgs_mouse",
                "ert": "wgs_mouse_ert",
                "minibwa": "wgs_mouse_minibwa",
            }),
            ("水稻", {
                "default": "wgs_rice",
                "ert": "wgs_rice_ert",
                "minibwa": "wgs_rice_minibwa",
            }),
        ],
    },
    {
        "id": "wes",
        "title": "外显子组",
        "blurb": "流程与 WGS 相同，调用限制在目标区域 BED 内",
        "layout": "cards",
        "cases": ["wes_human"],
    },
    {
        "id": "poly",
        "title": "多倍体",
        "blurb": "马铃薯四倍体，基因分型走 jointcaller",
        "layout": "matrix",
        "columns": [("default", "默认索引"), ("ert", "ERT"), ("minibwa", "minibwa")],
        "rows": [
            ("马铃薯", {
                "default": "wgs_potato",
                "ert": "wgs_potato_ert",
                "minibwa": "wgs_potato_minibwa",
            }),
        ],
    },
    {
        "id": "somatic",
        "title": "体细胞",
        "blurb": "肿瘤 / 正常配对，Mutect2 + FilterMutectCalls",
        "layout": "cards",
        "cases": ["mutect2_tn"],
    },
]

INDEX_CSS = """
.case-group{margin:0 0 28px}
.group-head{display:flex;justify-content:space-between;align-items:flex-end;gap:16px;margin:4px 2px 12px}
.group-head h2{font-size:18px;font-weight:700;line-height:1.3}
.group-blurb{font-size:13px;color:var(--muted);margin-top:2px}
.group-score{font-size:13px;font-weight:700;padding:4px 12px;border-radius:999px;white-space:nowrap}
.mx{width:100%;border-collapse:separate;border-spacing:8px;margin:0 0 8px}
.mx th{font-size:12px;color:var(--muted);font-weight:600;text-align:left;padding:0 4px}
.mx-rowlabel{font-weight:600;font-size:13px;width:72px;white-space:nowrap;padding-right:4px}
.mx-cell{display:flex;align-items:center;gap:8px;width:100%;text-align:left;background:var(--card);border:1px solid var(--border);border-radius:10px;padding:10px 12px;cursor:pointer;font:inherit;color:inherit}
.mx-cell:hover{border-color:var(--accent)}
.mx-pass{border-color:#bbf7d0}
.mx-fail{border-color:#fecaca;background:#fff7f7}
.mx-warn{border-color:#fde68a;background:#fffbeb}
.mx-missing{cursor:default;color:var(--muted);background:transparent}
.mx-icon{font-weight:700;width:16px;flex-shrink:0}
.mx-name{font-weight:600;font-size:13px}
.mx-meta{margin-left:auto;font-size:12px;color:var(--muted);white-space:nowrap}
.group-pass-details{margin-top:4px}
.group-pass-details>summary{cursor:pointer;font-size:13px;color:var(--muted);margin:4px 2px 8px}
"""

_BAD = ("FAIL", "ERROR", "WARN")


def _group_names(group):
    if group["layout"] == "matrix":
        names = []
        for _label, row in group["rows"]:
            for key, _col in group["columns"]:
                names.append(row[key])
        return names
    return list(group["cases"])


def _score(items):
    """(pass_count, total, worst_status) over cases that actually ran."""
    present = [it for it in items if it is not None]
    worst = "PASS"
    n_pass = 0
    for it in present:
        st = it.get("status", "INFO")
        if st == "PASS":
            n_pass += 1
        if RANK.get(st, 0) > RANK.get(worst, 0):
            worst = st
    if not present:
        worst = "SKIP"
    return n_pass, len(present), worst


def _score_html(n_pass, total, worst):
    if total == 0:
        label, bg, fg = "未跑", "#f1f5f9", "#64748b"
    elif worst in ("FAIL", "ERROR"):
        label, bg, fg = "%d/%d 通过" % (n_pass, total), "#fee2e2", "#dc2626"
    elif worst == "WARN":
        label, bg, fg = "%d/%d 通过" % (n_pass, total), "#fef3c7", "#d97706"
    else:
        label, bg, fg = "%d/%d 通过" % (n_pass, total), "#dcfce7", "#16a34a"
    return '<span class="group-score" style="background:%s;color:%s">%s</span>' % (bg, fg, label)


def _expand_bad(html, item):
    if item.get("status") not in _BAD:
        return html
    name = item["name"]
    html = html.replace(
        'id="detail-%s" style="display:none"' % name,
        'id="detail-%s" style="display:block"' % name, 1)
    html = html.replace(
        'id="expand-%s">▼' % name,
        'id="expand-%s" style="transform:rotate(180deg)">▼' % name, 1)
    return html


def _matrix_cell(item, expected_name):
    if item is None:
        return ('<div class="mx-cell mx-missing"><span class="mx-name">未跑</span>'
                '<span class="mx-meta">%s</span></div>' % hr.esc(expected_name))
    st = item.get("status", "INFO")
    cls = {"PASS": "mx-pass", "FAIL": "mx-fail", "ERROR": "mx-fail",
           "WARN": "mx-warn"}.get(st, "")
    return ('<button type="button" class="mx-cell %s" onclick="openCase(\'%s\')">'
            '<span class="mx-icon" style="color:%s">%s</span>'
            '<span class="mx-name">%s</span>'
            '<span class="mx-meta">%s</span></button>'
            % (cls, hr.esc(item["name"]), hr.CLR.get(st, "#64748b"),
               hr.ICON.get(st, "?"), hr.esc(item["name"]),
               hr.dur(item.get("wall_s"))))


def _render_matrix(group, by_name):
    cols = group["columns"]
    p = ['<table class="mx"><thead><tr><th></th>']
    for _key, label in cols:
        p.append('<th>%s</th>' % hr.esc(label))
    p.append('</tr></thead><tbody>')
    for row_label, row in group["rows"]:
        p.append('<tr><th class="mx-rowlabel">%s</th>' % hr.esc(row_label))
        for key, _label in cols:
            name = row[key]
            p.append('<td>%s</td>' % _matrix_cell(by_name.get(name), name))
        p.append('</tr>')
    p.append('</tbody></table>')
    return "\n".join(p)


def _render_cards(items, hide_pass):
    """Case cards. Bad ones stay visible and start expanded. Pass cards can be
    tucked under a disclosure so a 9-wide matrix is not followed by 9 headers."""
    bad, ok = [], []
    for it in items:
        html = _expand_bad(hr.render_case(it), it)
        (bad if it.get("status") in _BAD else ok).append(html)
    parts = list(bad)
    if ok and hide_pass:
        parts.append('<details class="group-pass-details"><summary>通过项详情（%d）</summary>'
                     % len(ok))
        parts.extend(ok)
        parts.append('</details>')
    else:
        parts.extend(ok)
    return "\n".join(parts)


def render_grouped(items):
    by_name = {it["name"]: it for it in items}
    claimed = set()
    sections = []
    for group in GROUPS:
        names = _group_names(group)
        claimed.update(names)
        present = [by_name.get(n) for n in names]
        ran = [it for it in present if it is not None]
        n_pass, total, worst = _score(ran)
        head = ('<div class="group-head"><div><h2>%s</h2>'
                '<div class="group-blurb">%s</div></div>%s</div>'
                % (hr.esc(group["title"]), hr.esc(group["blurb"]),
                   _score_html(n_pass, total, worst)))
        body = []
        if group["layout"] == "matrix":
            body.append(_render_matrix(group, by_name))
        body.append(_render_cards(ran, hide_pass=(group["layout"] == "matrix")))
        sections.append('<section class="case-group" id="group-%s">%s%s</section>'
                        % (group["id"], head, "\n".join(body)))

    rest = [it for it in items if it["name"] not in claimed]
    if rest:
        n_pass, total, worst = _score(rest)
        head = ('<div class="group-head"><div><h2>其他</h2>'
                '<div class="group-blurb">未归入上面几组的用例</div></div>%s</div>'
                % _score_html(n_pass, total, worst))
        sections.append('<section class="case-group" id="group-other">%s%s</section>'
                        % (head, _render_cards(rest, hide_pass=False)))

    opener = (
        '<script>function openCase(id){'
        'var d=document.getElementById("detail-"+id);'
        'var e=document.getElementById("expand-"+id);'
        'var card=document.getElementById("case-"+id);'
        'if(!d||!card)return;'
        'var box=card.closest("details");if(box)box.open=true;'
        'd.style.display="block";'
        'if(e)e.style.transform="rotate(180deg)";'
        'card.scrollIntoView({behavior:"smooth",block:"start"});'
        '}</script>'
    )
    return opener + "\n".join(sections)


def render_run_history(runs, reports_dir):
    p = ['<div class="subsection"><h4>运行历史（所有 suite 运行）</h4>']
    p.append('<table class="table table-compact"><thead><tr>'
              '<th>Run ID</th><th>分层</th><th>判定</th><th>耗时</th>'
              '<th>用例数</th><th>报告</th></tr></thead><tbody>')
    for rid, s, sp in runs:
        tier = s.get("tier", "")
        v = s.get("verdict", "?")
        items = s.get("items", [])
        wall = sum(i.get("wall_s", 0) for i in items)
        rel = os.path.relpath(os.path.dirname(sp), reports_dir)
        link = os.path.join(rel, "report.html")
        p.append('<tr><td><code>%s</code></td><td>%s</td><td>%s</td><td>%s</td>'
                  '<td>%d</td><td><a href="%s">打开 ↗</a></td></tr>'
                  % (hr.esc(rid), hr.esc(tier), hr.badge(v, ""),
                     hr.dur(wall), len(items), hr.esc(link)))
    p.append('</tbody></table></div>')
    return "\n".join(p)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reports-dir", default="reports",
                        help="suite reports root (default: reports)")
    parser.add_argument("-o", "--out", help="output html (default: <reports-dir>/index.html)")
    args = parser.parse_args(argv)

    reports_dir = os.path.abspath(args.reports_dir)
    suite = latest_per_case(reports_dir)
    runs = discover_suite_runs(reports_dir)

    overview = hr.render_overview(suite)
    cases = render_grouped(suite["items"])
    # Run history is not shown. Append render_run_history(runs, reports_dir)
    # to {cases} to list reports/suite/<run_id>/ again.

    out = args.out or os.path.join(reports_dir, "index.html")
    content = (hr.HTML_TEMPLATE
               .replace("{run_id}", hr.esc("latest aggregate"))
               .replace("{css}", hr.CSS + INDEX_CSS)
               .replace("{overview}", overview)
               .replace("{cases}", cases))
    with open(out, "w") as f:
        f.write(content)
    print("index report: %s (%d cases, %d suite runs)" %
          (out, len(suite["items"]), len(runs)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
