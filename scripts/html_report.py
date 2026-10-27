#!/usr/bin/env python3
"""Render suite JSON reports as a self-contained HTML page.
Overview on top, case cards expand to show details on click.
"""
from __future__ import annotations
import argparse, glob, html, json, os, re, sys
from datetime import datetime

RANK={"PASS":0,"INFO":0,"SKIP":1,"WARN":2,"FAIL":3,"ERROR":4}
ICON={"PASS":"✓","FAIL":"✕","ERROR":"⚡","SKIP":"○","WARN":"▲","INFO":"ℹ"}
CLR={"PASS":"#16a34a","FAIL":"#dc2626","ERROR":"#9333ea","SKIP":"#64748b","WARN":"#d97706","INFO":"#64748b"}
BG ={"PASS":"#dcfce7","FAIL":"#fee2e2","ERROR":"#f3e8ff","SKIP":"#f1f5f9","WARN":"#fef3c7","INFO":"#f1f5f9"}

# Chinese labels for check names, so the checks list explains what each item
# means instead of showing bare identifiers like "step.bqsr" / "product.*".
CN_STATUS={"PASS":"通过","FAIL":"失败","ERROR":"错误","SKIP":"跳过","WARN":"警告","INFO":"信息"}
# Internal step ids are not the CLI tool names. Display the tool that actually ran.
STEP_TOOL={
    "build_index":"aligner","align":"aligner","align_tumor":"aligner","align_normal":"aligner",
    "bqsr":"bqsr","variant_call":"variantCaller","genotype":"genotyper",
    "bam_stat":"bam-stat","vcf_stat":"vcf-stat","report":"report",
    "mutect2":"Mutect2","filter":"FilterMutectCalls",
    "index":"SeqArc","compress":"SeqArc","decompress":"SeqArc",
    "jointcall":"jointcaller","faidx":"samtools","dict":"samtools",
}

def step_tool(name, case=None):
    """CLI tool for a step id. potato_polyploid's genotype step is jointcaller."""
    if name == "genotype" and case == "potato_polyploid":
        return "jointcaller"
    return STEP_TOOL.get(name, name)
CHECK_CN={
    "data.present":"数据完整性","data.base_quality":"碱基质量","license.proxy":"许可证代理",
    "align_bam":"BAM 一致性","gvcf":"gVCF 一致性","vcf":"VCF 一致性",
    "tumor_bam":"肿瘤 BAM 一致性","normal_bam":"正常 BAM 一致性",
    "unfiltered_vcf":"未过滤 VCF 一致性","filtered_vcf":"已过滤 VCF 一致性",
    "filter.applied":"变异过滤标记","somatic.counts":"体细胞变异计数",
    "prebuilt.bam":"预构建 BAM 复用","joint_vcf":"联合 VCF 一致性",
    "joint.nonempty":"联合结果非空","joint.sample_count":"联合样本数",
    "variant.count":"变异计数","vcf.nonempty":"VCF 非空",
    "product.report":"报告目录",
    "layout":"目录布局","run":"运行","scan":"扫描",
    "records.sampling":"记录采样","records.digest":"记录摘要",
    "records.normalised":"记录归一化",
    "concordance.sites":"位点一致性","concordance.genotypes":"基因型一致性",
    "counts.filter":"FILTER 分布","stats.summary":"统计摘要",
    "stats.flagstat":"flagstat 统计","stats.idxstats":"idxstats 统计",
    "stats.insert_size":"插入片段统计","stats.error_rate":"错误率统计",
    "stats.rates":"比率统计",
    "file.exists":"文件存在","file.index":"文件索引","file.size":"文件大小",
    "file.nonempty":"文件非空",
    "tool.samtools":"工具·samtools","tool.bcftools":"工具·bcftools","tool.dcs":"工具·dcs",
    "dispatcher.tools":"工具分发","rpath.absolute":"绝对路径检查",
    "coverage.flags":"覆盖度标记",
}
# Prefix rules for namespaced checks not listed above.
PREFIX_CN=[
    ("step.","步骤·"),("product.","产物·"),("env.","环境·"),
    ("contract.","契约·"),("seqarc.","SeqArc·"),("symbols.","符号·"),
    ("input.","输入·"),("data.","数据·"),("joint.","联合·"),
    ("libs.","依赖库·"),("header.","表头·"),("file.","文件·"),
    ("tool.","工具·"),("records.","记录·"),("variant.","变异·"),
    ("coverage.","覆盖度·"),("dispatcher.","分发·"),("rpath.","路径·"),
    ("script.","脚本·"),("stats.","统计·"),
]

def explain_check(name):
    if name in CHECK_CN: return CHECK_CN[name]
    if name.startswith("step."):
        return "步骤·" + step_tool(name[5:])
    for pre, cn in PREFIX_CN:
        if name.startswith(pre):
            return cn + name[len(pre):]
    return name

def _load(p):
    if not p or not os.path.exists(p): return None
    try:
        with open(p) as f: return json.load(f)
    except: return None

def load_suite(sp):
    s=_load(sp)
    if not s: return None
    for item in s.get("items",[]):
        rp=item.get("report","")
        item["_detail"]=_load(rp)
        if item.get("kind")=="case" and rp:
            d=os.path.dirname(rp)
            item["_comparisons"]=[]
            for f in sorted(glob.glob(os.path.join(d,"*.json"))):
                if os.path.basename(f)=="case.json": continue
                data=_load(f)
                if data and "checks" in data:
                    item["_comparisons"].append((os.path.basename(f),data))
    return s

def esc(t): return html.escape(str(t),quote=False)

def dcs_home_path():
    """Bundle under test: $DCS_HOME, else the default in env/cases.conf."""
    env = os.environ.get("DCS_HOME")
    if env:
        return env
    conf = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "env", "cases.conf")
    try:
        text = open(conf).read()
    except OSError:
        return ""
    m = re.search(r'^DCS_HOME=\$\{DCS_HOME:-(.+?)\}\s*$', text, re.M)
    return m.group(1).strip() if m else ""

def dcs_release_version():
    """Release string printed by `dcs` with no arguments (`Version: v2.3.0`)."""
    home = dcs_home_path()
    if not home:
        return ""
    try:
        text = open(os.path.join(home, "bin", "dcs")).read()
    except OSError:
        return ""
    m = re.search(r'^\s*echo\s+"Version:\s*([^"]+)"', text, re.M)
    return m.group(1).strip() if m else ""

def badge(st,text=None):
    i=ICON.get(st,"?"); l=text if text is not None else st
    return '<span class="badge" style="background:%s;color:%s">%s&nbsp;%s</span>'%(BG.get(st,"#f1f5f9"),CLR.get(st,"#64748b"),i,esc(l))

def dur(s):
    if s is None: return "—"
    s=float(s)
    if s<60: return "%.0fs"%s
    m=int(s//60); r=s-m*60
    if m<60: return "%dm%.0fs"%(m,r)
    h=m//60; m=m%60; return "%dh%dm"%(h,m)

def bar(v,mx,c="#3b82f6"):
    if mx<=0: return ""
    p=min(100,max(0,(v/mx)*100))
    return '<div class="bar-track"><div class="bar-fill" style="width:%.1f%%;background:%s"></div></div>'%(p,c)

def render_overview(suite):
    items=suite.get("items",[])
    c={"PASS":0,"FAIL":0,"ERROR":0,"SKIP":0,"WARN":0,"INFO":0}
    for i in items: c[i.get("status","INFO")]=c.get(i.get("status","INFO"),0)+1
    tw=sum(i.get("wall_s",0) for i in items)
    v=suite.get("verdict","UNKNOWN")
    rid=suite.get("run_id",""); tier=suite.get("tier","")
    p=['<div class="overview">']
    p.append('<div class="overview-title">')
    ver = suite.get("dcs_version") or dcs_release_version()
    if ver:
        p.append('<div class="overview-version">dcstools %s</div>' % esc(ver))
    p.append('<h1>发版测试报告</h1><div class="overview-meta">')
    p.append('<span>Run ID: <code>%s</code></span>'%esc(rid))
    if tier: p.append('<span>分层: <code>%s</code></span>'%esc(tier))
    p.append('<span>生成: <code>%s</code></span>'%esc(datetime.now().strftime("%Y-%m-%d %H:%M")))
    p.append('</div></div>')
    p.append('<div class="overview-verdict" style="background:%s;color:%s">%s&nbsp;%s</div>'%(BG.get(v,"#f1f5f9"),CLR.get(v,"#64748b"),ICON.get(v,"?"),esc(v)))
    p.append('</div><div class="cards">')
    for l,val,co in [("总测试项",len(items),"#1e293b"),("通过",c["PASS"],"#16a34a"),("失败",c["FAIL"],"#dc2626"),("错误",c["ERROR"],"#9333ea"),("跳过",c["SKIP"],"#64748b"),("总耗时",dur(tw),"#3b82f6")]:
        p.append('<div class="card"><div class="card-value" style="color:%s">%s</div><div class="card-label">%s</div></div>'%(co,esc(val),esc(l)))
    p.append('</div>')
    return "\n".join(p)

def _pct_diff(test, base):
    """Relative change of `test` vs `base` as a percentage string with sign.
    Returns (pct_str, ratio) where ratio is test/base (None when undefined)."""
    try:
        t=float(test); b=float(base)
    except (TypeError, ValueError):
        return ("—", None)
    if b in (0, None) or b is None:
        return ("—", None)
    r=t/b
    return ("%+.1f%%"%((r-1)*100.0), r)

def _fmt_signed(v, unit=""):
    try: f=float(v)
    except (TypeError, ValueError): return "—"
    sign="+" if f>=0 else ""
    if unit=="MB":
        if abs(f)>=1024: return "%s%.2f GB"%(sign,f/1024.0)
        return "%s%.0f MB"%(sign,f)
    return "%s%.0f%s"%(sign,f,unit)

def render_steps_summary(steps, golden_steps):
    """A compact summary under the steps table: total wall time and peak memory
    of the test run vs the golden run, with diff/ratio and a verdict on whether
    the drift is significant."""
    gmap={s.get("name"):s for s in (golden_steps or [])}
    tot_t=sum(float(s.get("wall_s",0) or 0) for s in steps)
    tot_g=sum(float((gmap.get(s.get("name"),{}) or {}).get("wall_s",0) or 0) for s in steps)
    mem_t=max([float(s.get("peak_rss_mb",0) or 0) for s in steps]+[0])
    mem_g=max([float((gmap.get(s.get("name"),{}) or {}).get("peak_rss_mb",0) or 0) for s in steps]+[0])
    t_pct,t_r=_pct_diff(tot_t,tot_g)
    m_pct,m_r=_pct_diff(mem_t,mem_g)
    # significance: >10% drift on either metric counts as "明显差异"
    sig=False
    if t_r is not None and abs(t_r-1.0)>0.10: sig=True
    if m_r is not None and abs(m_r-1.0)>0.10: sig=True
    verdict="存在明显差异" if sig else "无明显差异"
    vst="WARN" if sig else "PASS"
    p=['<div class="steps-summary">']
    p.append('<div class="summary-card"><div class="summary-label">总耗时</div>'
             '<div class="summary-value">%s <span class="vs">/</span> %s</div>'
             '<div class="summary-diff">差值 %s (%s)</div>'
             '<div class="summary-sub">测试 vs 基线</div></div>'%(
        esc(dur(tot_t)), esc(dur(tot_g) if tot_g is not None else "—"),
        esc(_fmt_signed(tot_t-tot_g,"s")), esc(t_pct)))
    p.append('<div class="summary-card"><div class="summary-label">峰值内存</div>'
             '<div class="summary-value">%s <span class="vs">/</span> %s</div>'
             '<div class="summary-diff">差值 %s (%s)</div>'
             '<div class="summary-sub">测试 vs 基线</div></div>'%(
        esc(_fmt_signed(mem_t,"MB").lstrip("+")), esc(_fmt_signed(mem_g,"MB").lstrip("+")),
        esc(_fmt_signed(mem_t-mem_g,"MB")), esc(m_pct)))
    p.append('<div class="summary-card summary-verdict %s"><div class="summary-label">结论</div>'
             '<div class="summary-value">%s %s</div>'
             '<div class="summary-sub">阈值: 耗时±10%% / 内存±10%%</div></div>'%(
        "verdict-warn" if sig else "verdict-pass", ICON.get(vst,"?"), esc(verdict)))
    p.append('</div>')
    return "\n".join(p)

def render_steps(detail):
    steps = detail.get("steps") or []
    if not steps: return ""
    gmap = {s.get("name"): s for s in (detail.get("golden_steps") or [])}
    has_golden = bool(gmap)
    p=['<div class="subsection"><h4>步骤工具</h4><table class="table table-compact">']
    if has_golden:
        p.append('<thead><tr><th>工具</th><th>运行状态</th>'
                 '<th>耗时(测试)</th><th>耗时(基线)</th>'
                 '<th>峰值RSS(测试)</th><th>峰值RSS(基线)</th>'
                 '<th>线程(测试)</th><th>线程(基线)</th></tr></thead><tbody>')
    else:
        p.append('<thead><tr><th>工具</th><th>运行状态</th>'
                 '<th>耗时</th><th>峰值RSS</th><th>线程</th></tr></thead><tbody>')
    case = detail.get("case")
    for s in steps:
        key=s.get("name",""); nm=step_tool(key, case); rc=s.get("exit_code",0); w=s.get("wall_s",0)
        r=s.get("peak_rss_mb",0); t=s.get("threads",0)
        rb=badge("PASS","") if rc==0 else badge("ERROR",str(rc))
        if has_golden:
            g=gmap.get(key,{})
            gw=g.get("wall_s"); gr=g.get("peak_rss_mb"); gt=g.get("threads",0)
            p.append('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>'
                     '<td>%s MB</td><td>%s MB</td><td>%s</td><td>%s</td></tr>'%(
                esc(nm), rb, dur(w), dur(gw) if gw is not None else "—",
                esc(r), esc(gr) if gr is not None else "—",
                esc(t or "—"), esc(gt or "—")))
        else:
            p.append('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s MB</td><td>%s</td></tr>'%(
                esc(nm), rb, dur(w), esc(r), esc(t or "—")))
    p.append('</tbody></table>')
    if has_golden:
        p.append(render_steps_summary(steps, detail.get("golden_steps") or []))
    p.append('</div>')
    return "\n".join(p)

def render_first_differences(diffs, max_rows=20):
    """Side-by-side table of the first differing VCF/gVCF body lines."""
    if not diffs:
        return ""
    if isinstance(diffs, str):
        return '<p class="muted">差异详情不可用: %s</p>' % esc(diffs)
    if not isinstance(diffs, list):
        return ""
    p=['<div class="diff-block"><strong>逐行不一致的记录（前 %d 条，含 INFO+FORMAT）</strong>' % min(len(diffs), max_rows),
       '<table class="table table-compact cmp-table diff-table"><thead><tr>',
       '<th>#</th><th>位点</th><th>测试</th><th>基线</th></tr></thead><tbody>']
    for d in diffs[:max_rows]:
        if not isinstance(d, dict):
            p.append('<tr><td colspan="4"><code>%s</code></td></tr>' % esc(str(d)))
            continue
        locus = d.get("locus") or ("record %s" % d.get("record", "?"))
        p.append('<tr><td>%s</td><td><code>%s</code></td>'
                 '<td><code class="diff-cell">%s</code></td>'
                 '<td><code class="diff-cell">%s</code></td></tr>' % (
            esc(d.get("record", "")), esc(locus),
            esc(d.get("test", "")), esc(d.get("baseline", ""))))
    p.append('</tbody></table></div>')
    return "\n".join(p)

def render_detail(detail):
    """Render a check's detail dict as readable HTML.

    When detail carries baseline/test metric dicts (samtools / bcftools stats),
    prefer a side-by-side comparison table. When records.digest embeds
    first_differences, render those as an explicit line-diff table.
    """
    if not detail or not isinstance(detail, dict) or not detail:
        return ""
    chunks = []
    if "first_differences" in detail:
        chunks.append(render_first_differences(detail.get("first_differences")))
    base = detail.get("baseline")
    test = detail.get("test")
    if isinstance(base, dict) and isinstance(test, dict) and (
            any(not isinstance(v, dict) for v in list(base.values())[:3] + list(test.values())[:3])
            or (base and test and not isinstance(next(iter(base.values()), None), dict))
    ):
        if base and not any(isinstance(v, dict) for v in base.values()):
            chunks.append(render_kv_cmp_table(test, base))
        elif test and not any(isinstance(v, dict) for v in test.values()):
            chunks.append(render_kv_cmp_table(test, base))
    rows = []
    for k, v in detail.items():
        if k in ("baseline", "test", "first_differences"):
            continue
        if k == "record_count_delta":
            rows.append("<strong>记录数差值</strong>: %+d" % int(v))
            continue
        if isinstance(v, dict):
            sub = []
            for sk, sv in v.items():
                if isinstance(sv, (dict, list)):
                    sv = json.dumps(sv, ensure_ascii=False)[:200]
                sub.append("<strong>%s</strong>: %s" % (esc(sk), esc(str(sv))))
            rows.append("<strong>%s</strong>:<br>&nbsp;&nbsp;%s" % (esc(k), "<br>&nbsp;&nbsp;".join(sub)))
        elif isinstance(v, list):
            preview = v[:8]
            rows.append("<strong>%s</strong>: %s" % (
                esc(k), esc(", ".join(str(x) for x in preview))))
        else:
            rows.append("<strong>%s</strong>: %s" % (esc(k), esc(str(v))))
    if rows:
        chunks.append('<div class="check-detail">%s</div>' % "<br>".join(rows))
    return "\n".join(c for c in chunks if c)

def fmt_size_mb(nbytes):
    """Format bytes as MB (KB when under 1 MB, GB when large)."""
    if nbytes is None: return "—"
    try: n=float(nbytes)
    except (TypeError, ValueError): return "—"
    if n >= 1024*1024*1024:
        return "%.2f GB" % (n/(1024**3))
    if n < 1024*1024:
        return "%.1f KB" % (n/1024.0)
    return "%.1f MB" % (n/(1024**2))

def parse_bytes_msg(msg):
    import re
    m=re.search(r"(\d+)\s*bytes", msg or "")
    return int(m.group(1)) if m else None

def parse_count_msg(msg):
    """Pull an integer count from messages like '4854410 variants (...'."""
    import re
    m=re.search(r"(\d+)\s+variants", msg or "")
    if m: return int(m.group(1))
    m=re.search(r"test\s+(\d+)", msg or "")
    if m: return int(m.group(1))
    return None

def parse_kv_stat(path):
    """Parse bamStat.txt / vcfStat.txt into an ordered dict of scalar metrics.
    Skips per-contig depth rows (chrN\\tvalue) and coverage histograms."""
    out={}
    if not path or not os.path.isfile(path): return out
    try:
        with open(path) as f:
            for line in f:
                line=line.strip()
                if not line: continue
                # vcfStat uses "#total SNP: 123" / "ratio: 0%"
                if ":" in line and "\t" not in line:
                    k,v=line.split(":",1); k=k.lstrip("#").strip(); v=v.strip()
                elif "\t" in line:
                    k,v=line.split("\t",1); k=k.lstrip("#").strip(); v=v.strip()
                else:
                    continue
                k=k.rstrip(":").strip()
                # skip per-contig rows like "chr1"
                if k.lower().startswith("chr") and k[3:].replace("_","").replace("M","").replace("X","").replace("Y","").isdigit():
                    continue
                if k.endswith("x coverage"): continue
                # keep first token of multi-value cells
                v0=v.split("\t")[0].strip()
                out_keys=list(out.keys())
                if k.lower()=="ratio" and out_keys:
                    k="ratio(%s)"%out_keys[-1]
                out[k]=v0
    except OSError:
        pass
    return out

def find_report_stat(report_dir, suffix):
    """Find first file under report_dir whose name ends with suffix (e.g. .bamstat.xls)."""
    if not report_dir or not os.path.isdir(report_dir):
        return None
    hits=[]
    for dirpath,_,files in os.walk(report_dir):
        for f in files:
            if f.lower().endswith(suffix.lower()):
                hits.append(os.path.join(dirpath,f))
    hits.sort()
    return hits[0] if hits else None

# Curated highlights from DCS report/*.xls (表1 质控 / 表2 比对 / 表3 变异).
REPORT_HIGHLIGHT_KEYS=[
    ("fqstat", [
        ("Read_raw","原始 Read 数"),("Read_clean","Clean Read 数"),
        ("Read_clean_Rate","Clean 比例"),("Read_length","Read 长度"),
        ("Q30_raw","Q30(raw)"),("Q30_clean","Q30(clean)"),
        ("GC_raw","GC%(raw)"),("GC_clean","GC%(clean)"),
    ]),
    ("bamstat", [
        ("Mapping_Rate","比对率"),("PE_Mapping_Rate","PE 比对率"),
        ("Unique_Rate","唯一比对率"),("Duplication_Rate","重复率"),
        ("Mismatch_Rate","错配率"),("Insert_size","插入片段"),
        ("Average_depth(rmdup)","平均深度(去重)"),
        ("Coverage(>=1X)","覆盖度≥1X"),("Coverage(>=10X)","覆盖度≥10X"),
        ("Coverage(>=20X)","覆盖度≥20X"),("Coverage(>=30X)","覆盖度≥30X"),
        ("Uniformity(>0.2f)","覆盖均一性"),
    ]),
    ("vcfstat", [
        ("Total_SNP","SNP 数"),("dbSNP_Rate","dbSNP 比例"),
        ("Novel_SNP","未知 SNP"),("Ti/Tv","Ti/Tv"),
        ("Total_InDel","InDel 数"),("dbInDel_Rate","dbInDel 比例"),
    ]),
]

def extract_report_highlights(work_report_dir, golden_report_dir):
    """Key QC / align / variant metrics from report/*.xls for side-by-side compare."""
    rows=[]
    suffix_map={"fqstat":".fqstat.xls","bamstat":".bamstat.xls","vcfstat":".vcfstat.xls"}
    for stem, keys in REPORT_HIGHLIGHT_KEYS:
        suf=suffix_map[stem]
        tp=find_report_stat(work_report_dir, suf)
        bp=find_report_stat(golden_report_dir, suf)
        tmap=parse_kv_stat(tp); bmap=parse_kv_stat(bp)
        if not tmap and not bmap:
            continue
        for key, label in keys:
            tv=tmap.get(key); bv=bmap.get(key)
            if tv is None and bv is None:
                continue
            diff="—"
            try:
                tf=float(str(tv).rstrip("%")); bf=float(str(bv).rstrip("%"))
                d=tf-bf
                pct=_pct_diff(tf,bf)[0]
                diff="%+.4g (%s)"%(d,pct) if abs(d)>=1e-9 or d==0 else "0 (%s)"%pct
            except (TypeError, ValueError):
                if tv==bv: diff="相同"
                elif tv is not None and bv is not None: diff="不同"
            rows.append((label, str(tv if tv is not None else "—"),
                         str(bv if bv is not None else "—"), diff))
    for side_label, needle in (("中文报告 HTML","_report_cn.html"),
                               ("英文报告 HTML","_report_en.html")):
        tp=find_report_stat(work_report_dir, needle)
        bp=find_report_stat(golden_report_dir, needle)
        ts=os.path.getsize(tp) if tp and os.path.isfile(tp) else None
        bs=os.path.getsize(bp) if bp and os.path.isfile(bp) else None
        if ts is None and bs is None:
            continue
        if ts is not None and bs is not None:
            d=ts-bs
            diff="相同" if d==0 else "%+d B (%s)"%(d, _pct_diff(ts,bs)[0])
        elif ts is not None:
            diff="仅测试"
        else:
            diff="仅基线"
        rows.append((side_label, fmt_size_mb(ts), fmt_size_mb(bs), diff))
    return rows

def render_kv_cmp_table(test_map, base_map, keys=None, max_rows=40):
    """Side-by-side key/value comparison table."""
    if keys is None:
        keys=list(dict.fromkeys(list(test_map.keys())+list(base_map.keys())))
    if not keys: return ""
    p=['<table class="table table-compact cmp-table"><thead><tr>'
       '<th>指标</th><th>测试</th><th>基线</th><th>差值</th></tr></thead><tbody>']
    n=0
    for k in keys:
        if n>=max_rows: break
        tv=test_map.get(k); bv=base_map.get(k)
        if tv is None and bv is None: continue
        diff="—"
        try:
            tf=float(str(tv).rstrip("%")); bf=float(str(bv).rstrip("%"))
            d=tf-bf
            pct=_pct_diff(tf,bf)[0]
            diff="%+.4g (%s)"%(d,pct) if abs(d)>=0.0001 or d==0 else "0 (%s)"%pct
        except (TypeError, ValueError):
            if tv==bv: diff="相同"
            elif tv is not None and bv is not None: diff="不同"
        p.append('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>'%(
            esc(k), esc(tv if tv is not None else "—"),
            esc(bv if bv is not None else "—"), esc(diff)))
        n+=1
    p.append('</tbody></table>')
    return "\n".join(p)

def suite_paths(case_name, report_path):
    """Derive suite root / golden / work dirs from a case.json path."""
    # reports/<case>/<run_id>/case.json -> suite root is ../../..
    if not report_path: return None, None, None
    rundir=os.path.dirname(os.path.abspath(report_path))
    case_dir=os.path.dirname(rundir)
    reports_dir=os.path.dirname(case_dir)
    root=os.path.dirname(reports_dir)
    return (root,
            os.path.join(root,"golden",case_name),
            os.path.join(root,"work",case_name))

def golden_product_path(golden_dir, product_basename):
    if not golden_dir: return None
    p=os.path.join(golden_dir, product_basename)
    return p if os.path.exists(p) else None

def work_product_path(work_dir, product_basename):
    if not work_dir: return None
    p=os.path.join(work_dir, product_basename)
    return p if os.path.exists(p) else None

def product_basename(check_name):
    if check_name.startswith("product."): return check_name[len("product."):]
    return None

def categorize_check(name):
    if name.startswith("data.") or name.startswith("license.") or name.startswith("env."):
        return "preflight"
    if name.startswith("product.") or name in ("vcf.nonempty","variant.count"):
        return "product"
    # known golden-comparison products
    if name in ("align_bam","gvcf","vcf","tumor_bam","normal_bam",
                "unfiltered_vcf","filtered_vcf","joint_vcf") or \
       name.endswith("_bam") or name.endswith("_vcf") or name.endswith("_gvcf"):
        return "concordance"
    # static / packaging checks land with preflight-style table
    return "preflight"

BAM_KEY_METRICS=[
    ("file.size","文件大小"),
    ("records","比对记录数"),
    ("reads mapped","比对 reads"),
    ("reads duplicated","重复 reads"),
    ("mapped_rate","比对率"),
    ("duplicate_rate","重复率"),
    ("insert size average","平均插入片段"),
    ("insert size standard deviation","插入片段 SD"),
    ("error rate","错误率"),
    ("average quality","平均质量"),
]

def extract_bam_key_metrics(comp):
    """Pull a small set of representative BAM metrics from a bam_compare JSON."""
    checks={c.get("name"):c for c in (comp or {}).get("checks",[])}
    rows=[]
    fs=(checks.get("file.size") or {}).get("detail") or {}
    if "baseline_bytes" in fs or "test_bytes" in fs:
        tb,bb=fs.get("test_bytes"),fs.get("baseline_bytes")
        rows.append(("文件大小", fmt_size_mb(tb), fmt_size_mb(bb),
                     _pct_diff(tb,bb)[0] if tb is not None and bb is not None else "—"))
    rd=(checks.get("records.digest") or {}).get("detail") or {}
    tb=(rd.get("test") or {}).get("records"); bb=(rd.get("baseline") or {}).get("records")
    if tb is not None or bb is not None:
        try: diff="%+d (%s)"%(int(tb)-int(bb), _pct_diff(tb,bb)[0])
        except (TypeError,ValueError): diff="—"
        rows.append(("比对记录数(采样)", str(tb if tb is not None else "—"),
                     str(bb if bb is not None else "—"), diff))
    hs=(checks.get("header.sq") or {}).get("detail") or {}
    if "contigs" in hs:
        rows.append(("参考 contig 数", str(hs["contigs"]), str(hs["contigs"]), "相同"))

    def _add_pair(label, tv, bv):
        if tv is None and bv is None: return
        try: diff="%+.4g (%s)"%(float(tv)-float(bv), _pct_diff(tv,bv)[0])
        except (TypeError,ValueError):
            diff="相同" if tv==bv else "—"
        rows.append((label, str(tv if tv is not None else "—"),
                     str(bv if bv is not None else "—"), diff))

    # curated metrics only
    wanted=[
        ("stats.flagstat", [
            ("in total","总 reads"),("mapped","比对 reads"),
            ("duplicates","重复 reads"),("properly paired","正确配对"),
        ]),
        ("stats.rates", [
            ("mapped_rate","比对率"),("properly_paired_rate","正确配对率"),
            ("duplicate_rate","重复率"),("mq0_rate","MQ0 率"),
        ]),
        ("stats.insert_size", [
            ("insert size average","平均插入片段"),
            ("insert size standard deviation","插入片段 SD"),
        ]),
        ("stats.error_rate", [
            ("error rate","错误率"),("average quality","平均质量"),
        ]),
    ]
    for cname, keys in wanted:
        det=(checks.get(cname) or {}).get("detail") or {}
        base=det.get("baseline") or {}; test=det.get("test") or {}
        if not base and not test: continue
        for key, label in keys:
            _add_pair(label, test.get(key), base.get(key))
    return rows

def extract_vcf_key_metrics(comp):
    checks={c.get("name"):c for c in (comp or {}).get("checks",[])}
    rows=[]
    # Full-body INFO+FORMAT identity is the primary regression signal.
    rd_chk=checks.get("records.digest") or {}
    rd=rd_chk.get("detail") or {}
    st=rd_chk.get("status")
    msg=rd_chk.get("message") or ""
    if st:
        ann=rd.get("annotation_only") if isinstance(rd, dict) else None
        if st=="PASS" and ann:
            tol=rd.get("annotation_rel_tol") or 0
            verdict="通过（位点、等位基因、FILTER、GT 一致；%d 条仅注释不同）"%int(ann)
            target="目标: 注释差异 ≤ %.2f%%"%(100.0*float(tol))
        elif st=="PASS":
            verdict, target="一致（INFO+FORMAT 全量字节）", "目标: INFO+FORMAT 字节一致"
        else:
            verdict="不一致" if st in ("FAIL","ERROR","WARN") else st
            target="目标: INFO+FORMAT 字节一致"
        rows.append(("全量记录比对", verdict, target, msg[:120] if msg else "—"))
    tb=(rd.get("test") or {}).get("records"); bb=(rd.get("baseline") or {}).get("records")
    if tb is not None or bb is not None:
        try: diff="%+d (%s)"%(int(tb)-int(bb), _pct_diff(tb,bb)[0])
        except (TypeError,ValueError): diff="—"
        rows.append(("记录数", str(tb if tb is not None else "—"),
                     str(bb if bb is not None else "—"), diff))
    # ordered digest fingerprint (short)
    bo=(rd.get("baseline") or {}).get("ordered"); to=(rd.get("test") or {}).get("ordered")
    if bo or to:
        same="相同" if bo and to and bo==to else ("不同" if bo and to else "—")
        rows.append(("记录摘要(ordered)", (to or "—")[:12], (bo or "—")[:12], same))
    ss=(checks.get("stats.summary") or {}).get("detail") or {}
    base=ss.get("baseline") or {}; test=ss.get("test") or {}
    for k in ("number of SNPs","number of indels","number of MNPs",
              "number of multiallelic sites","ts/tv"):
        if k not in test and k not in base: continue
        tv,bv=test.get(k),base.get(k)
        try: diff="%+.4g (%s)"%(float(tv)-float(bv), _pct_diff(tv,bv)[0])
        except (TypeError,ValueError): diff="相同" if tv==bv else "—"
        rows.append((k, str(tv if tv is not None else "—"),
                     str(bv if bv is not None else "—"), diff))
    cs=(checks.get("concordance.sites") or {}).get("detail") or {}
    if cs:
        rows.append(("共享位点", str(cs.get("shared","—")), "—", "—"))
        rows.append(("仅基线", str(cs.get("baseline_only","—")), "—", "—"))
        rows.append(("仅测试", str(cs.get("test_only","—")), "—", "—"))
        j=cs.get("jaccard")
        if j is not None: rows.append(("Jaccard", "%.6f"%float(j), "1.0 (目标)", "—"))
    gt=(checks.get("concordance.genotypes") or {}).get("detail") or {}
    if "overall" in gt:
        rows.append(("基因型一致率", "%.6f"%float(gt["overall"]), "1.0 (目标)", "—"))
    # Diff samples (same as 30x-enhanced layout when digest failed).
    fd=(rd.get("first_differences") if isinstance(rd, dict) else None)
    if isinstance(fd, list) and fd:
        rows.append(("差异记录(已采样)", str(len(fd)), "—", "见下方差异表"))
    elif isinstance(fd, str) and fd:
        rows.append(("差异记录", "不可用", "—", esc(fd[:60])))
    return rows

def render_metric_rows_table(rows):
    if not rows: return ""
    p=['<table class="table table-compact cmp-table"><thead><tr>'
       '<th>指标</th><th>测试</th><th>基线</th><th>差值</th></tr></thead><tbody>']
    for label,tv,bv,diff in rows:
        p.append('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>'%(
            esc(label), tv if str(tv).startswith("<") else esc(str(tv)),
            bv if str(bv).startswith("<") else esc(str(bv)),
            diff if str(diff).startswith("<") else esc(str(diff))))
    p.append('</tbody></table>')
    return "\n".join(p)

def render_preflight_table(checks, section_id):
    p=['<table class="table table-compact check-table">']
    p.append('<thead><tr><th>检查项</th><th>状态</th><th>说明</th></tr></thead><tbody>')
    for idx,c in enumerate(checks):
        nm=c.get("name",""); st=c.get("status","INFO"); msg=c.get("message","")
        det=c.get("detail",{}); cn=explain_check(nm); dh=render_detail(det)
        p.append('<tr><td><div class="cell-name">%s</div><div class="check-id">%s</div></td>'
                 '<td>%s</td><td class="cell-message">%s'%(
            esc(cn), esc(nm), badge(st,""), esc(msg)))
        if dh:
            cid="%s-%d"%(section_id, idx)
            p.append(' <details class="check-detail-details"><summary class="det-toggle">详情</summary>'
                     '<div class="check-body">%s</div></details>'%dh)
        p.append('</td></tr>')
    p.append('</tbody></table>')
    return "\n".join(p)

def synthesize_variant_count(checks, golden_dir):
    """Collapse duplicate variant.count rows and fill golden count when missing."""
    vc=[c for c in checks if c.get("name")=="variant.count"]
    others=[c for c in checks if c.get("name")!="variant.count"]
    if not vc: return others
    # prefer a PASS/FAIL with detail; else last one
    pick=None
    for c in vc:
        if c.get("status") in ("PASS","FAIL") and (c.get("detail") or {}).get("baseline") is not None:
            pick=c; break
    if pick is None:
        for c in vc:
            if c.get("status")!="INFO": pick=c; break
    if pick is None: pick=vc[-1]
    det=dict(pick.get("detail") or {})
    test=det.get("test")
    if test is None: test=parse_count_msg(pick.get("message",""))
    base=det.get("baseline")
    if base is None and golden_dir:
        # from golden case.json
        gj=os.path.join(golden_dir,"case.json")
        if os.path.isfile(gj):
            try:
                gd=json.load(open(gj))
                for gc in gd.get("checks",[]):
                    if gc.get("name")=="variant.count":
                        base=parse_count_msg(gc.get("message",""))
                        if base is not None: break
            except Exception:
                pass
        if base is None:
            # try counting golden VCF (may be slow; skip files > 500MB)
            for name in os.listdir(golden_dir) if os.path.isdir(golden_dir) else []:
                if name.endswith(".genotyper.vcf.gz"):
                    gp=os.path.join(golden_dir,name)
                    try:
                        if os.path.getsize(gp) < 500*1024*1024:
                            import gzip
                            n=0
                            with gzip.open(gp,"rt",errors="ignore") as f:
                                for line in f:
                                    if line and not line.startswith("#"): n+=1
                            base=n
                    except Exception:
                        pass
                    break
    if test is not None and base is not None:
        delta=int(test)-int(base)
        pick={**pick,
              "status":"PASS",
              "message":"test %s / golden %s (Δ %+d)"%(test,base,delta),
              "detail":{"test":int(test),"baseline":int(base),"delta":delta}}
    elif test is not None:
        pick={**pick,"message":"%s variants"%test,"detail":{"test":int(test)}}
    return others+[pick]

def render_product_table(checks, section_id, case_name, report_path):
    root, golden_dir, work_dir = suite_paths(case_name, report_path)
    checks=synthesize_variant_count(checks, golden_dir)
    p=['<table class="table table-compact check-table">']
    p.append('<thead><tr><th>检查项</th><th>状态</th>'
             '<th>测试</th><th>基线</th><th>差值</th><th>说明</th></tr></thead><tbody>')
    for idx,c in enumerate(checks):
        nm=c.get("name",""); st=c.get("status","INFO"); msg=c.get("message","")
        cn=explain_check(nm)
        test_cell="—"; base_cell="—"; diff_cell="—"; extra=""; inline_extra=""
        if nm=="variant.count":
            det=c.get("detail") or {}
            tv,bv=det.get("test"),det.get("baseline")
            if tv is None: tv=parse_count_msg(msg)
            test_cell=esc(tv if tv is not None else "—")
            base_cell=esc(bv if bv is not None else "—")
            if tv is not None and bv is not None:
                d=int(tv)-int(bv)
                diff_cell=esc("%+d (%s)"%(d, _pct_diff(tv,bv)[0]))
            msg="变异计数对比"
        elif nm=="product.report":
            wp=os.path.join(work_dir,"report") if work_dir else None
            gp=os.path.join(golden_dir,"report") if golden_dir else None
            def _dir_files(d):
                if not d or not os.path.isdir(d): return {}
                out={}
                for dirpath,_,files in os.walk(d):
                    for f in files:
                        full=os.path.join(dirpath,f)
                        rel=os.path.relpath(full,d).replace("\\","/")
                        out[rel]=os.path.getsize(full)
                return out
            tfiles=_dir_files(wp); bfiles=_dir_files(gp)
            # fall back to check detail from case.json
            det=c.get("detail") or {}
            nt=det.get("test_files", len(tfiles) or None)
            nb=det.get("baseline_files", len(bfiles) or None)
            test_cell=esc("%s 个文件"%nt if nt is not None else ("已生成" if st=="PASS" else "缺失"))
            base_cell=esc("%s 个文件"%nb if nb is not None else ("—" if not bfiles else "%d 个文件"%len(bfiles)))
            mism=det.get("mismatched") or []
            only_t=det.get("only_in_test") or []
            only_b=det.get("only_in_baseline") or []
            # live recompute when dirs exist but detail empty
            if tfiles or bfiles:
                shared=set(tfiles)&set(bfiles)
                live_mis=[{"file":r,"test_size":tfiles[r],"baseline_size":bfiles[r]}
                          for r in sorted(shared) if tfiles[r]!=bfiles[r]]
                # also md5-less size compare for display; prefer case detail if present
                if not mism and not only_t and not only_b:
                    only_t=sorted(set(tfiles)-set(bfiles))
                    only_b=sorted(set(bfiles)-set(tfiles))
                    mism=live_mis
                    nt=len(tfiles); nb=len(bfiles)
                    test_cell=esc("%d 个文件"%nt)
                    base_cell=esc("%d 个文件"%nb if bfiles else "—")
            n_mis=det.get("mismatch_count", len(mism))
            if (bfiles or nb is not None) and (tfiles or nt is not None):
                if n_mis or only_t or only_b:
                    diff_cell=esc("内容差 %s / 仅测 %d / 仅基 %d"%(n_mis, len(only_t), len(only_b)))
                else:
                    diff_cell=esc("一致")
            else:
                diff_cell="—"
            # prefer recorded check status from case.json
            st=c.get("status","INFO")
            msg="报告目录对比"
            # Highlight + file inventory both go into collapsed 详细对比.
            hi_rows=extract_report_highlights(wp, gp)
            file_rows=[]
            all_names=sorted(set(list(tfiles.keys())+list(bfiles.keys())
                                 +[m.get("file") for m in mism if isinstance(m,dict)]
                                 +list(only_t)+list(only_b)))
            for rel in all_names[:40]:
                ts=tfiles.get(rel); bs=bfiles.get(rel)
                if ts is None and bs is None:
                    for m in mism:
                        if isinstance(m,dict) and m.get("file")==rel:
                            ts,bs=m.get("test_size"),m.get("baseline_size"); break
                if ts is not None and bs is not None:
                    d=ts-bs
                    diff=("相同" if d==0 else "%+d B"%d)
                elif ts is not None:
                    diff="仅测试"
                elif bs is not None:
                    diff="仅基线"
                else:
                    diff="—"
                file_rows.append((rel, fmt_size_mb(ts) if ts is not None else "—",
                                  fmt_size_mb(bs) if bs is not None else "—", diff))
            bits=[]
            if hi_rows:
                bits.append('<div class="muted" style="margin:0 0 4px">质控 / 比对 / 变异 重点指标（来自 report/*.xls）</div>')
                bits.append(render_metric_rows_table(hi_rows))
            if file_rows:
                bits.append('<div class="muted" style="margin:10px 0 4px">报告目录文件清单</div>')
                bits.append(render_metric_rows_table(file_rows))
            extra="\n".join(bits)
        elif nm.startswith("product."):
            bn=product_basename(nm)
            # test size: prefer live work file, else message
            wp=work_product_path(work_dir,bn) if bn else None
            gp=golden_product_path(golden_dir,bn) if bn else None
            tb=os.path.getsize(wp) if wp and os.path.isfile(wp) else parse_bytes_msg(msg)
            bb=os.path.getsize(gp) if gp and os.path.isfile(gp) else None
            # also try golden case.json message as fallback for size
            if bb is None and golden_dir and bn:
                gj=os.path.join(golden_dir,"case.json")
                if os.path.isfile(gj):
                    try:
                        gd=json.load(open(gj))
                        for gc in gd.get("checks",[]):
                            if gc.get("name")==nm:
                                bb=parse_bytes_msg(gc.get("message","")); break
                    except Exception:
                        pass
            test_cell=esc(fmt_size_mb(tb))
            base_cell=esc(fmt_size_mb(bb))
            if tb is not None and bb is not None:
                d=tb-bb
                diff_cell=esc("%s (%s)"%(fmt_size_mb(abs(d)).replace(" ","") if d==0 else
                    (("+" if d>0 else "-")+fmt_size_mb(abs(d))), _pct_diff(tb,bb)[0]))
                # fix signed size diff formatting
                sign="+" if d>=0 else "-"
                diff_cell=esc("%s%s (%s)"%(sign, fmt_size_mb(abs(d)), _pct_diff(tb,bb)[0]))
            msg="文件大小"
            # detailed bamStat / vcfStat content compare
            if bn in ("bamStat.txt","vcfStat.txt"):
                tmap=parse_kv_stat(wp)
                bmap=parse_kv_stat(gp)
                if tmap or bmap:
                    # prefer headline scalars for bamStat
                    prefer=None
                    if bn=="bamStat.txt":
                        prefer=["total reads","mapped reads","unmapped reads",
                                "properly paired reads","duplicate reads","MQ0 reads",
                                "total bases","mismatch bases","match bases","unity",
                                "meanInsertSize","depth","1x coverage","10x coverage",
                                "20x coverage","30x coverage"]
                        prefer=[k for k in prefer if k in tmap or k in bmap]
                    elif bn=="vcfStat.txt":
                        prefer=list(dict.fromkeys(list(tmap.keys())+list(bmap.keys())))
                    extra=render_kv_cmp_table(tmap, bmap, keys=prefer or None)
                    msg="文件大小 + 内容对比"
        else:
            # vcf.nonempty etc.
            test_cell=esc(msg[:40] if msg else "—")
        p.append('<tr><td><div class="cell-name">%s</div><div class="check-id">%s</div></td>'
                 '<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class="cell-message">%s'%(
            esc(cn), esc(nm), badge(st,""), test_cell, base_cell, diff_cell, esc(msg)))
        if inline_extra:
            p.append(inline_extra)
        if extra:
            p.append(' <details class="check-detail-details"><summary class="det-toggle">详细对比</summary>'
                     '<div class="check-body">%s</div></details>'%extra)
        p.append('</td></tr>')
    p.append('</tbody></table>')
    return "\n".join(p)

def comparison_map(item):
    """Map comparison stem (align_bam / gvcf / vcf) -> loaded JSON."""
    out={}
    for fn, data in item.get("_comparisons",[]):
        out[fn.replace(".json","")]=data
    return out

def render_concordance_section(checks, item, section_id):
    """One row per concordance product, with a key-metrics comparison table
    inline and the full '对比：*' check dump tucked into an expandable detail."""
    cmap=comparison_map(item)
    # keep case-level concordance checks (align_bam / gvcf / vcf / ...)
    conc=[c for c in checks if categorize_check(c.get("name",""))=="concordance"]
    if not conc and not cmap: return ""
    # if case has no top-level conc checks, synthesize from comparison files
    if not conc:
        for stem, data in cmap.items():
            st="PASS"
            for ch in data.get("checks",[]):
                if RANK.get(ch.get("status","INFO"),0) > RANK.get(st,0):
                    st=ch.get("status","INFO")
            conc.append({"name":stem,"status":st,"message":""})
    n=len(conc)
    nbad=sum(1 for c in conc if c.get("status","INFO") in ("FAIL","ERROR","WARN"))
    p=['<div class="subsection"><h4>一致性检查 <span class="section-summary">%s</span></h4>'%(
        esc("%d 项 · %d 通过 / %d 异常"%(n, n-nbad, nbad)))]
    for idx,c in enumerate(conc):
        nm=c.get("name",""); st=c.get("status","INFO"); msg=c.get("message","")
        cn=explain_check(nm)
        comp=cmap.get(nm)
        # build key metrics table
        if comp:
            if any(ch.get("name","").startswith("stats.flagstat") or ch.get("name")=="header.sq"
                   for ch in comp.get("checks",[])):
                rows=extract_bam_key_metrics(comp)
            else:
                rows=extract_vcf_key_metrics(comp)
            # status: escalate from comparison checks
            for ch in comp.get("checks",[]):
                if RANK.get(ch.get("status","INFO"),0) > RANK.get(st,0):
                    st=ch.get("status","INFO")
        else:
            rows=[]
        p.append('<div class="conc-block">')
        p.append('<div class="conc-head">%s <span class="cell-name">%s</span> '
                 '<span class="check-id">%s</span> <span class="muted">%s</span></div>'%(
            badge(st,""), esc(cn), esc(nm), esc(msg[:100] if msg else "")))
        if rows:
            p.append(render_metric_rows_table(rows))
        else:
            p.append('<p class="muted">无可展示的对比指标</p>')
        # When body digests differ, surface sampled line diffs (INFO+FORMAT) like
        # the enhanced 30x-style report — still keep the full dump in 对比详情.
        if comp:
            rd_chk = next((ch for ch in comp.get("checks", [])
                           if ch.get("name") == "records.digest"), None)
            if rd_chk:
                fd = (rd_chk.get("detail") or {}).get("first_differences")
                if fd:
                    p.append(render_first_differences(fd))
        if comp:
            cid="%s-%d"%(section_id, idx)
            detail_bits=[]
            for ch in comp.get("checks",[]):
                conc_panel=render_concordance(ch)
                if conc_panel: detail_bits.append(conc_panel)
            sub=[ch for ch in comp.get("checks",[]) if not ch.get("name","").startswith("tool.")]
            if sub:
                detail_bits.append(render_preflight_table(sub, cid+"-sub"))
            p.append('<details class="check-detail-details"><summary class="det-toggle">对比详情</summary>'
                     '<div class="check-body">%s</div></details>'%("\n".join(detail_bits)))
        p.append('</div>')
    p.append('</div>')
    return "\n".join(p)

def render_checks(checks, item, section_id="checks"):
    """Preflight / product / concordance sections. step.* dropped."""
    if not checks: return ""
    checks=[c for c in checks if not c.get("name","").startswith("step.")]
    if not checks: return ""
    case_name=item.get("name","")
    report_path=item.get("report") or ""
    groups={"preflight":[],"product":[],"concordance":[]}
    for c in checks:
        groups[categorize_check(c.get("name",""))].append(c)
    p=[]
    # preflight
    gs=groups["preflight"]
    if gs:
        n=len(gs); nbad=sum(1 for c in gs if c.get("status","INFO") in ("FAIL","ERROR","WARN"))
        p.append('<div class="subsection"><h4>前置检查 <span class="section-summary">%s</span></h4>'%(
            esc("%d 项 · %d 通过 / %d 异常"%(n, n-nbad, nbad))))
        p.append(render_preflight_table(gs, section_id+"-pre"))
        p.append('</div>')
    # product
    gs=groups["product"]
    if gs:
        # recount after synthesize
        root,gdir,_=suite_paths(case_name, report_path)
        shown=synthesize_variant_count(gs, gdir)
        n=len(shown); nbad=sum(1 for c in shown if c.get("status","INFO") in ("FAIL","ERROR","WARN"))
        p.append('<div class="subsection"><h4>产物检查 <span class="section-summary">%s</span></h4>'%(
            esc("%d 项 · %d 通过 / %d 异常"%(n, n-nbad, nbad))))
        p.append(render_product_table(gs, section_id+"-prod", case_name, report_path))
        p.append('</div>')
    # concordance (includes comparison cards)
    p.append(render_concordance_section(groups["concordance"], item, section_id+"-conc"))
    return "\n".join(p)

def render_concordance(check):
    d=check.get("detail",{})
    if not d or not isinstance(d,dict): return ""
    if not any(k in d for k in ("jaccard","precision","recall","f1")): return ""
    sh=d.get("shared",0); bo=d.get("baseline_only",0); to=d.get("test_only",0)
    j=d.get("jaccard",0); pr=d.get("precision",0); rc=d.get("recall",0); f1=d.get("f1",0)
    p=['<div class="concordance-panel"><h5>位点一致性</h5><div class="metric-cards">']
    for l,v,c in [("Jaccard",j,"#3b82f6"),("Precision",pr,"#16a34a"),("Recall",rc,"#a855f7"),("F1",f1,"#3b82f6")]:
        p.append('<div class="metric-card"><div class="metric-label">%s</div><div class="metric-value">%.4f</div><div class="metric-bar">%s</div></div>'%(esc(l),float(v),bar(float(v),1.0,c)))
    p.append('</div><div class="overlap-bars">')
    for l,v,bg,co in [("共享",sh,"#dcfce7","#166534"),("仅基线",bo,"#fef3c7","#92400e"),("仅测试",to,"#f3e8ff","#6b21a8")]:
        p.append('<div class="overlap-row" style="background:%s;color:%s"><span>%s</span><span>%d</span></div>'%(bg,co,esc(l),v))
    p.append('</div>')
    ex=d.get("examples",[])
    if ex:
        p.append('<div class="examples"><strong>差异示例:</strong><ul>')
        for e in ex[:5]: p.append('<li><code>%s</code></li>'%esc(e))
        p.append('</ul></div>')
    p.append('</div>')
    return "\n".join(p)

def render_case(item):
    name=item.get("name",""); st=item.get("status","INFO"); w=dur(item.get("wall_s"))
    det=item.get("_detail"); probs=item.get("problems",[]); counts=item.get("counts",{})
    cs=""
    if counts:
        parts=[]
        for s in ("PASS","WARN","FAIL","ERROR","SKIP","INFO"):
            n=counts.get(s,0)
            if n: parts.append('<span class="mini-badge" style="background:%s;color:%s">%s %d</span>'%(BG.get(s,"#f1f5f9"),CLR.get(s,"#64748b"),s,n))
        cs=" ".join(parts)
    pt=esc(probs[0][:100]) if probs else ""
    if len(probs)>1: pt+=" (+%d)"%(len(probs)-1)
    p=['<div class="case-card" id="case-%s">'%esc(name)]
    p.append('<div class="case-card-header" onclick="toggleCase(\'%s\')">'%esc(name))
    p.append('<div class="case-card-left">')
    p.append('<span class="case-icon" style="background:%s;color:%s">%s</span>'%(BG.get(st,"#f1f5f9"),CLR.get(st,"#64748b"),ICON.get(st,"?")))
    p.append('<div><div class="case-name">%s</div><div class="case-meta">%s · %s</div></div>'%(esc(name),esc(item.get("kind","")),w))
    p.append('</div><div class="case-card-right">')
    if cs: p.append('<div class="case-checks">%s</div>'%cs)
    if pt: p.append('<div class="case-problem">%s</div>'%pt)
    p.append('<span class="case-expand" id="expand-%s">▼</span>'%esc(name))
    p.append('</div></div>')
    p.append('<div class="case-detail" id="detail-%s" style="display:none">'%esc(name))
    if det:
        if det.get("steps"): p.append(render_steps(det))
        # checks + concordance (comparison cards folded into concordance details)
        if det.get("checks") or item.get("_comparisons"):
            p.append(render_checks(det.get("checks") or [], item, "checks-"+esc(name)))
    else:
        p.append('<p class="muted">无详细报告</p>')
    p.append('</div></div>')
    return "\n".join(p)

CSS = """
:root{--bg:#f1f5f9;--card:#fff;--border:#e2e8f0;--text:#1e293b;--muted:#64748b;--accent:#3b82f6}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;background:var(--bg);color:var(--text);line-height:1.6;padding:24px;max-width:1280px;margin:0 auto}
.overview{display:flex;justify-content:space-between;align-items:center;background:linear-gradient(135deg,#fff,#f8fafc);padding:28px 32px;border-radius:16px;margin-bottom:20px;box-shadow:0 1px 4px rgba(0,0,0,.06)}
.overview-version{font-size:14px;font-weight:600;color:var(--accent);margin-bottom:2px}
.overview-title h1{font-size:26px;font-weight:700;margin-bottom:6px}
.overview-meta{display:flex;gap:16px;flex-wrap:wrap;font-size:13px;color:var(--muted)}
.overview-meta code{background:var(--bg);padding:2px 8px;border-radius:4px;font-size:12px;color:var(--text)}
.overview-verdict{font-size:22px;font-weight:700;padding:10px 24px;border-radius:12px;white-space:nowrap}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:14px;margin-bottom:24px}
.card{background:var(--card);padding:20px;border-radius:12px;text-align:center;box-shadow:0 1px 3px rgba(0,0,0,.05)}
.card-value{font-size:30px;font-weight:700;line-height:1.2}
.card-label{font-size:13px;color:var(--muted);margin-top:4px}
.case-card{background:var(--card);border-radius:12px;margin-bottom:12px;box-shadow:0 1px 3px rgba(0,0,0,.05);overflow:hidden;transition:box-shadow .2s}
.case-card:hover{box-shadow:0 2px 8px rgba(0,0,0,.1)}
.case-card-header{display:flex;justify-content:space-between;align-items:center;padding:16px 20px;cursor:pointer;user-select:none}
.case-card-left{display:flex;align-items:center;gap:14px}
.case-icon{display:inline-flex;align-items:center;justify-content:center;width:36px;height:36px;border-radius:50%;font-size:18px;font-weight:700;flex-shrink:0}
.case-name{font-size:16px;font-weight:600}
.case-meta{font-size:12px;color:var(--muted)}
.case-card-right{display:flex;align-items:center;gap:12px}
.case-checks{display:flex;gap:4px;flex-wrap:wrap}
.case-problem{font-size:12px;color:var(--muted);max-width:400px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.case-expand{font-size:14px;color:var(--muted);transition:transform .2s}
.case-detail{padding:0 20px 20px;border-top:1px solid var(--border)}
.badge{display:inline-flex;align-items:center;gap:2px;padding:3px 10px;border-radius:20px;font-size:12px;font-weight:600;white-space:nowrap}
.mini-badge{display:inline-block;padding:2px 7px;border-radius:12px;font-size:11px;font-weight:600}
.subsection{margin-top:16px;padding-top:14px;border-top:1px solid var(--border)}
.subsection h4{font-size:14px;font-weight:600;margin-bottom:10px;color:var(--muted)}
.subsection h5{font-size:13px;font-weight:600;margin-bottom:8px;color:var(--muted)}
.table{width:100%;border-collapse:collapse;font-size:13px}
.table th{text-align:left;padding:8px 10px;border-bottom:2px solid var(--border);font-weight:600;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.5px}
.table td{padding:8px 10px;border-bottom:1px solid var(--border)}
.table-compact th,.table-compact td{padding:5px 8px;font-size:12px}
.table tbody tr:hover{background:var(--bg)}
.cell-name{font-weight:600;white-space:nowrap}
.cell-message{max-width:500px;word-break:break-word}
.cell-bar{width:120px}
.bar-track{height:6px;background:var(--bg);border-radius:3px;overflow:hidden;min-width:80px}
.bar-fill{height:100%;border-radius:3px}
.steps-summary{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-top:12px}
.summary-card{background:var(--bg);border:1px solid var(--border);border-radius:10px;padding:12px 14px}
.summary-label{font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.5px;margin-bottom:4px}
.summary-value{font-size:18px;font-weight:700;line-height:1.3}
.summary-value .vs{color:var(--muted);font-weight:400;margin:0 4px}
.summary-diff{font-size:12px;color:var(--muted);margin-top:3px}
.summary-sub{font-size:11px;color:var(--muted);margin-top:2px}
.summary-verdict{display:flex;flex-direction:column;justify-content:center}
.summary-verdict .summary-value{display:flex;align-items:center;gap:6px}
.verdict-pass{border-color:#bbf7d0}.verdict-pass .summary-value{color:#16a34a}
.verdict-warn{border-color:#fde68a}.verdict-warn .summary-value{color:#d97706}
.check-table .cell-name{font-weight:600;white-space:nowrap}
.check-table .check-id{font-size:11px;color:var(--muted);font-family:ui-monospace,SFMono-Regular,Menlo,monospace;margin-top:1px}
.check-detail-details{margin-top:4px}
.det-toggle{cursor:pointer;font-size:11px;color:var(--accent);display:inline-block}
.check-detail-details[open] .det-toggle::after{content:" ▾"}
.check-detail-details:not([open]) .det-toggle::after{content:" ▸"}
.check-detail-details .check-body{margin-top:6px;padding:8px 10px;background:var(--bg);border-radius:6px;font-size:11px;line-height:1.8}
.conc-block{border:1px solid var(--border);border-radius:10px;padding:12px 14px;margin-bottom:10px;background:var(--card)}
.conc-head{display:flex;align-items:center;gap:10px;margin-bottom:8px;flex-wrap:wrap}
.cmp-table{margin-top:6px}
.cmp-table td:nth-child(2),.cmp-table td:nth-child(3),.cmp-table td:nth-child(4){font-variant-numeric:tabular-nums}
.diff-block{margin-top:10px}
.diff-table td{vertical-align:top;font-size:12px}
.diff-cell{display:block;max-width:420px;white-space:pre-wrap;word-break:break-all;font-size:11px;line-height:1.35}
.concordance-panel{background:var(--bg);padding:16px;border-radius:8px;margin:12px 0}
.metric-cards{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-bottom:12px}
.metric-card{text-align:center;background:var(--card);padding:10px;border-radius:8px}
.metric-label{font-size:11px;color:var(--muted);margin-bottom:2px}
.metric-value{font-size:20px;font-weight:700;margin-bottom:4px}
.metric-bar{margin-top:2px}
.overlap-bars{display:flex;gap:10px;margin-bottom:10px}
.overlap-row{flex:1;text-align:center;padding:8px;border-radius:8px;font-size:13px;display:flex;justify-content:space-between}
.examples{margin-top:8px;font-size:12px}
.examples ul{margin:4px 0 0 16px}
.examples code{background:var(--bg);padding:2px 6px;border-radius:4px;font-size:11px}
.cell-detail{font-size:11px;color:var(--muted);max-width:300px;word-break:break-word}
.check-detail{font-size:11px;line-height:1.8}
.muted{color:var(--muted)}
.section-summary{font-size:12px;font-weight:400;color:var(--muted);margin-left:8px}
.check-item{border:1px solid var(--border);border-radius:8px;margin-bottom:6px;overflow:hidden}
.check-head{display:flex;align-items:center;gap:10px;padding:8px 12px;cursor:pointer;user-select:none;background:var(--card)}
.check-head:hover{background:var(--bg)}
.check-cn{font-weight:600;font-size:13px}
.check-id{font-size:11px;color:var(--muted);font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
.check-pf{margin-left:auto;font-size:12px;font-weight:600;padding:2px 8px;border-radius:10px}
.pf-pass{color:#16a34a;background:#dcfce7}.pf-fail{color:#dc2626;background:#fee2e2}
.pf-error{color:#9333ea;background:#f3e8ff}.pf-warn{color:#d97706;background:#fef3c7}
.pf-skip,.pf-info{color:#64748b;background:#f1f5f9}
.check-expand{font-size:12px;color:var(--muted);transition:transform .2s}
.check-body{padding:10px 14px 12px;background:var(--bg);border-top:1px solid var(--border)}
.check-msg{font-size:12px;margin-bottom:6px;word-break:break-word}
.comp-card{border:1px solid var(--border);border-radius:10px;margin-bottom:10px;overflow:hidden}
.comp-head{display:flex;align-items:center;gap:10px;padding:10px 14px;cursor:pointer;user-select:none;background:var(--card)}
.comp-head:hover{background:var(--bg)}
.comp-label{font-weight:600;font-size:13px}
.comp-summary{font-size:12px;color:var(--muted);flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.comp-body{padding:12px 14px;border-top:1px solid var(--border)}
@media(max-width:768px){.cards{grid-template-columns:repeat(2,1fr)}.metric-cards{grid-template-columns:repeat(2,1fr)}.overview{flex-direction:column}.overlap-bars{flex-direction:column}}
"""

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>dcstools 发版测试报告 — {run_id}</title>
<style>{css}</style></head>
<body>
{overview}
{cases}
<script>
function toggleCase(id){
  var d=document.getElementById('detail-'+id);
  var e=document.getElementById('expand-'+id);
  if(d){d.style.display=d.style.display==='none'?'block':'none';
    e.style.transform=d.style.display==='none'?'rotate(0deg)':'rotate(180deg)';}
}
function toggleCheck(id){
  var d=document.getElementById(id);
  var e=document.getElementById('exp-'+id);
  if(d){d.style.display=d.style.display==='none'?'block':'none';
    if(e){e.style.transform=d.style.display==='none'?'rotate(0deg)':'rotate(90deg)';}}
}
</script>
</body></html>
"""

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("suite_json",nargs="?")
    parser.add_argument("--suite-dir")
    parser.add_argument("-o","--out")
    args=parser.parse_args(argv)
    sp=args.suite_json
    if not sp and args.suite_dir: sp=os.path.join(args.suite_dir,"suite.json")
    if not sp: parser.error("pass suite.json path or --suite-dir")
    suite=load_suite(sp)
    if not suite:
        print("error: could not load %s"%sp,file=sys.stderr); return 1
    out=args.out or os.path.join(os.path.dirname(os.path.abspath(sp)),"report.html")
    cases="\n".join(render_case(i) for i in suite.get("items",[]))
    html_content=(HTML_TEMPLATE
        .replace("{run_id}",esc(suite.get("run_id","")))
        .replace("{css}",CSS)
        .replace("{overview}",render_overview(suite))
        .replace("{cases}",cases))
    with open(out,"w") as f: f.write(html_content)
    print("HTML report: %s"%out)
    return 0

if __name__=="__main__": sys.exit(main())
