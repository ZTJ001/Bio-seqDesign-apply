#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ViennaRNA 适配器：RNA 二级结构 / 最小自由能(MFE)预测。

契约：
  --check              检查 ViennaRNA Python 模块是否可用
  --run (默认)          FASTA(CDS/RNA/DNA) -> 每条序列的 dot-bracket 结构 + MFE + GC + AUG下游窗口MFE

知识库 ai_hook.rna_fold_predict：AUG 下游 10–15 bp ΔG < -10 kcal/mol 告警。
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_fasta, PY_OPT, check_report, write_json, write_fasta, ROOT


def _has_rna():
    code = "import RNA; s,m = RNA.fold('GGGAACCC'); assert isinstance(m,float)"
    try:
        return subprocess.run([str(PY_OPT), "-c", code], capture_output=True).returncode == 0
    except OSError:
        return False


def do_check():
    ok = _has_rna()
    check_report([
        ("Python 环境", PY_OPT.exists(), str(PY_OPT)),
        ("ViennaRNA 模块", ok, "pip install ViennaRNA（在 codon_opt/.venv 中）"),
    ], "ViennaRNA")
    return 0


def window_mfe(seq, start=3, size=15):
    """AUG 下游窗口(默认从第4 nt 起 15 nt)的 MFE。"""
    try:
        import RNA
        w = seq[start:start + size]
        if len(w) < 6:
            return None, None
        struct, mfe = RNA.fold(w)
        return struct, round(float(mfe), 3)
    except Exception:
        return None, None


def analyze(seq):
    import RNA
    rna = seq.replace("T", "U")
    struct, mfe = RNA.fold(rna)
    gc = (rna.count("G") + rna.count("C")) / len(rna) * 100 if rna else 0
    w_struct, w_mfe = window_mfe(rna)
    alerts = []
    if w_mfe is not None and w_mfe < -10:
        alerts.append("AUG 下游 15nt ΔG={} kcal/mol (< -10)，可能形成强茎环影响起始翻译".format(w_mfe))
    return {
        "length_nt": len(seq),
        "gc_percent": round(gc, 2),
        "mfe_kcal_per_mol": round(float(mfe), 3),
        "structure": struct,
        "aug_downstream_window_mfe": w_mfe,
        "aug_downstream_structure": w_struct,
        "alerts": alerts,
    }


def main():
    ap = argparse.ArgumentParser(description="ViennaRNA MFE 二级结构适配器")
    ap.add_argument("--check", action="store_true", help="仅检查部署状态")
    ap.add_argument("-i", "--input", help="输入 FASTA（蛋白或 DNA/RNA）")
    ap.add_argument("-o", "--output", default="mfe_report.json", help="输出 JSON 报告")
    ap.add_argument("--fasta-out", help="可选：输出带结构的 FASTA")
    args = ap.parse_args()

    if args.check:
        return do_check()

    if not _has_rna():
        print("[错误] ViennaRNA 不可用，请先运行: python viennarna_adapter.py --check", file=sys.stderr)
        return 2
    if not args.input:
        ap.error("运行模式需要 -i/--input FASTA")

    # 用 codon_opt 的解释器执行实际计算（ViennaRNA 装在该环境）
    code = r'''
import json, sys
sys.path.insert(0, r"{script_dir}")
import viennarna_adapter as va
from common import read_fasta, write_json, write_fasta
recs = read_fasta(r"{inp}")
res = []
for name, seq in recs:
    item = {{"name": name, "seq": seq}}
    item.update(va.analyze(seq))
    res.append(item)
write_json(r"{out}", {{"tool": "ViennaRNA", "count": len(res), "results": res}})
if r"{fa}":
    write_fasta(r"{fa}", [(r["name"] + "|MFE=" + str(r["mfe_kcal_per_mol"]), r["seq"]) for r in res])
print("已分析", len(res), "条序列 ->", r"{out}")
'''
    here = Path(__file__).resolve().parent
    filled = code.format(script_dir=str(here), inp=os.path.abspath(args.input),
                         out=os.path.abspath(args.output), fa=os.path.abspath(args.fasta_out) if args.fasta_out else "")
    rc = subprocess.run([str(PY_OPT), "-c", filled])
    return rc.returncode


if __name__ == "__main__":
    sys.exit(main())
