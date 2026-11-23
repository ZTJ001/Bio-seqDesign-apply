#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Jcat 适配器（本地 Python 等价实现）：按宿主密码子偏好适配 + CAI 打分。

背景：Jcat 官方是 Java web/jar 工具，本机无 Java 运行时。本适配器用项目知识库的
5 宿主 Kazusa 偏好密码子表实现等价核心功能：
  - score：计算输入 CDS 的宿主偏好命中率与简化 CAI（相对 adaptiveness）
  - adapt：按每个氨基酸的宿主偏好密码子做同义替换（不改蛋白），输出优化 CDS

契约：--check / --run；-i FASTA -> JSON + FASTA；--host 指定宿主（默认人源）。
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (read_fasta, CODON_TABLES, STOP, translate, host_preferred,
                    check_report, write_json, write_fasta, KB_FILE)

AA_TO_CODONS = {}
for _aa, _codons in CODON_TABLES.items():
    if _aa != "*":
        AA_TO_CODONS[_aa] = _codons
CODON_TO_AA = {}
for _aa, _codons in CODON_TABLES.items():
    for _c in _codons:
        CODON_TO_AA[_c] = _aa


def do_check():
    ok = KB_FILE.exists()
    print("[Jcat(Python)] 部署检查")
    print("  [{}] 知识库(密码子表): {}".format("就绪" if ok else "缺失", KB_FILE))
    print("  [说明] 纯 Python 等价实现，无需 Java / 官方 jar")
    return 0


def _adaptiveness(host):
    """从知识库构建每个氨基酸内密码子的相对 adaptiveness（偏好=1，其他按 0.5 简化）。"""
    pref, _ = host_preferred(host)
    adapt = {}
    for aa, codons in AA_TO_CODONS.items():
        best = pref.get(aa)
        for c in codons:
            adapt[c] = 1.0 if c == best else 0.5
    return adapt


def score_cds(cds, host):
    adapt = _adaptiveness(host)
    pref, _ = host_preferred(host)
    codons = [cds[i:i + 3] for i in range(0, len(cds) - 2, 3) if cds[i:i + 3] not in STOP]
    matched = sum(1 for c in codons if pref.get(CODON_TO_AA.get(c, "")) == c)
    import math
    cai = math.exp(sum(math.log(adapt.get(c, 0.5)) for c in codons) / len(codons)) if codons else 0.0
    return {
        "host": host,
        "codons": len(codons),
        "preferred_matched": matched,
        "preferred_hit_rate": round(matched / len(codons), 4) if codons else 0.0,
        "simplified_cai": round(cai, 4),
    }


def adapt_cds(cds, host, keep_stop=True):
    pref, _ = host_preferred(host)
    out = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i + 3]
        aa = CODON_TO_AA.get(codon, "")
        if codon in STOP:
            out.append("TAA" if keep_stop else "")
        elif aa in pref:
            out.append(pref[aa])
        else:
            out.append(codon)
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description="Jcat(Python) 宿主密码子适配与 CAI 打分")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("-i", "--input", help="输入 FASTA")
    ap.add_argument("-o", "--output", default="jcat_report.json")
    ap.add_argument("--fasta-out", help="输出适配后 CDS FASTA")
    ap.add_argument("--host", default="homo_sapiens")
    ap.add_argument("--score-only", action="store_true", help="只打分，不做替换")
    args = ap.parse_args()

    if args.check:
        return do_check()
    if not args.input:
        ap.error("运行模式需要 -i/--input FASTA")

    records = read_fasta(args.input)
    results = []
    adapted = []
    for name, seq in records:
        before = score_cds(seq, args.host)
        item = {"name": name, "length_nt": len(seq), "protein": translate(seq), "before": before}
        if not args.score_only:
            new_seq = adapt_cds(seq, args.host)
            item["adapted"] = score_cds(new_seq, args.host)
            item["changed_nt"] = sum(a != b for a, b in zip(seq, new_seq))
            item["protein_unchanged"] = translate(new_seq) == translate(seq)
            adapted.append((name + "|jcat_adapted", new_seq))
        results.append(item)

    write_json(args.output, {"tool": "Jcat(Python)", "host": args.host, "count": len(results), "results": results})
    if args.fasta_out and adapted:
        write_fasta(args.fasta_out, adapted)
    print("Jcat 分析完成 {} 条 -> {}".format(len(results), args.output))
    for r in results:
        print("  {}: 命中率 {:.1%} -> {:.1%}, CAI {:.3f} -> {:.3f}".format(
            r["name"], r["before"]["preferred_hit_rate"],
            r.get("adapted", r["before"])["preferred_hit_rate"],
            r["before"]["simplified_cai"], r.get("adapted", r["before"])["simplified_cai"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
