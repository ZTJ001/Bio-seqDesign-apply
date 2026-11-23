#!/usr/bin/env python3
"""CodonBERT 密码子优化封装（供 Agent / CLI 调用）。

默认使用本地微调权重（若存在），否则回退官方预训练权重。
输入蛋白 FASTA 或 DNA FASTA，输出优化后的密码子 DNA FASTA，并计算指标
（人源偏好命中率、GC、CAI 近似、蛋白一致性）。

用法:
    # 从蛋白 FASTA 优化
    python optimize_model.py -i protein.fasta -o optimized.fasta
    # 指定权重
    python optimize_model.py -i x.fasta -o y.fasta -m models/codonbert_finetuned.pt
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUNDLE = HERE.parent
CB_DIR = HERE / "codonbert"
PY = Path(os.environ.get("BIO_SEQDESIGN_PYTHON") or sys.executable)
PRETRAINED = CB_DIR / "models" / "kidney_1_1_CodonBert_model_20230726_320_model_param.pt"
FINETUNED = BUNDLE / "models" / "codonbert_finetuned.pt"
KB = BUNDLE / "molecular_build_rules.yaml"


def set_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-i", "--input", help="蛋白 FASTA（氨基酸字母）")
    ap.add_argument("-o", "--output", required=True, help="优化结果 FASTA")
    ap.add_argument("-m", "--model", help="权重路径，默认自动选择微调/预训练")
    return ap.parse_args()


def resolve_model(cli_model):
    if cli_model:
        return Path(cli_model)
    if FINETUNED.exists():
        return FINETUNED, "local_finetuned"
    return PRETRAINED, "official_pretrained"


def prepare_input(args, tmp_fasta):
    if not args.input:
        sys.exit("需提供 -i 输入 FASTA 或使用 --from-db")
    # 若输入是 DNA（含 ATGC 且长为3倍数），先翻译为蛋白
    text = Path(args.input).read_text(encoding="utf-8")
    seq_lines = [l.strip() for l in text.splitlines() if l and not l.startswith(">")]
    joined = "".join(seq_lines).upper()
    is_dna = set(joined) <= set("ATGCNU") and len(joined) % 3 == 0
    if not is_dna:
        Path(tmp_fasta).write_text(text, encoding="utf-8")
        return
    import yaml
    d = yaml.safe_load(KB.read_text(encoding="utf-8"))
    gc = d["codon_optimization"]["genetic_code"]
    c2a = {}
    for aa, codons in gc.items():
        if aa.startswith("_"):
            continue
        for c in codons:
            c2a[c] = aa
    for c in gc["_stop"]:
        c2a[c] = "*"
    out, name = [], ""
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                out.append(f">{name}\n" + "".join(cur).rstrip("*") + "\n")
            name, cur = line[1:].split()[0], []
        elif line.strip():
            s = line.strip().upper().replace("U", "T")
            cur += [c2a.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3)]
    if name:
        out.append(f">{name}\n" + "".join(cur).rstrip("*") + "\n")
    Path(tmp_fasta).write_text("".join(out), encoding="utf-8")


def metrics(fasta_path):
    import yaml
    d = yaml.safe_load(KB.read_text(encoding="utf-8"))
    gc = d["codon_optimization"]["genetic_code"]
    c2a = {}
    for aa, codons in gc.items():
        if aa.startswith("_"):
            continue
        for c in codons:
            c2a[c] = aa
    for c in gc["_stop"]:
        c2a[c] = "*"
    pref = d["codon_usage_bias"]["homo_sapiens"]["preferred_codons"]
    res = {}
    name, seq = None, []
    for line in Path(fasta_path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name:
                res[name] = "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip().upper().replace("U", "T"))
    if name:
        res[name] = "".join(seq)
    out = {}
    for n, s in res.items():
        codons = [s[i:i + 3] for i in range(0, len(s) - 2, 3)]
        coding = codons[:-1]
        hit = sum(1 for c in coding if pref.get(c2a.get(c, ""), "") == c)
        comp = sum(1 for c in coding if c2a.get(c) in pref)
        out[n] = {
            "length_nt": len(s),
            "gc": round((s.count("G") + s.count("C")) / len(s), 4) if s else 0,
            "human_pref_hit": round(hit / comp, 4) if comp else 0,
            "stop_codon": codons[-1] if codons else "",
        }
    return out


def main():
    args = set_args()
    model, model_tag = resolve_model(args.model)
    if not model.exists():
        sys.exit(f"模型权重不存在: {model}")
    tmp_in = HERE / "_model_input.fasta"
    prepare_input(args, tmp_in)

    out_abs = str(Path(args.output).resolve())
    cmd = [str(PY), str(CB_DIR / "predict.py"), "-m", str(model),
           "-f", str(tmp_in.resolve()), "-o", out_abs]
    print(f"[模型] {model_tag}: {model.name}", file=sys.stderr)
    subprocess.run(cmd, cwd=str(CB_DIR), check=True)

    m = metrics(args.output)
    print(json.dumps({"model": model_tag, "metrics": m}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
