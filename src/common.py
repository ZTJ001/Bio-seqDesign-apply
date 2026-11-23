#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""模型适配器的公共工具：路径解析、FASTA 读写、--check 契约、JSON 输出。

所有适配器遵循统一契约：
  --check  仅检查部署状态，输出 [就绪]/[缺失] 与安装指引，退出码 0
  --run    执行推理（默认行为），FASTA 进 -> JSON + FASTA 出
"""
import json
import os
import sys
from pathlib import Path


def project_dir():
    """交付包根目录。"""
    return Path(__file__).resolve().parent.parent

ROOT = project_dir()
THIRD_PARTY = ROOT / "third_party"
PY_OPT = Path(os.environ.get("BIO_SEQDESIGN_PYTHON") or sys.executable)
PY_ENGINE = PY_OPT
KB_FILE = project_dir() / "molecular_build_rules.yaml"

CODON_TABLES = {
    "F": ["TTT", "TTC"], "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"], "M": ["ATG"], "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"], "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"], "*": ["TAA", "TAG", "TGA"], "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"], "N": ["AAT", "AAC"], "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"], "E": ["GAA", "GAG"], "C": ["TGT", "TGC"],
    "W": ["TGG"], "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
}
STOP = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    """返回 [(name, seq), ...]，去除空白并把 U 统一为 T，序列转大写。"""
    records, name, chunks = [], None, []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                # 兼容引擎的注释分段 FASTA（5UTR:/CDS:/polyA: 等）——只取 CDS 段拼接
                if ":" in line and any(line.upper().startswith(p) for p in ("CDS:", "5UTR:", "3UTR:", "POLYA:")):
                    if line.upper().startswith("CDS:"):
                        chunks.append(line.split(":", 1)[1].strip())
                else:
                    chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return [(n, s.replace("U", "T").upper()) for n, s in records]


def translate(seq):
    """DNA 翻译为蛋白，遇到首个终止密码子停止（不含 *）。"""
    table = {}
    for aa, codons in CODON_TABLES.items():
        for c in codons:
            table[c] = aa
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = table.get(seq[i:i + 3], "X")
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


def load_kb():
    import yaml
    with open(KB_FILE, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def host_preferred(host="homo_sapiens"):
    """取宿主偏好密码子映射 aa->codon；未知宿主回退人源。"""
    kb = load_kb()
    bias = kb.get("codon_usage_bias", {})
    if host not in bias:
        alias = kb.get("codon_optimization", {}).get("host_alias", {})
        host = alias.get(host, "homo_sapiens")
    block = bias.get(host, bias.get("homo_sapiens", {}))
    return block.get("preferred_codons", {}), host


def write_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2)


def write_fasta(path, records):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(">{}\n{}\n".format(name, seq))


def check_report(items, name):
    """items: [(label, ok, detail)]，打印统一的部署检查。"""
    print("[{}] 部署检查".format(name))
    for label, ok, detail in items:
        print("  [{}] {}: {}".format("就绪" if ok else "缺失", label, detail))
    return all(x[1] for x in items)
