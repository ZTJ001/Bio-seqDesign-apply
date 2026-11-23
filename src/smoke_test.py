#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""交付包自检：clone 后运行本脚本，确认环境齐全、两条优化路径可用。

用法： python src/smoke_test.py
退出码 0 = 全部通过；非 0 = 存在失败项。
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
API = HERE / "agent_api.py"
EXAMPLE = ROOT / "examples" / "example_protein.fasta"


def check(label, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {label}" + (f" — {detail}" if detail else ""))
    return bool(cond)


def main():
    ok = True
    print("=== Bio-seqDesign 交付包自检 ===")

    # 1. 依赖
    for mod in ("numpy", "sklearn", "joblib", "yaml", "Bio"):
        try:
            __import__(mod)
            ok &= check(f"依赖 {mod}", True)
        except Exception as e:
            ok &= check(f"依赖 {mod}", False, str(e))
    try:
        import torch  # noqa: F401
        ok &= check("依赖 torch", True)
    except Exception as e:
        ok &= check("依赖 torch", False, str(e))
    try:
        import RNA  # noqa: F401
        ok &= check("依赖 ViennaRNA(RNA)", True)
    except Exception as e:
        ok &= check("依赖 ViennaRNA(RNA)", False, str(e) + "（MFE 指标将不可用，但不阻断优化）")

    # 2. 关键文件
    for rel in ("molecular_build_rules.yaml", "config/target_reference.json",
                "models/effect_predictor.joblib", "models/binding_predictor.joblib",
                "models/codonbert_finetuned.pt", "examples/example_protein.fasta"):
        ok &= check(f"文件 {rel}", (ROOT / rel).is_file())

    with tempfile.TemporaryDirectory(prefix="bioseq_smoke_") as td:
        # 3. agent API 规则路径 + 目标蛋白预测
        req = {"input_file": str(EXAMPLE), "purpose": "reagent", "optimizer": "rule",
               "mutations": "C59A", "output_dir": str(Path(td) / "rule")}
        p = subprocess.run([sys.executable, str(API)],
                           input=json.dumps(req), text=True, capture_output=True)
        try:
            res = json.loads(p.stdout)
        except Exception:
            res = {}
        ok &= check("agent API 规则路径 + 目标蛋白预测",
                    res.get("status") == "ok" and len(res.get("records", [])) == 1,
                    (res.get("error") or {}).get("message", "") if res.get("status") != "ok" else "")
        if res.get("status") == "ok":
            rec = res["records"][0]
            ok &= check("  蛋白一致性 100%", rec["optimization"]["protein_identity_pct"] >= 100.0,
                        str(rec["optimization"]["protein_identity_pct"]))
            pred = rec["prediction"]
            ok &= check("  目标蛋白预测可用", bool(pred.get("available")) and len(pred.get("rows", [])) > 0)

        # 4. agent API CodonBERT 路径
        req2 = {"input_file": str(EXAMPLE), "purpose": "reagent", "optimizer": "codonbert",
                "output_dir": str(Path(td) / "cb")}
        p2 = subprocess.run([sys.executable, str(API)],
                            input=json.dumps(req2), text=True, capture_output=True)
        try:
            res2 = json.loads(p2.stdout)
        except Exception:
            res2 = {}
        ok &= check("agent API CodonBERT 路径",
                    res2.get("status") == "ok" and len(res2.get("records", [])) == 1,
                    (res2.get("error") or {}).get("message", "") if res2.get("status") != "ok" else "")

    print("===" + ("全部通过，可直接使用" if ok else "存在失败项，请检查环境") + "===")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
