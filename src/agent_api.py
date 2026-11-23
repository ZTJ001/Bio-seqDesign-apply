#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""稳定的程序化调用接口（供任意 agent / 脚本通过 subprocess 调用）。

设计契约（对调用方承诺，底层报告结构变化不影响调用方）：
  - 输入：一个 JSON 对象，可通过 --request <file.json> 或 stdin 传入
  - 输出：**stdout 仅输出一个 JSON 对象**；所有日志/错误走 stderr
  - 顶层字段 status 恒为 "ok" 或 "error"

最小请求：
  {"sequence": "MKVIFVLLL..."}
或：
  {"fasta": ">my_prot\\nMKVIFVLLL..."}
或：
  {"input_file": "/abs/path/input.fasta"}

可选字段（全部有合理默认）：
  name, purpose, optimizer(auto|codonbert|rule), host, system, vector, promoter,
  mutations("M1T,C59A"), glycosylation, secreted, disulfide, toxic, stable_line,
  skip_predict, output_dir

命令行等价：
  python src/agent_api.py --request req.json
  echo '{"sequence":"MKV..."}' | python src/agent_api.py
  python src/agent_api.py --sequence MKV... --purpose reagent
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKFLOW = HERE / "design_workflow.py"
API_VERSION = "1.0"


def _emit(obj):
    sys.stdout.write(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")
    sys.stdout.flush()


def _build_request(args):
    """合并 JSON 请求与 CLI 覆盖参数，返回规范化的 dict。"""
    req = {}
    if args.request:
        req.update(json.loads(Path(args.request).read_text(encoding="utf-8")))
    elif not sys.stdin.isatty() and not args.sequence and not args.fasta:
        raw = sys.stdin.read().strip()
        if raw:
            req.update(json.loads(raw))
    # CLI 覆盖（显式给出时）
    for key in ("sequence", "fasta", "input_file", "name", "purpose", "optimizer",
                "host", "system", "vector", "promoter", "mutations", "output_dir"):
        val = getattr(args, key, None)
        if val is not None:
            req[key] = val
    for key in ("glycosylation", "secreted", "disulfide", "toxic", "stable_line", "skip_predict"):
        if getattr(args, key, False):
            req[key] = True
    return req


def _write_input(req, workdir):
    """把请求里的序列落地为 FASTA 文件，返回路径。"""
    if req.get("input_file"):
        p = Path(req["input_file"])
        if not p.is_file():
            raise ValueError(f"input_file 不存在: {p}")
        return p
    if req.get("fasta"):
        p = workdir / "input.fasta"
        p.write_text(req["fasta"] if req["fasta"].strip().startswith(">")
                     else ">query\n" + req["fasta"].strip() + "\n", encoding="utf-8")
        return p
    if req.get("sequence"):
        name = str(req.get("name") or "query").split()[0]
        p = workdir / "input.fasta"
        p.write_text(f">{name}\n{req['sequence'].strip()}\n", encoding="utf-8")
        return p
    raise ValueError("请求需提供 sequence / fasta / input_file 之一")


def _run(req):
    workdir = Path(req.get("output_dir") or tempfile.mkdtemp(prefix="bioseq_apply_"))
    workdir.mkdir(parents=True, exist_ok=True)
    fasta = _write_input(req, workdir)

    cmd = [sys.executable, str(WORKFLOW), "-i", str(fasta), "-o", str(workdir)]
    simple = {"purpose": "--purpose", "optimizer": "--optimizer", "host": "--host",
              "system": "--system", "vector": "--vector", "promoter": "--promoter",
              "mutations": "--mutations"}
    for key, flag in simple.items():
        if req.get(key):
            cmd += [flag, str(req[key])]
    flags = {"glycosylation": "--glycosylation", "secreted": "--secreted",
             "disulfide": "--disulfide", "toxic": "--toxic",
             "stable_line": "--stable-line", "skip_predict": "--skip-predict"}
    for key, flag in flags.items():
        if req.get(key):
            cmd.append(flag)

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"工作流失败(exit={proc.returncode}): {proc.stderr[-800:]}")

    report_path = workdir / "design_report.json"
    if not report_path.is_file():
        raise RuntimeError("工作流未生成 design_report.json: " + proc.stdout[-800:])
    report = json.loads(report_path.read_text(encoding="utf-8"))

    records = []
    warnings = []
    for rec in report.get("records", []):
        sel = rec.get("selection", {})
        dec = rec.get("decision", {})
        sm = rec.get("sequence_metrics", {})
        quality = rec.get("quality", {})
        val = rec.get("validation", {})
        pred = rec.get("prediction", {})

        mfe = None
        for r in quality.get("vienna_records", []) or []:
            if isinstance(r, dict) and r.get("mfe_kcal_per_mol") is not None:
                mfe = r["mfe_kcal_per_mol"]
                break
        cai = None
        for r in quality.get("jcat_records", []) or []:
            if isinstance(r, dict):
                for k in ("cai", "CAI", "score"):
                    if r.get(k) is not None:
                        cai = r[k]
                        break
            if cai is not None:
                break

        identity = sm.get("protein_identity_pct")
        identity_ok = (identity is None) or (float(identity) >= 100.0)
        if sm.get("identity_warning"):
            warnings.append(f"{rec.get('name')}: {sm['identity_warning']}")

        records.append({
            "name": rec.get("name"),
            "protein_length": rec.get("protein_length"),
            "recommendation": {
                "system": sel.get("system"), "vector": sel.get("vector"),
                "host": sel.get("host"), "promoter": sel.get("promoter"),
                "vector_in_whitelist": sel.get("vector_in_whitelist"),
                "rule_id": dec.get("rule_id"), "rule_name": dec.get("rule_name"),
                "confidence": dec.get("confidence"), "rationale": dec.get("rationale"),
            },
            "optimization": {
                "optimizer": rec.get("optimizer"),
                "length_nt": sm.get("length_nt"), "gc_percent": sm.get("gc_percent"),
                "stop_codon": sm.get("stop_codon"),
                "protein_identity_pct": identity, "identity_ok": identity_ok,
                "aa_differences": sm.get("aa_differences", []),
            },
            "quality": {
                "mfe_kcal_per_mol": mfe, "cai": cai,
                "validation_available": val.get("available"),
                "validation_summary": val.get("summary"),
            },
            "prediction": {
                "available": bool(pred.get("available")),
                "reason": pred.get("reason"),
                "note": pred.get("note"),
                "model": pred.get("model", []),
                "rows": pred.get("rows", []),
            },
        })

    return {
        "status": "ok",
        "api_version": API_VERSION,
        "input": {"name": report.get("input"), "n_records": len(records)},
        "records": records,
        "warnings": warnings,
        "artifacts": {
            "final_fasta": str(workdir / "final_optimized.fasta"),
            "report_json": str(report_path),
            "report_md": str(workdir / "design_report.md"),
        },
    }


def main():
    ap = argparse.ArgumentParser(description="Bio-seqDesign agent API（stdout 仅输出 JSON）")
    ap.add_argument("--request", help="JSON 请求文件路径")
    ap.add_argument("--sequence", help="蛋白序列（裸序列）")
    ap.add_argument("--fasta", help="FASTA 文本")
    ap.add_argument("--input-file", dest="input_file", help="FASTA 文件路径")
    ap.add_argument("--name")
    ap.add_argument("--purpose")
    ap.add_argument("--optimizer")
    ap.add_argument("--host")
    ap.add_argument("--system")
    ap.add_argument("--vector")
    ap.add_argument("--promoter")
    ap.add_argument("--mutations")
    ap.add_argument("--output-dir", dest="output_dir")
    ap.add_argument("--glycosylation", action="store_true")
    ap.add_argument("--secreted", action="store_true")
    ap.add_argument("--disulfide", action="store_true")
    ap.add_argument("--toxic", action="store_true")
    ap.add_argument("--stable-line", dest="stable_line", action="store_true")
    ap.add_argument("--skip-predict", dest="skip_predict", action="store_true")
    args = ap.parse_args()

    try:
        req = _build_request(args)
        _emit(_run(req))
    except Exception as exc:  # noqa: BLE001 - 对调用方承诺永远返回 JSON
        sys.stderr.write(f"[agent_api] {type(exc).__name__}: {exc}\n")
        _emit({"status": "error", "api_version": API_VERSION,
               "error": {"type": type(exc).__name__, "message": str(exc)}})
        sys.exit(1)


if __name__ == "__main__":
    main()
