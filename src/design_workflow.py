#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""端到端分子设计工作流：蛋白 -> 体系/载体推荐 -> 密码子优化 -> 质量指标 -> 校验 -> 表达/活性预测。

链路（对应知识库 molecular_build_rules.yaml#expression_system_recommendation）：
  1. 读蛋白/DNA FASTA，翻译校验
  2. 抽取蛋白特征（长度/MW/信号肽启发式/跨膜/二硫键潜力/糖基化基序/电荷）
  3. 按决策表推荐表达体系(mrna/mammalian/prokaryotic)、载体、宿主与启动子
  4. 密码子优化（--optimizer codonbert | rule | auto）
  5. 质量指标：GC / 宿主偏好命中率 / MFE(ViennaRNA) / CAI(Jcat) / 终止 / 蛋白一致性
  6. SOP 校验引擎复核（ERROR/WARNING/INFO）
  7. 表达/活性预测：
     - 目标蛋白及点突变：加载 effect_predictor + binding_predictor .joblib
     - 其它蛋白：明确标注「通用预测模型不可用」并给出扩展位（不编造结果）
  8. 输出 final_optimized.fasta + design_report.json + design_report.md

用法（用装有全部依赖的 Python 解释器运行，见 requirements.txt）：
  python design_workflow.py -i protein.fasta -o out_dir \
      [--purpose functional_assay] [--glycosylation] [--secreted] [--toxic] \
      [--mutations "M1T,C59A"] [--optimizer auto] [--system mammalian --vector pcDNA3.4]
"""
import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import (  # noqa: E402
    BUNDLE, SRC, MODELS, CONFIG, PY_OPT, KB_FILE, CODONBERT_DIR,
)

VALID_SYSTEMS = ("mrna", "mammalian", "prokaryotic")
AA = "ACDEFGHIKLMNPQRSTVWY"
HYDRO = set("AVILMFWYC")


# --------------------------------------------------------------------------- #
# FASTA / 翻译工具
# --------------------------------------------------------------------------- #
def read_fasta(path):
    out, name, buf = [], None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif line.strip() and not line.startswith(("5UTR:", "CDS:", "3UTR:", "polyA:")):
            buf.append(line.strip().upper().replace("U", "T"))
        elif line.strip().startswith("CDS:"):
            buf.append(line.strip()[4:].upper().replace("U", "T"))
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def codon_tables():
    import yaml
    kb = yaml.safe_load(KB_FILE.read_text(encoding="utf-8"))
    gc = kb["codon_optimization"]["genetic_code"]
    c2a = {c: aa for aa, codons in gc.items()
           if not aa.startswith("_") for c in codons}
    c2a.update({c: "*" for c in gc["_stop"]})
    a2c = {aa: codons for aa, codons in gc.items() if not aa.startswith("_")}
    return c2a, a2c, kb


def translate(seq, c2a):
    return "".join(c2a.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def to_protein(raw, c2a):
    """判定输入是蛋白还是 DNA，返回蛋白序列（去末端 *）。"""
    s = raw.replace("*", "").strip()
    if len(s) % 3 == 0 and set(s) <= set("ATGCN"):
        return translate(s, c2a).rstrip("*")
    return re.sub("[^A-Z]", "", s)


# --------------------------------------------------------------------------- #
# 蛋白特征提取（启发式）
# --------------------------------------------------------------------------- #
def extract_features(protein):
    n = len(protein)
    nterm = protein[:30]
    hydro_run = re.search(r"[AVILMFWYC]{7,}", nterm) is not None
    tm = len(re.findall(r"[AVILMFWYC]{18,}", protein))
    glyco = len(re.findall(r"N[^P][ST]", protein))
    return {
        "length_aa": n,
        "molecular_weight_kda": round(n * 110 / 1000, 1),
        "signal_peptide": round(sum(a in HYDRO for a in nterm) / max(len(nterm), 1), 2) > 0.5 and hydro_run,
        "has_transmembrane": tm > 0,
        "transmembrane_helices": tm,
        "disulfide_potential": protein.count("C") // 2,
        "glycosylation_motifs": glyco,
        "charged_fraction": round(sum(a in "DEKR" for a in protein) / max(n, 1), 3),
    }


# --------------------------------------------------------------------------- #
# 决策表匹配
# --------------------------------------------------------------------------- #
def recommend_system(features, user, kb):
    rules = kb["expression_system_recommendation"]["decision_rules"]
    ctx = dict(features)
    ctx.update({k: v for k, v in user.items() if v is not None})
    ctx.setdefault("secreted", features["signal_peptide"])
    ctx.setdefault("has_transmembrane", features["has_transmembrane"])
    ctx.setdefault("requires_glycosylation", bool(features["glycosylation_motifs"] >= 2))
    ctx.setdefault("requires_disulfide", features["disulfide_potential"] >= 2)
    ctx.setdefault("purpose", "functional_assay")
    ctx.setdefault("toxic", False)
    ctx.setdefault("stable_line", False)

    def hit(cond):
        for key, val in cond.items():
            if key == "when_all" or key == "when_any":
                continue
            if key == "length_aa_min" and ctx["length_aa"] < val:
                return False
            if key == "length_aa_max" and ctx["length_aa"] > val:
                return False
            if key in ctx:
                cv = ctx[key]
                if isinstance(val, list):
                    if cv not in val:
                        return False
                elif cv != val:
                    return False
        return True

    for rule in sorted(rules, key=lambda r: r["priority"]):
        if not hit(rule.get("when_all", {})):
            continue
        any_cond = rule.get("when_any")
        if any_cond and not any(hit({k: v}) for k, v in any_cond.items()):
            continue
        rec = dict(rule["recommend"])
        return {
            "rule_id": rule["id"], "rule_name": rule["name"],
            "confidence": rule["confidence"], "rationale": rule["rationale"],
            "recommendation": rec, "context": ctx,
            "notes": rule.get("notes", ""),
        }
    raise RuntimeError("决策表无匹配规则，请检查 knowledge base。")


def check_vector_whitelist(system, vector, kb):
    allowed = kb["system_list"][system]["sequence_design"]["backbone"]["allowed_vector"]
    return vector in allowed, allowed


# --------------------------------------------------------------------------- #
# 密码子优化
# --------------------------------------------------------------------------- #
def reverse_translate(protein, host, a2c, kb):
    """规则法反向翻译：每个氨基酸取宿主偏好密码子，末端补 TAA。"""
    pref = kb["codon_usage_bias"].get(host, {}).get("preferred_codons", {})
    out = []
    for aa in protein:
        if aa in pref:
            out.append(pref[aa])
        elif aa in a2c:
            out.append(a2c[aa][0])
        else:
            out.append("GCC")
    return "".join(out) + "TAA"


def optimize_codonbert(protein_fasta, raw_out):
    script = SRC / "optimize_model.py"
    cmd = [str(PY_OPT), str(script), "-i", str(protein_fasta), "-o", str(raw_out)]
    proc = subprocess.run(cmd, cwd=str(SRC), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("CodonBERT 优化失败: " + proc.stdout + proc.stderr)
    return proc.stdout


def choose_optimizer(requested, system, host):
    if requested != "auto":
        return requested
    return "rule" if host == "escherichia_coli" else "codonbert"


# --------------------------------------------------------------------------- #
# 质量指标（子进程调用适配器，避免主环境缺依赖）
# --------------------------------------------------------------------------- #
def run_script(script_rel, args, cwd=None):
    cmd = [str(PY_OPT), str(SRC / script_rel)] + args
    proc = subprocess.run(cmd, cwd=str(SRC), capture_output=True, text=True)
    if proc.returncode != 0:
        return {"available": False, "error": (proc.stdout + proc.stderr)[-500:]}
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return {"available": True, "raw": proc.stdout[-500:]}


def quality_metrics(cds, host, tmpdir, name):
    cds_fasta = tmpdir / f"{name}_cds.fasta"
    cds_fasta.write_text(f">{name}\n{cds}\n", encoding="utf-8")

    mfe = run_script("viennarna_adapter.py",
                     ["-i", str(cds_fasta), "-o", str(tmpdir / f"{name}_mfe.json"), "--fasta-out", str(tmpdir / f"{name}_struct.fasta")])
    cai = run_script("jcat_adapter.py",
                     ["-i", str(cds_fasta), "-o", str(tmpdir / f"{name}_jcat.json"),
                      "--host", host, "--score-only"])
    out = {"vienna": mfe, "jcat": cai}
    # 从报告文件读取结构化结果
    mfe_path = tmpdir / f"{name}_mfe.json"
    cai_path = tmpdir / f"{name}_jcat.json"
    if mfe_path.exists():
        try:
            d = json.load(open(mfe_path, encoding="utf-8"))
            out["vienna_records"] = d.get("results", [])
        except Exception:
            pass
    if cai_path.exists():
        try:
            d = json.load(open(cai_path, encoding="utf-8"))
            out["jcat_records"] = d.get("results", [])
        except Exception:
            pass
    return out


# --------------------------------------------------------------------------- #
# SOP 校验
# --------------------------------------------------------------------------- #
def wrap_construct(cds, system):
    """把裸 CDS 包装成校验引擎可识别的整段构建体：Kozak + CDS + 3UTR/polyA。"""
    utr5 = "GCCACC" if system in ("mrna", "mammalian") else "AAAGGAGG"
    utr3 = "AATAAA" + "A" * 60 if system in ("mrna", "mammalian") else "AATAAA"
    return ">" + "construct" + "\n" + utr5 + cds + utr3 + "\n"


def validate_construct(system, cds, vector, promoter, host, tmpdir, name):
    construct = tmpdir / f"{name}_construct.fasta"
    construct.write_text(wrap_construct(cds, system), encoding="utf-8")
    cmd = [str(PY_OPT), "-c",
           "import sys;sys.path.insert(0,%r);from check_engine import run_validation;"
           "import json;text=open(%r).read();"
           "print(json.dumps(run_validation(%r,text,vector=%r,promoter=%r,host=%r),ensure_ascii=False))"
           % (str(SRC), str(construct), system, vector, promoter, host)]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(SRC))
    if proc.returncode != 0:
        return {"available": False, "error": (proc.stdout + proc.stderr)[-500:]}
    try:
        d = json.loads(proc.stdout)
        return {"available": True, "summary": d.get("summary", {}),
                "error_items": [x for x in d.get("items", []) if x.get("level") == "ERROR"][:10],
                "codon_preference": d.get("codon_preference")}
    except json.JSONDecodeError:
        return {"available": False, "error": proc.stdout[-500:]}


# --------------------------------------------------------------------------- #
# 目标蛋白 表达/结合预测
# --------------------------------------------------------------------------- #
def load_featurizer():
    from featurizer import featurize
    return featurize


def parse_mutations(text):
    out = []
    for token in (text or "").split(","):
        token = token.strip()
        m = re.match(r"([A-Z])(\d+)([A-Z])", token)
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def predict_target(mutations):
    import numpy as np
    import joblib
    featurize = load_featurizer()
    expr_path = MODELS / "effect_predictor.joblib"
    bind_path = MODELS / "binding_predictor.joblib"
    expr = joblib.load(expr_path)
    bind = joblib.load(bind_path)
    X = np.array([featurize(w, s, m) for w, s, m in mutations], dtype=float)
    expr_score = expr["regressor"].predict(X)
    expr_harmful = expr["classifier"].predict_proba(X)[:, 1]
    bind_delta = bind["regressor"].predict(X)
    bind_lost = bind["classifier"].predict_proba(X)[:, 1]
    rows = []
    for i, (w, s, m) in enumerate(mutations):
        rows.append({
            "mutation": f"{w}{s}{m}",
            "expression_score": round(float(expr_score[i]), 3),
            "harmful_probability": round(float(expr_harmful[i]), 3),
            "binding_delta": round(float(bind_delta[i]), 3),
            "binding_lost_probability": round(float(bind_lost[i]), 3),
            "interpretation": (
                "结合丧失风险高" if bind_lost[i] >= 0.5 else
                "表达强有害风险高" if expr_harmful[i] >= 0.5 else "预测可耐受"
            ),
        })
    return {
        "available": True, "model": ["effect_predictor", "binding_predictor"],
        "note": "目标蛋白特异模型（酵母表达/受体结合 DMS 训练）；仅对 目标蛋白点突变有效，不可外推其它蛋白。",
        "rows": rows,
    }


def auto_detect_target(protein):
    """目标蛋白成熟肽 133aa 参考：高度相似则识别为 目标蛋白，用于自动提取点突变。"""
    ref = _TARGET_REF
    if not ref or abs(len(protein) - len(ref)) > 2:
        return None, []
    # 长度对齐（文件1 含 M1，成熟肽 133aa）
    a, b = protein, ref
    if len(a) == len(b) + 1 and a[0] == "M":
        a = a[1:]
    if len(a) != len(b):
        return None, []
    ident = sum(x == y for x, y in zip(a, b)) / len(b)
    if ident < 0.80:
        return None, []
    # DMS 文件1 编号含起始 M1，成熟肽第 i 位(0-based)对应 file1 的 site = i + 2
    muts = [(y, i + 2, x) for i, (x, y) in enumerate(zip(a, b)) if x != y]
    return round(ident, 3), muts


def load_target_ref():
    p = CONFIG / "target_reference.json"
    try:
        return json.load(open(p, encoding="utf-8")).get("mature_sequence", "")
    except Exception:
        return ""


_TARGET_REF = ""


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def main():
    global _TARGET_REF
    ap = argparse.ArgumentParser(description="端到端分子设计工作流（推荐→优化→指标→校验→预测）")
    ap.add_argument("-i", "--input", required=True, help="蛋白或 DNA FASTA")
    ap.add_argument("-o", "--output", required=True, help="输出目录")
    ap.add_argument("--purpose", choices=("vaccine", "nucleic_acid_therapy", "therapeutic_protein",
                                          "functional_assay", "structural", "reagent"))
    ap.add_argument("--glycosylation", action="store_true", help="依赖糖基化")
    ap.add_argument("--secreted", action="store_true", help="分泌蛋白")
    ap.add_argument("--disulfide", action="store_true", help="依赖二硫键")
    ap.add_argument("--toxic", action="store_true", help="对宿主有毒性")
    ap.add_argument("--stable-line", action="store_true", help="需要稳定株")
    ap.add_argument("--host", help="覆盖推荐宿主")
    ap.add_argument("--system", choices=VALID_SYSTEMS, help="覆盖推荐体系")
    ap.add_argument("--vector", help="覆盖推荐载体")
    ap.add_argument("--promoter", default="", help="启动子（原核）")
    ap.add_argument("--optimizer", choices=("auto", "codonbert", "rule"), default="auto")
    ap.add_argument("--mutations", default="", help="目标蛋白点突变，逗号分隔，如 M1T,C59A")
    ap.add_argument("--skip-predict", action="store_true", help="跳过表达/活性预测")
    args = ap.parse_args()

    import yaml
    kb = yaml.safe_load(KB_FILE.read_text(encoding="utf-8"))
    c2a, a2c, _ = codon_tables()
    _TARGET_REF = load_target_ref()

    out_dir = Path(args.output).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tmpdir = out_dir / "_tmp"
    tmpdir.mkdir(exist_ok=True)

    user = {
        "purpose": args.purpose,
        "requires_glycosylation": True if args.glycosylation else None,
        "secreted": True if args.secreted else None,
        "requires_disulfide": True if args.disulfide else None,
        "toxic": True if args.toxic else None,
        "stable_line": True if args.stable_line else None,
        "target_host": args.host,
    }

    report = {
        "workflow": "design_workflow", "version": "1.0",
        "input": str(Path(args.input).resolve()),
        "decision_table": "molecular_build_rules.yaml#expression_system_recommendation",
        "records": [],
    }

    final_fasta = out_dir / "final_optimized.fasta"
    with open(final_fasta, "w", encoding="utf-8") as fp_final:
        for idx, (header, raw) in enumerate(read_fasta(args.input)):
            name = header.split()[0] or f"seq{idx + 1}"
            protein = to_protein(raw, c2a)
            if not protein or set(protein) - set(AA):
                raise RuntimeError(f"{name}: 蛋白序列为空或含非法字符")
            features = extract_features(protein)
            rec = {"name": name, "protein_length": len(protein), "features": features}

            decision = recommend_system(features, user, kb)
            rec["decision"] = decision
            system = args.system or decision["recommendation"]["system"]
            vector = args.vector or decision["recommendation"]["vector"]
            host = args.host or decision["recommendation"]["host"] or "homo_sapiens"
            promoter = args.promoter or decision["recommendation"].get("promoter", "")
            host = {"HEK293": "homo_sapiens", "HEK293T": "homo_sapiens", "CHO": "cricetulus_griseus",
                    "BL21(DE3)": "escherichia_coli"}.get(host, host)
            ok, allowed = check_vector_whitelist(system, vector, kb)
            rec["selection"] = {"system": system, "vector": vector, "host": host,
                                "promoter": promoter, "vector_in_whitelist": ok,
                                "allowed_vectors": allowed}

            optimizer = choose_optimizer(args.optimizer, system, host)
            rec["optimizer"] = optimizer
            if optimizer == "codonbert":
                pf = tmpdir / f"{name}.protein.fasta"
                pf.write_text(f">{name}\n{protein}\n", encoding="utf-8")
                raw_out = tmpdir / f"{name}.raw.fasta"
                model_stdout = optimize_codonbert(pf, raw_out)
                try:
                    rec["codonbert_metrics"] = json.loads(model_stdout)
                except json.JSONDecodeError:
                    rec["codonbert_metrics"] = {"raw": model_stdout[-300:]}
                cds_records = read_fasta(raw_out)
            else:
                cds_records = [(name, reverse_translate(protein, host, a2c, kb))]

            cds = cds_records[0][1].upper().replace("U", "T")
            stops = {"TAA", "TAG", "TGA"}
            if cds[-3:] not in stops:
                cds += "TAA"
                rec["stop_codon_appended"] = True
            coding = cds[:-3]
            translated = translate(coding, c2a).rstrip("*")
            identity = round(sum(a == b for a, b in zip(translated, protein)) / len(protein) * 100, 2)
            gc = round((cds.count("G") + cds.count("C")) / len(cds) * 100, 2)
            rec["sequence_metrics"] = {
                "length_nt": len(cds), "gc_percent": gc, "stop_codon": cds[-3:],
                "protein_identity_pct": identity,
            }
            if identity < 100.0:
                diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(translated, protein)) if a != b]
                rec["sequence_metrics"]["identity_warning"] = (
                    "蛋白一致性非 100%：CodonBERT 生成式输出存在氨基酸改变，"
                    "该序列不是纯同义优化，用于实验前必须人工核对/回改。"
                )
                rec["sequence_metrics"]["aa_differences"] = [
                    {"position": p, "input": a, "output": b} for p, a, b in diffs
                ]
            fp_final.write(f">{name} system={system} vector={vector} host={host} optimizer={optimizer}\n{cds}\n")

            rec["quality"] = quality_metrics(cds, host, tmpdir, name)
            rec["validation"] = validate_construct(system, cds, vector, promoter, host, tmpdir, name)

            if not args.skip_predict:
                ident, auto_muts = auto_detect_target(protein)
                muts = parse_mutations(args.mutations) or auto_muts
                if muts:
                    rec["prediction"] = predict_target(muts)
                    rec["prediction"]["target_identity"] = ident
                    if args.mutations:
                        rec["prediction"]["mutation_source"] = "user"
                    else:
                        rec["prediction"]["mutation_source"] = "auto-diff-vs-target-ref"
                else:
                    rec["prediction"] = {
                        "available": False,
                        "reason": "非目标蛋白 或未提供突变；通用「序列→表达/活性」模型尚无训练数据，不编造预测值。",
                        "extension": "需按体系积累带标签数据（序列+实测表达+活性）后训练通用预测器；目标蛋白模型方法论可复用。",
                    }
            report["records"].append(rec)

    report_path = out_dir / "design_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    md = ["# 端到端分子设计报告", "", f"- 输入：`{report['input']}`",
          f"- 决策表：`{report['decision_table']}`", ""]
    for r in report["records"]:
        sel = r["selection"]
        md += [f"## {r['name']}", "",
               f"- 蛋白长度：{r['protein_length']} aa；MW ≈ {r['features']['molecular_weight_kda']} kDa",
               f"- 推荐规则：**{r['decision']['rule_id']} {r['decision']['rule_name']}**（置信度 {r['decision']['confidence']}）",
               f"- 选择：体系 `{sel['system']}` / 载体 `{sel['vector']}` / 宿主 `{sel['host']}` / 启动子 `{sel['promoter'] or '—'}`",
               f"- 载体白名单：{'通过' if sel['vector_in_whitelist'] else '未通过（需复核）'}",
               f"- 优化器：{r['optimizer']}；优化后 {r['sequence_metrics']['length_nt']} nt，GC {r['sequence_metrics']['gc_percent']}%，"
               f"终止 {r['sequence_metrics']['stop_codon']}，蛋白一致性 {r['sequence_metrics']['protein_identity_pct']}%"]
        if r["sequence_metrics"].get("identity_warning"):
            diffs = r["sequence_metrics"]["aa_differences"]
            md.append(f"- **⚠ {r['sequence_metrics']['identity_warning']}** 差异位点："
                      + ", ".join(f"{d['position']}{d['input']}→{d['output']}" for d in diffs))
        md += [f"- 理由：{r['decision']['rationale']}"]
        if r.get("decision", {}).get("notes"):
            md.append(f"- 备注：{r['decision']['notes']}")
        val = r["validation"]
        if val.get("available"):
            md.append(f"- SOP 校验：{val['summary']}")
        pred = r.get("prediction", {})
        if pred.get("available"):
            md.append("- 表达/活性预测（目标蛋白特异）：")
            for row in pred["rows"][:10]:
                md.append(f"  - {row['mutation']}: score {row['expression_score']}, "
                          f"结合 delta {row['binding_delta']}, {row['interpretation']}")
        else:
            md.append(f"- 表达/活性预测：不可用 — {pred.get('reason', '未运行')}")
        md.append("")
    md_path = out_dir / "design_report.md"
    md_path.write_text("\n".join(md), encoding="utf-8")

    print(json.dumps({
        "final_fasta": str(final_fasta), "report_json": str(report_path),
        "report_md": str(md_path), "n_records": len(report["records"]),
        "selections": [r["selection"] for r in report["records"]],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
