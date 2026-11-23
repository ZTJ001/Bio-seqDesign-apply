---
name: bio-seq-design-apply
description: "End-to-end molecular design for protein sequences: recommend expression system and vector, perform codon optimization, compute quality metrics, run construct validation, and predict expression/binding effects for a designated target protein. Use when a user provides a protein sequence and asks for codon optimization, expression-system/vector selection, expression or activity prediction, or a design report."
---

# Bio-seqDesign 应用技能（公开交付版）

端到端分子设计工具：**蛋白序列 → 体系/载体推荐 → 密码子优化 → 质量指标 → SOP 校验 → 表达与活性预测 → 报告**。

全部本地离线运行。不含任何实验原始数据、明细数据库或未脱敏知识库。

## 环境

```bash
cd <本仓库根目录>
python -m pip install -r requirements.txt
python src/smoke_test.py        # 一键自检，全部 PASS 即可用
```

如已装好依赖，可用 `BIO_SEQDESIGN_PYTHON=/path/to/python` 指定解释器。

## 推荐调用方式：agent API（稳定 JSON 接口）

所有语言/agent 统一通过 `src/agent_api.py` 调用。**stdout 只有一个 JSON 对象**，日志走 stderr。

```bash
python src/agent_api.py --sequence "MKVIFVLLL..." --purpose reagent --optimizer auto
```

或把请求写成 JSON 文件：

```json
{
  "sequence": "MKVIFVLLL...",
  "name": "my_protein",
  "purpose": "therapeutic_protein",
  "glycosylation": true,
  "secreted": true,
  "optimizer": "auto"
}
```

```bash
python src/agent_api.py --request req.json
```

返回结构（稳定契约，字段详见 `AGENT_CALL.md`）：

```json
{
  "status": "ok",
  "api_version": "1.0",
  "records": [{
    "name": "...",
    "recommendation": {"system": "...", "vector": "...", "host": "...", "rule_id": "...", "rationale": "..."},
    "optimization": {"optimizer": "...", "gc_percent": 64.8, "protein_identity_pct": 100.0, "identity_ok": true},
    "quality": {"mfe_kcal_per_mol": ..., "cai": ..., "validation_summary": {...}},
    "prediction": {"available": true/false, "reason": "...", "rows": [...]}
  }],
  "warnings": [...],
  "artifacts": {"final_fasta": "...", "report_json": "...", "report_md": "..."}
}
```

## 直接 CLI（生成人类可读报告）

```bash
python src/design_workflow.py -i examples/example_protein.fasta -o out_demo \
  --purpose reagent --optimizer auto
```

输出 `out_demo/` 下 `final_optimized.fasta`、`design_report.json`、`design_report.md`。

## 参数速查

| 参数 | 说明 |
|---|---|
| `--purpose` | vaccine / nucleic_acid_therapy / therapeutic_protein / functional_assay / structural / reagent |
| `--optimizer` | auto（真核 CodonBERT / 原核规则法）、codonbert、rule |
| `--glycosylation` `--secreted` `--disulfide` `--toxic` `--stable-line` | 蛋白特性，驱动体系推荐 |
| `--host` `--system` `--vector` | 覆盖推荐 |
| `--mutations` | 目标蛋白点突变，逗号分隔，如 `M1T,C59A` |

## 能力边界（务必如实告知用户）

1. **体系/载体推荐**是规则决策表，可解释但不是训练模型。
2. **表达/活性预测仅覆盖目标蛋白及其点突变**（酵母体系 DMS 训练）；其它蛋白 `prediction.available=false` 并给出原因，不得编造数值。
3. **蛋白一致性非 100% 是阻断项**（CodonBERT 生成式特性），`warnings` 会列出差异位点，须人工核对。
4. 预测仅计算辅助，不替代实验验证；详见 `USAGE_AGREEMENT.md`。
