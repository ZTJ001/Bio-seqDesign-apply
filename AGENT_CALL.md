# Agent 调用说明

本仓库提供两层调用接口，任何 agent / 程序均可离线调用。**所有能力都在本地运行，不联网、不上传序列、不依赖任何外部数据库。**

## 接口一：agent API（推荐，稳定 JSON 契约）

入口：`src/agent_api.py`。用当前环境 Python 调用，`stdout` 仅输出一个 JSON 对象，便于程序解析；日志与错误走 `stderr`。

### 请求字段

| 字段 | 必填 | 说明 |
|---|---|---|
| `sequence` / `fasta` / `input_file` | 三选一 | 蛋白序列（裸序列）/ FASTA 文本 / FASTA 文件路径 |
| `name` | 否 | 序列名 |
| `purpose` | 否 | vaccine / nucleic_acid_therapy / therapeutic_protein / functional_assay / structural / reagent |
| `optimizer` | 否 | auto（默认）/ codonbert / rule |
| `glycosylation` `secreted` `disulfide` `toxic` `stable_line` | 否 | 布尔，蛋白特性，驱动体系推荐 |
| `host` `system` `vector` `promoter` | 否 | 手动覆盖推荐 |
| `mutations` | 否 | 目标蛋白点突变，逗号分隔，如 `M1T,C59A` |
| `skip_predict` | 否 | 布尔，跳过预测 |
| `output_dir` | 否 | 产物目录；缺省自动创建临时目录 |

### 三种传入方式

```bash
# A. 直接给序列
python src/agent_api.py --sequence "MKVIFVLLL..." --purpose reagent

# B. JSON 文件
python src/agent_api.py --request req.json

# C. stdin 管道（跨平台，无需 --request）
echo '{"sequence":"MKVIFVLLL..."}' | python src/agent_api.py
```

### 返回字段（稳定，底层变化不影响调用方）

```
status           "ok" | "error"
api_version      "1.0"
input            {name, n_records}
records[]        {
  name, protein_length,
  recommendation  {system, vector, host, promoter, vector_in_whitelist,
                   rule_id, rule_name, confidence, rationale},
  optimization    {optimizer, length_nt, gc_percent, stop_codon,
                   protein_identity_pct, identity_ok, aa_differences[]},
  quality         {mfe_kcal_per_mol, cai, validation_available, validation_summary},
  prediction      {available, reason, note, model[], rows[]}
}
warnings[]       需人工注意的告警（如蛋白一致性非 100%）
artifacts        {final_fasta, report_json, report_md}
```

出错时返回 `{"status":"error", "error":{"type":..., "message":...}}` 并以非 0 退出码结束。

### 调用示例（Python）

```python
import json, subprocess, sys
req = {"sequence": "MKVIFVLLL...", "purpose": "reagent"}
# 通过 stdin 传入 JSON（跨平台，Windows/Linux/macOS 通用）
p = subprocess.run([sys.executable, "src/agent_api.py"],
                   input=json.dumps(req), text=True, capture_output=True)
result = json.loads(p.stdout)
assert result["status"] == "ok"
for rec in result["records"]:
    print(rec["recommendation"]["system"], rec["optimization"]["gc_percent"])
```

## 接口二：直接 CLI（人类可读报告）

```bash
python src/design_workflow.py -i examples/example_protein.fasta -o out_demo --purpose reagent
```

生成 `final_optimized.fasta` + `design_report.json` + `design_report.md`。

## WorkBuddy 技能装载

仓库根的 `SKILL.md` 是 WorkBuddy 技能描述。把本仓库作为技能目录装载后，agent 会自动获知上述能力与边界。

## 边界（务必遵守）

- 表达/活性预测**仅对目标蛋白及其点突变有效**；其它蛋白返回 `available=false`，不得编造数值。
- 体系/载体推荐为规则决策表，非训练模型。
- 蛋白一致性非 100% 为阻断项，见 `warnings`，须人工核对。
- 禁止反编译模型或反向推断训练数据，见 `USAGE_AGREEMENT.md`。
