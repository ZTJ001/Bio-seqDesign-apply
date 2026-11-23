# Bio-seqDesign-apply

端到端分子设计工具：**蛋白序列 → 表达体系/载体推荐 → 密码子优化 → 质量指标 → SOP 校验 → 表达与活性预测 → 报告**。

本地离线运行，纯 Python。可用于任意蛋白的密码子优化与体系/载体推荐；表达与活性预测目前仅覆盖目标蛋白及其点突变。

> **授权与合规**：使用本仓库即表示你已阅读并同意 [`USAGE_AGREEMENT.md`](USAGE_AGREEMENT.md) 与 [`LICENSE`](LICENSE)。
> 本仓库仅分发推理能力，不包含、也不会提供任何实验原始数据、明细数据库或未脱敏的知识库。

---

## 1. 快速开始（3 步）

```bash
# 1. 安装依赖（建议 Python 3.10+，用虚拟环境）
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 2. 一键自检：确认环境与两条优化路径可用
python src/smoke_test.py

# 3. 跑一个示例
python src/design_workflow.py -i examples/example_protein.fasta -o out_demo --purpose reagent
```

自检全部 `[PASS]` 并以 `全部通过` 结束即可。输出在 `out_demo/`：`final_optimized.fasta`、`design_report.json`、`design_report.md`。

依赖说明：`ViennaRNA` 仅用于 MFE 指标，缺失时该指标不可用但不阻断优化；其余依赖为推理与优化所必需。

如已装好依赖，可用环境变量指定解释器（所有子进程复用）：
```bash
export BIO_SEQDESIGN_PYTHON=/path/to/your/python
```

## 2. 给 AI Agent / 程序调用

`src/agent_api.py` 是稳定 JSON 接口：**stdout 仅输出一个 JSON 对象**，日志走 stderr，任何语言都能用 subprocess 调用。

```bash
# 直接给序列
python src/agent_api.py --sequence "MKVIFVLLL..." --purpose reagent

# 通过 stdin 传 JSON（跨平台）
echo '{"sequence":"MKVIFVLLL...","purpose":"reagent"}' | python src/agent_api.py

# JSON 请求文件
python src/agent_api.py --request req.json
```

```python
import json, subprocess, sys
p = subprocess.run([sys.executable, "src/agent_api.py"],
                   input=json.dumps({"sequence": "MKVIFVLLL...", "purpose": "reagent"}),
                   text=True, capture_output=True)
result = json.loads(p.stdout)     # {"status":"ok", "records":[...], "warnings":[...], "artifacts":{...}}
```

完整字段契约见 [`AGENT_CALL.md`](AGENT_CALL.md)。WorkBuddy agent 可直接装载本仓库的 [`SKILL.md`](SKILL.md)。

## 3. 能力边界（重要，请如实告知用户）

1. **体系/载体推荐**是规则决策表，按蛋白长度、糖基化、二硫键、分泌信号、毒性、用途给出建议，可解释但不是训练模型。
2. **表达/活性预测仅覆盖目标蛋白及其点突变**，训练体系为酵母。其他蛋白返回 `prediction.available=false` 并说明原因，**不会编造数值**。
   - `effect_predictor`：突变对酵母表达 score 的影响（测试集 Pearson r≈0.91）
   - `binding_predictor`：突变对受体结合的影响（Pearson r≈0.85）
3. **蛋白一致性非 100% 是阻断项**：CodonBERT 是生成式模型，可能改变个别氨基酸；`warnings` 会列出差异位点，实验前须人工核对。
4. 所有预测仅为计算辅助，不替代实验验证；禁止反编译模型或反向推断训练数据（见使用协议）。

## 4. 目录结构

```
.
├── src/
│   ├── design_workflow.py     # 端到端工作流 CLI
│   ├── agent_api.py           # 稳定 JSON 调用接口（推荐给 agent）
│   ├── smoke_test.py          # 自检
│   ├── check_engine.py        # SOP 校验引擎
│   ├── optimize_model.py      # CodonBERT 优化封装
│   ├── featurizer.py          # 预测模型特征构造
│   ├── viennarna_adapter.py / jcat_adapter.py   # MFE / CAI
│   └── codonbert/             # 第三方 CodonBERT 推理源码（见其 LICENSE）
├── config/target_reference.json # 目标蛋白参考序列（仅序列）
├── models/                      # 3 个推理权重
├── examples/                  # 示例输入
├── molecular_build_rules.yaml # 脱敏规则库（仅运行时规则）
├── SKILL.md / AGENT_CALL.md   # agent 集成
└── LICENSE / USAGE_AGREEMENT.md
```

第三方组件 CodonBERT 的原始许可证见 `src/codonbert/LICENSE`。
