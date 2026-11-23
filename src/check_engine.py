#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分子构建校验引擎 check_engine.py
=================================
依据知识库: molecular_build_rules.yaml (v1.0, 2026-11-05)
支持系统  : mrna / mammalian / prokaryotic
校验分级  : ERROR(阻断,必须修改) / WARNING(风险复核) / INFO(提示备忘)

用法示例:
    python3 check_engine.py mrna  test_mrna_seq.fasta
    python3 check_engine.py mammalian  input.fasta  --vector pcDNA3.1 --host HEK293
    python3 check_engine.py prokaryotic input.fasta  --vector pET-28a(+) --host BL21(DE3)

输入格式: FASTA, 支持两种模式
  1) 注释分段模式(测试用, 推荐):
     >name
     5UTR:...GCGCGC...
     CDS:...AUGGCU...
     3UTR:...AUUUA...
     polyA:...AAAA...
  2) 整段模式(真实输入):
     >name
     GGGGCCACC AUGGCU...UAAUAA AUUUAAA...AAAA
     (引擎自动定位: 5'UTR=首个AUG前, CDS=AUG..终止密码子, 3'UTR=终止后)
"""
import sys
import re
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("[FATAL] 缺少 PyYAML, 请先安装: python3 -m pip install PyYAML")

# ---------- 路径 ----------
BASE_DIR = Path(__file__).resolve().parent
# 规则库与 check_engine.py 同目录（仓库嵌套结构下 BASE_DIR.parent 会指向外层目录，
# 因此同时兼容同目录与父目录两种布局，自动选取实际存在的那个）
RULES_FILE = BASE_DIR / "molecular_build_rules.yaml"
if not RULES_FILE.exists():
    RULES_FILE = BASE_DIR.parent / "molecular_build_rules.yaml"

# ---------- 基础常量 ----------
START_CODONS = {"AUG", "ATG"}
STOP_CODONS = {"UAA", "UAG", "UGA"}
# Type IIS 酶识别位点 (mrna SOP 要求识别位点与切割位点不一致)
TYPE_IIS_SITES = {
    "BsaI": ("GGTCTC", "GAGACC"),   # 正向识别 / 反向互补
    "BsmBI": ("CGTCTC", "GAGACG"),
    "BbsI": ("GAAGAC", "GTCTTC"),
}
AU_RICH_MOTIF = "AUUUA"             # 3'UTR 不稳定元件
POLYA_SIGNAL = "AAUAAA"             # poly(A) 信号


def load_rules():
    """加载 molecular_build_rules.yaml 规则库"""
    if not RULES_FILE.exists():
        sys.exit(f"[FATAL] 规则库不存在: {RULES_FILE}")
    with open(RULES_FILE, encoding="utf-8") as f:
        return yaml.safe_load(f)


def parse_input(text):
    """解析输入: 返回 {'5utr': str, 'cds': str, '3utr': str, 'polyA': str, 'full': str}"""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip() and not ln.startswith(">")]
    body = "".join(lines)

    seg = {"5utr": "", "cds": "", "3utr": "", "polyA": "", "full": body}

    # 模式1: 注释分段模式 (按标签定位, 每段到下一个标签前结束, 避免误吞标签名)
    label_re = re.compile(r"(5UTR|CDS|3UTR|polyA)\s*[:：]\s*")
    matches = list(label_re.finditer(body))
    if matches and len(matches) >= 3:
        for idx, m in enumerate(matches):
            key = m.group(1).upper()
            start = m.end()
            end = matches[idx + 1].start() if idx + 1 < len(matches) else len(body)
            val = re.sub(r"(5UTR|CDS|3UTR|polyA)\s*[:：]\s*$", "", body[start:end]).strip().upper()
            if key == "5UTR":
                seg["5utr"] = val
            elif key == "CDS":
                seg["cds"] = val
            elif key == "3UTR":
                seg["3utr"] = val
            elif key == "POLYA":  # polyA 大写化后为 POLYA
                seg["polyA"] = val
        seg["full"] = seg["5utr"] + seg["cds"] + seg["3utr"] + seg["polyA"]
        if seg["cds"]:
            return seg

    # 模式2: 自动定位
    seq = body.upper()
    start = seq.find("AUG")
    if start < 0:
        start = seq.find("ATG")
    if start < 0:
        return seg  # 找不到起始密码子, 整段按CDS候选
    seg["5utr"] = seq[:start]
    # CDS: 从AUG起, 找框内终止密码子
    cds_end = len(seq)
    for i in range(start, len(seq) - 2, 3):
        if seq[i:i+3] in STOP_CODONS:
            cds_end = i
            break
    seg["cds"] = seq[start:cds_end]
    tail = seq[cds_end:]
    # polyA: 末尾连续 A
    m = re.search(r"A+$", tail)
    if m and len(m.group(0)) >= 5:
        seg["polyA"] = m.group(0)
        tail = tail[:m.start()]
    seg["3utr"] = tail
    return seg


def gc_content(s):
    """计算 GC 百分比"""
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s) * 100.0


def revcomp(s):
    """反向互补 (DNA/RNA 通用, U 视作 T 的配对)"""
    comp = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(s.upper()))


# ==================== 校验器 ====================
class Validator:
    """按规则库分级校验, 产出报告条目"""

    def __init__(self, rules):
        self.rules = rules
        self.items = []

    def add(self, level, rule_id, pos, detail, suggestion):
        self.items.append({
            "level": level, "rule": rule_id, "pos": pos,
            "detail": detail, "suggestion": suggestion,
        })

    def _system(self, sys_name):
        return self.rules["system_list"][sys_name]

    # ---- 5'UTR ----
    def check_five_utr(self, sys_name, seg):
        """error: 不含AUG; 末尾精确匹配Kozak GCCACC; G含量尽量少"""
        sys_cfg = self._system(sys_name)["sequence_design"]
        five_cfg = sys_cfg.get("five_utr", {})
        if not five_cfg:
            return
        utr = seg["5utr"]
        pos = "5'UTR (1-{}nt)".format(len(utr)) if utr else "5'UTR(空)"
        # 1) 是否含起始密码子
        hit = [m.start() for m in re.finditer(r"AUG|ATG", utr)]
        if hit:
            self.add("ERROR", "five_utr.no_start_codon",
                     pos + " AUG@" + ",".join(str(h + 1) for h in hit),
                     "5'UTR 内检测到起始密码子 AUG/ATG, 可能导致翻译起始位点提前/错乱",
                     "删除 5'UTR 内 AUG/ATG, 或调整序列使首个 AUG 成为真实翻译起点")
        else:
            self.add("INFO", "five_utr.no_start_codon", pos,
                     "5'UTR 未含起始密码子, 符合要求", "无需处理")
        # 2) Kozak 序列
        kozak_len = 6
        tail = utr[-kozak_len:] if len(utr) >= kozak_len else utr
        if tail == "GCCACC":
            self.add("INFO", "five_utr.kozak", pos + " 末端6nt=" + tail,
                     "Kozak 序列 GCCACC 已精确匹配", "无需处理")
        else:
            self.add("ERROR", "five_utr.kozak", pos + " 末端6nt=" + (tail or "空"),
                     "5'UTR 末端未匹配 Kozak 序列 GCCACC",
                     "在 5'UTR 末端(紧邻 CDS 起始 AUG 上游)添加 GCCACC")
        # 3) G 含量 (仅 mrna 规则库含此检查提示)
        if sys_name == "mrna" and utr:
            g_pct = utr.count("G") / len(utr) * 100
            if g_pct > 40:
                self.add("WARNING", "five_utr.gc", pos + f" G%={g_pct:.1f}",
                         "5'UTR G 含量偏高, 可能形成稳定二级结构阻碍翻译",
                         "替换 G 为 A/U, 降低 5'UTR G 含量")

    # ---- CDS ----
    def check_cds(self, sys_name, seg):
        sys_cfg = self._system(sys_name)["sequence_design"]
        cds_cfg = sys_cfg.get("cds", {})
        cds = seg["cds"]
        pos = "CDS"
        if not cds:
            self.add("ERROR", "cds.missing", pos, "未识别到 CDS 序列", "补充编码区序列")
            return
        # 起始密码子
        if cds[:3] in START_CODONS:
            self.add("INFO", "cds.start_codon", pos + " 起始3nt=" + cds[:3],
                     "起始密码子 {} 符合偏好 AUG".format(cds[:3]), "无需处理")
        else:
            self.add("ERROR", "cds.start_codon", pos + " 起始3nt=" + cds[:3],
                     "CDS 不以起始密码子 AUG 开头",
                     "确保 CDS 第一个密码子为 AUG")
        # 终止密码子 + 双终止
        end3 = cds[-3:]
        end6 = cds[-6:]
        if end3 in STOP_CODONS:
            note = ""
            if end6[:3] == end3 and end3 == "UAA":
                note = " (已使用双终止 UAAUAA, 推荐)"
            elif end6[:3] in STOP_CODONS:
                note = " (存在双终止密码子, 推荐)"
            self.add("INFO", "cds.stop_codon", pos + " 终止3nt=" + end3,
                     "终止密码子 {} 符合偏好 UAA{}".format(end3, note),
                     "若为单终止, 建议在末尾追加一个 UAA 形成双终止")
        else:
            self.add("WARNING", "cds.stop_codon", pos + " 末尾3nt=" + end3,
                     "CDS 末尾未识别到终止密码子",
                     "在 CDS 末尾添加 UAA(推荐双终止 UAAUAA)")
        # 框内提前终止检查 (保留末尾6nt双终止区不检查)
        for i in range(3, max(3, len(cds) - 6), 3):
            if cds[i:i+3] in STOP_CODONS:
                self.add("ERROR", "cds.internal_stop", pos + " 框内终止@" + str(i + 1),
                         "CDS 内部存在提前终止密码子 {}".format(cds[i:i+3]),
                         "修正序列消除内部终止密码子(密码子优化时规避)")
                break
        # GC 范围 (mammalian/prokaryotic 为 error)
        gcr = cds_cfg.get("gc_range")
        if gcr:
            gc = gc_content(cds)
            lo, hi = gcr["min"], gcr["max"]
            if lo <= gc <= hi:
                self.add("INFO", "cds.gc_range", pos + f" GC%={gc:.1f}",
                         f"CDS GC 含量 {gc:.1f}% 在 {lo}-{hi}% 范围内", "无需处理")
            else:
                self.add("ERROR", "cds.gc_range", pos + f" GC%={gc:.1f}",
                         f"CDS GC 含量 {gc:.1f}% 超出 {lo}-{hi}%",
                         "进行密码子优化(ai_hook=codon_optimize), 调整 GC 至范围内")
        # 密码子优化提示 (info + ai_hook)
        if "codon_optimization" in cds_cfg:
            self.add("INFO", "cds.codon_optimization", pos,
                     "建议基于宿主进行密码子优化(Jcat/常用密码子), ai_hook=codon_optimize",
                     "合成前确认密码子适应指数(CAI)")
        # 二级结构 (AUG 下游 10-15bp)
        if "secondary_structure" in cds_cfg:
            win = cds[3:3 + 15] if len(cds) >= 18 else cds[3:]
            gcw = gc_content(win)
            self.add("WARNING", "cds.secondary_structure", pos + " AUG下游15nt窗口",
                     f"AUG 下游窗口 GC%={gcw:.1f}; 需用 ViennaRNA 计算 ΔG (ai_hook=rna_fold_predict)",
                     "若 ΔG < -10 kcal/mol 则优化该窗口密码子, 避免稳定茎环")
        # Type IIS 酶切位点 (mrna); 序列U→T后按DNA形式搜索酶识别位点
        if "type_iis_restriction" in sys_cfg:
            cds_dna = cds.replace("U", "T")
            for enz, (fwd, rev) in TYPE_IIS_SITES.items():
                for mm in re.finditer(fwd, cds_dna):
                    self.add("ERROR", "cds.type_iis", pos + " {}识别位点@{}".format(enz, mm.start() + 1),
                             "CDS 内存在 {} 识别位点 {}; SOP要求识别位点与切割位点不一致, 需人工复核切割位置".format(enz, fwd),
                             "若切割位点落入编码区或与识别位点冲突, 请移除/突变该位点")
        # 信号肽 (mammalian warning)
        if "signal_peptide" in cds_cfg and cds_cfg["signal_peptide"].get("level") == "warning":
            self.add("WARNING", "cds.signal_peptide", pos + " N端20aa",
                     "分泌型/膜蛋白必须有信号肽; 胞内蛋白不应包含信号肽 — 请人工确认蛋白类型",
                     "按目标蛋白定位决定是否保留 N 端信号肽")
        # 共表达 (mammalian warning)
        if "co_expression" in cds_cfg:
            self.add("WARNING", "cds.co_expression", pos,
                     "如需要共表达: P2A会残留碱基, IRES下游表达量低 — 请确认是否使用及选择",
                     "需要共表达时选择合适策略; 不需要则忽略")

    # ---- 3'UTR / polyA ----
    def check_three_utr_polya(self, sys_name, seg):
        sys_cfg = self._system(sys_name)["sequence_design"]
        utr3 = seg["3utr"]
        pos = "3'UTR (1-{}nt)".format(len(utr3)) if utr3 else "3'UTR(空)"
        if sys_name == "mrna":
            cfg = sys_cfg.get("three_utr", {})
        else:
            cfg = sys_cfg.get("three_utr_polyA", {})
        if not cfg:
            return
        if len(utr3) > 500:
            self.add("WARNING", "three_utr.length", pos,
                     f"3'UTR 长度 {len(utr3)}nt 偏长, 可能引入过多调控元件",
                     "精简 3'UTR 长度, 去除冗余序列")
        au = [m.start() for m in re.finditer(AU_RICH_MOTIF, utr3)]
        if au:
            self.add("WARNING", "three_utr.au_rich", pos + " AUUUA@" + ",".join(str(a + 1) for a in au),
                     "3'UTR 含 AU-rich 不稳定元件 AUUUA x{} 个, 可能降低 mRNA/蛋白稳定性".format(len(au)),
                     "移除或替换 AU-rich 基序")
        else:
            self.add("INFO", "three_utr.au_rich", pos, "未检测到 AU-rich 不稳定元件", "无需处理")
        self.add("INFO", "three_utr.mirna", pos,
                 "miRNA 结合位点预测需专用数据库/工具(ai_hook), 位点宜少、序列不宜过长",
                 "交付前做 miRNA 结合位点统计, 位点多则精简")
        # polyA 信号
        if POLYA_SIGNAL in utr3:
            self.add("INFO", "three_utr.polya_signal", pos + " AAUAAA@" + str(utr3.find(POLYA_SIGNAL) + 1),
                     "检测到 poly(A) 信号 AAUAAA", "无需处理")
        else:
            self.add("WARNING", "three_utr.polya_signal", pos,
                     "3'UTR 未检测到典型 poly(A) 信号 AAUAAA",
                     "确认是否需要添加 poly(A) 信号(哺乳/通用)")
        # polyA 尾巴
        pa = seg.get("polyA", "")
        if pa:
            self.add("INFO", "polyA.tail", "polyA 尾部",
                     "polyA 长度 {}nt; 测序比对时监控 A 长度是否丢失".format(len(pa)),
                     "测序时重点观察 A 尾巴完整性")
        else:
            self.add("WARNING", "polyA.tail", "3'端",
                     "未检测到 polyA 尾巴(连续A≥5nt)",
                     "若为 mRNA 模板应保留 polyA; 若为质粒模板可忽略")

    # ---- 5'加帽 (mrna) ----
    def check_cap(self, sys_name, seg):
        sys_cfg = self._system(sys_name)["sequence_design"]
        cap = sys_cfg.get("five_cap")
        if not cap:
            return
        utr = seg["5utr"]
        agg = [m.start() for m in re.finditer(r"AGG", utr)]
        if agg:
            self.add("WARNING", "five_cap.agg", "5'UTR AGG@" + ",".join(str(a + 1) for a in agg),
                     "共转录加帽需人工复核 5' 端 AGG 起始位点(加帽效率)",
                     "确认加帽方式(共转录/酶法)与 AGG 位点适配")
        else:
            self.add("INFO", "five_cap.agg", "5'UTR",
                     "5'UTR 未含 AGG; 共转录加帽时注意 AGG 起始位点选择", "按加帽方案复核")

    # ---- SD/RBS (prokaryotic) ----
    def check_sd_rbs(self, sys_name, seg):
        sys_cfg = self._system(sys_name)["sequence_design"]
        sd = sys_cfg.get("sd_rbs")
        if not sd:
            return
        utr = seg["5utr"]
        cds = seg["cds"]
        # AUG 上游 7nt 区域应含 AAGGAG
        upstream = utr[-13:] if len(utr) >= 13 else utr
        preferred = sd["preferred"]
        if preferred in upstream:
            pos7 = upstream.rfind(preferred)
            dist = len(utr) - (len(upstream) - pos7)  # 到AUG的距离
            self.add("INFO", "sd_rbs", "5'UTR RBS区",
                     f"RBS {preferred} 位于 AUG 上游约 {dist}nt (SOP偏好7nt)",
                     "若距离偏离7nt较多, 微调间隔序列")
        else:
            self.add("ERROR", "sd_rbs", "5'UTR RBS区",
                     f"未在 AUG 上游找到优选 RBS {preferred}(应位于上游7nt)",
                     "在 AUG 上游约7nt处添加 AAGGAG")

    # ---- 启动子 (prokaryotic) ----
    def check_promoter(self, sys_name, promoter, host_note):
        sys_cfg = self._system(sys_name)["sequence_design"]
        pc = sys_cfg.get("promoter")
        if not pc:
            return
        high = pc.get("high_expression", [])
        if promoter in high:
            self.add("WARNING", "promoter.high", "启动子",
                     f"启动子 {promoter} 为高表达启动子, 可能导致包涵体或代谢负担过重",
                     "若蛋白有毒/难折叠, 改用诱导型(lac/araBAD)或低表达(trc)启动子")
        ind = pc.get("inducible", [])
        if promoter in ind:
            self.add("INFO", "promoter.inducible", "启动子",
                     f"启动子 {promoter} 为诱导型, 可通过IPTG/阿拉伯糖调控表达时间, 适合毒性蛋白",
                     "确认诱导剂与诱导时机")
        low = pc.get("low_expression", [])
        if promoter in low:
            self.add("INFO", "promoter.low", "启动子",
                     f"启动子 {promoter} 为低表达, 适合毒性蛋白或需要正确折叠的蛋白",
                     "确认表达量是否满足需求")

    # ---- 通用: 命名/载体/模板 ----
    def check_meta(self, sys_name, seq_name, vector):
        sys_cfg = self._system(sys_name)["sequence_design"]
        self.add("INFO", "meta.naming", "序列头部",
                 f"命名需'编号+全称'完整展示 (当前序列名: {seq_name or '未命名'})",
                 "完善序列/载体注释: 同时含编号与全称")
        allowed = sys_cfg.get("backbone", {}).get("allowed_vector", [])
        if vector:
            if vector in allowed:
                self.add("INFO", "meta.backbone", "载体",
                         f"载体 {vector} 在允许列表 {allowed} 内", "无需处理")
            else:
                self.add("WARNING", "meta.backbone", "载体",
                         f"载体 {vector} 不在 SOP 允许列表 {allowed} 内",
                         "确认载体选择或补充 SOP 允许清单")
        else:
            self.add("INFO", "meta.backbone", "载体",
                     f"未指定载体; SOP允许: {allowed}",
                     "构建前确认目标载体")

    # ---- 汇总 ----
    def validate(self, sys_name, seg, seq_name="", vector="", promoter="", host_note=""):
        sys_cfg = self._system(sys_name)
        self.check_meta(sys_name, seq_name, vector)
        self.check_cap(sys_name, seg)
        self.check_five_utr(sys_name, seg)
        self.check_cds(sys_name, seg)
        self.check_three_utr_polya(sys_name, seg)
        if sys_name == "prokaryotic":
            self.check_sd_rbs(sys_name, seg)
            self.check_promoter(sys_name, promoter, host_note)
        return self.items


# ==================== 密码子偏好打分与优化 ====================
def _resolve_host(system_name, host, rules):
    """按系统解析宿主: 默认映射 / --host 覆盖(支持细胞名别名)"""
    co = rules["codon_optimization"]
    alias = co.get("host_alias", {})
    if host:
        key = alias.get(host, host)
        if key in rules["codon_usage_bias"] and "preferred_codons" in rules["codon_usage_bias"][key]:
            return key
        raise ValueError("未知宿主: {} (可用别名: {})".format(host, ",".join(sorted(alias))))
    cfg = co["system_host_map"].get(system_name)
    if not cfg:
        raise ValueError("系统 {} 未配置宿主映射".format(system_name))
    return cfg["default"]


def _codon_maps(rules):
    """从遗传密码表构建 codon->aa 反向映射 (含终止 *)"""
    gc = rules["codon_optimization"]["genetic_code"]
    codon2aa = {}
    for aa, codons in gc.items():
        if aa.startswith("_") or not isinstance(codons, list):
            continue
        for c in codons:
            codon2aa[c] = aa
    for c in gc["_stop"]:
        codon2aa[c] = "*"
    return codon2aa


def _as_dna(seq):
    return seq.upper().replace("U", "T")


def _split_codons(cds):
    """按读框切为三联体; 末尾不足3nt的残基也保留, 保证优化前后长度守恒"""
    return [cds[i:i + 3] for i in range(0, len(cds), 3)]


def _translate(cds, codon2aa):
    """翻译 CDS 为氨基酸串 (终止后截断, 用于校验优化不改蛋白)"""
    aa_seq = []
    for codon in _split_codons(cds):
        aa = codon2aa.get(codon, "?")
        if aa == "*":
            break
        aa_seq.append(aa)
    return "".join(aa_seq)


def score_codon_preference(system_name, seq_text, host=None, rules=None):
    """对输入序列的 CDS 按宿主偏好密码子打分 (偏好使用率, 非严格 CAI)。

    返回 dict: host/score/level/matched/comparable/nonpreferred_detail 等。
    """
    rules = rules or load_rules()
    if system_name not in ("mrna", "mammalian", "prokaryotic"):
        raise ValueError("系统必须为 mrna / mammalian / prokaryotic")
    resolved = _resolve_host(system_name, host, rules)
    pref = rules["codon_usage_bias"][resolved]["preferred_codons"]
    codon2aa = _codon_maps(rules)
    seg = parse_input(seq_text)
    dna = _as_dna(seg["cds"])
    codons = _split_codons(dna)
    matched = comparable = 0
    detail = []
    for i, codon in enumerate(codons):
        aa = codon2aa.get(codon)
        if aa is None or aa == "*" or i == 0:  # 未知/终止/起始不计入
            continue
        target = pref.get(aa)
        if not target:
            continue
        comparable += 1
        if codon == target:
            matched += 1
        else:
            detail.append({"pos": i * 3 + 1, "codon": codon, "aa": aa, "preferred": target})
    score = round(matched / comparable, 4) if comparable else None
    if score is None:
        level = "n/a"
    elif score >= 0.80:
        level = "excellent"
    elif score >= 0.60:
        level = "good"
    elif score >= 0.40:
        level = "fair"
    else:
        level = "poor"
    return {
        "system": system_name,
        "host": resolved,
        "host_alias_used": bool(host),
        "score": score,
        "level": level,
        "matched": matched,
        "comparable": comparable,
        "total_codons": len(codons),
        "gc_percent": round(gc_content(dna), 2),
        "nonpreferred_count": len(detail),
        "nonpreferred_detail": detail[:50],
    }


def optimize_sequence(system_name, seq_text, host=None, rules=None):
    """按宿主偏好密码子对 CDS 同义替换优化, 保持氨基酸序列/UTR/polyA 不变。

    返回 dict, 含 optimized_cds / optimized_sequence / changes / before / after。
    """
    rules = rules or load_rules()
    if system_name not in ("mrna", "mammalian", "prokaryotic"):
        raise ValueError("系统必须为 mrna / mammalian / prokaryotic")
    resolved = _resolve_host(system_name, host, rules)
    pref = rules["codon_usage_bias"][resolved]["preferred_codons"]
    stop_pref = rules["codon_usage_bias"][resolved].get("stop_codon", {}).get("preferred", "TAA")
    codon2aa = _codon_maps(rules)
    seg = parse_input(seq_text)
    dna = _as_dna(seg["cds"])
    codons = _split_codons(dna)
    opt = []
    changes = []
    for i, codon in enumerate(codons):
        aa = codon2aa.get(codon)
        if aa is None:
            opt.append(codon)
            continue
        if i == 0 and aa == "M":  # 起始密码子保留
            opt.append(codon)
            continue
        if aa == "*":  # 终止: 原核统一首选 TAA, 真核保留原终止
            new = stop_pref if system_name == "prokaryotic" else codon
            opt.append(new)
            if new != codon:
                changes.append({"pos": i * 3 + 1, "from": codon, "to": new, "aa": "终止"})
            continue
        target = pref.get(aa)
        if target and target != codon:
            changes.append({"pos": i * 3 + 1, "from": codon, "to": target, "aa": aa})
            opt.append(target)
        else:
            opt.append(codon)
    optimized_cds_dna = "".join(opt)
    before = score_codon_preference(system_name, seq_text, host=host, rules=rules)
    after_dna = score_codon_preference(system_name, "CDS:" + optimized_cds_dna, host=host, rules=rules)
    protein_before = _translate(dna, codon2aa)
    protein_after = _translate(optimized_cds_dna, codon2aa)
    # 保持输入 RNA/DNA 书写风格
    as_rna = "U" in seg["full"] and "T" not in seg["full"]
    def style(s):
        return s.replace("T", "U") if as_rna else s
    optimized_cds = style(optimized_cds_dna)
    optimized_sequence = (
        style(_as_dna(seg["5utr"])) + optimized_cds +
        style(_as_dna(seg["3utr"])) + style(_as_dna(seg["polyA"]))
    )
    return {
        "system": system_name,
        "host": resolved,
        "strategy": rules["codon_optimization"]["optimization"]["strategy"],
        "original_cds": style(dna),
        "optimized_cds": optimized_cds,
        "optimized_sequence": optimized_sequence,
        "protein": protein_after[:80] + ("..." if len(protein_after) > 80 else ""),
        "protein_unchanged": protein_before == protein_after,
        "changes_count": len(changes),
        "changes": changes[:100],
        "score_before": before["score"],
        "score_after": after_dna["score"],
        "gc_before": before["gc_percent"],
        "gc_after": after_dna["gc_percent"],
        "cautions": rules["codon_optimization"]["optimization"]["cautions"],
    }


# ==================== 引物设计 ====================
def design_primers(seg, homology_arm=20):
    """同源重组引物: 5'端为同源臂(15-25bp), 3'端为插入片段特异性序列"""
    if not (15 <= homology_arm <= 25):
        homology_arm = 20
    ins = seg["full"]
    f_primer = seg["5utr"][:20] or ins[:20]       # 上游扩增段
    # 下游: 序列末端20nt 反向互补
    tail = (seg["3utr"] or ins)[-20:] if len(seg["full"]) >= 20 else seg["full"]
    r_bind = revcomp(seg["full"][-20:])
    primers = {
        "homology_arm_len": homology_arm,
        "note": "同源臂序列需按实际载体末端序列替换(此处以插入片段端占位演示)",
        "fwd": "[同源臂20bp] + " + f_primer,
        "rev": "[同源臂20bp(反向互补)] + " + r_bind,
        "fwd_seq": f_primer,
        "rev_seq": r_bind,
        "conc": "10 μM, 引物与样本一一对应",
        "plate": "8连排按顺序摆放并及时编号",
    }
    return primers


# ==================== SOP 备忘 ====================
def sop_memo(system_name, rules):
    cp = rules["common_protocol"]
    lines = []
    lines.append("### 实验 SOP 备忘 ({} / common_protocol)".format(system_name))
    lines.append("**PCR:** " + "; ".join(
        "{} ({})".format(m["name"], m["speed"]) for m in cp["pcr"]["mastermix_option"]))
    lines.append("   - " + cp["pcr"].get("note", "") + " (最后加酶)")
    lines.append("**酶切:** " + cp["restriction_digest"]["buffer_note"] + "; " +
                 cp["restriction_digest"]["enzyme_note"] + "; " +
                 cp["restriction_digest"]["instrument"])
    lines.append("**胶回收/纯化:** " + cp["gel_extraction"]["comb"] + ", 洗脱 " +
                 cp["gel_extraction"]["elution_volume"])
    lig = cp["ligation"]
    lines.append("**连接(同源重组):** " + lig["homologous_recombination"]["mix"] +
                 ", 模板量" + lig["homologous_recombination"]["min_input"] +
                 ", " + lig["homologous_recombination"]["instrument"])
    lines.append("**连接(T4):** " + lig["t4_ligase"]["enzyme"] + ", " +
                 lig["t4_ligase"]["instrument"])
    tf = cp["transformation"]
    lines.append("**转化感受态(通用):** " + "/".join(tf["competent_general"]["strain_list"]) +
                 " — " + tf["competent_general"]["note"])
    if "prokaryotic_special_strain" in tf:
        lines.append("**表达菌株(原核):** BL21(DE3)=常规/T7, Rosetta=稀有密码子, " +
                     "BL21-CodonPlus=稀有tRNA, Origami=二硫键, AD494=周质分泌; " +
                     "毒性蛋白用低拷贝载体+严谨/温控启动子")
    cs = cp["colony_screen_verify"]
    lines.append("**挑菌:** " + cs["pick_colony"]["rule"])
    lines.append("**菌落PCR:** " + cs["colony_pcr"]["plate_8well"] + "; " +
                 cs["colony_pcr"]["mix_prep"] + "; 体系 " + cs["colony_pcr"]["system_volume"])
    lines.append("**测序:** 每构建提交 {} 个阳性克隆".format(cp["seq_submit_count"].get(system_name, 4)))
    cap = cp["culture_and_preserve"]
    lines.append("**培养/保存:** " + cap["ep_label"]["rule"] + "; 小提质粒洗脱 " +
                 cap["mini_prep"]["elution"] + "; 大提(重力柱)洗脱 " +
                 cap["maxi_prep_gravity"]["elution_range"] + "; " + cap["passage_note"]["rule"])
    ia = cp["info_archive"]
    lines.append("**归档:** ①" + ia["step1_primer"] + " ②" + ia["step2_after_build"] +
                 " ③" + ia["step3_seq_done"] + " ④" + ia["step4_glycerol_stock"])
    gs = cp["glycerol_stock_storage"]
    lines.append("**甘油菌库存:** " + gs["work_stock"] + " / " + gs["long_term_stock"] +
                 "; " + gs["synthesis_note"])
    return "\n".join(lines)


# ==================== 报告渲染 ====================
def render_report(system_name, seq_name, items, primers, sop_text, codon_score=None):
    order = {"ERROR": 0, "WARNING": 1, "INFO": 2}
    items.sort(key=lambda x: order[x["level"]])
    errs = sum(1 for i in items if i["level"] == "ERROR")
    warns = sum(1 for i in items if i["level"] == "WARNING")
    infos = sum(1 for i in items if i["level"] == "INFO")
    out = []
    out.append("# 分子构建校验报告")
    out.append("")
    out.append("- 系统: **{}** | 序列: **{}**".format(system_name, seq_name or "未命名"))
    out.append("- 汇总: 🔴 ERROR **{}** 项 | 🟡 WARNING **{}** 项 | ℹ️ INFO **{}** 项".format(errs, warns, infos))
    out.append("")
    out.append("## 【严重错误 ERROR】(阻断, 必须修改)")
    for i in items:
        if i["level"] == "ERROR":
            out.append("- **{}** | 位置: {}  \n  原因: {}  \n  整改: {}".format(
                i["rule"], i["pos"], i["detail"], i["suggestion"]))
    out.append("")
    out.append("## 【风险警告 WARNING】(风险复核)")
    for i in items:
        if i["level"] == "WARNING":
            out.append("- **{}** | 位置: {}  \n  原因: {}  \n  建议: {}".format(
                i["rule"], i["pos"], i["detail"], i["suggestion"]))
    out.append("")
    out.append("## 【信息提示 INFO】(备忘)")
    for i in items:
        if i["level"] == "INFO":
            out.append("- **{}** | 位置: {} — {}".format(i["rule"], i["pos"], i["detail"]))
    out.append("")
    out.append("## 同源重组引物表")
    out.append("- 同源臂长度: {} bp (SOP: 15–25 bp)".format(primers["homology_arm_len"]))
    out.append("- 上游引物(5'→3'): `{}`".format(primers["fwd"]))
    out.append("- 下游引物(5'→3'): `{}`".format(primers["rev"]))
    out.append("- 扩增结合段(供核对): F=`{}` R=`{}`".format(primers["fwd_seq"], primers["rev_seq"]))
    out.append("- 引物浓度: {}; {}".format(primers["conc"], primers["plate"]))
    out.append("- ⚠️ " + primers["note"])
    out.append("")
    out.append("## 实验 SOP 备忘")
    out.append(sop_text)
    # 密码子偏好打分 (补充输出)
    if codon_score:
        out.append("")
        out.append("## 密码子偏好打分（按宿主偏好密码子）")
        out.append("- 宿主: **{}** | 得分: **{}** ({})".format(
            codon_score["host"],
            "{:.0%}".format(codon_score["score"]) if codon_score["score"] is not None else "n/a",
            codon_score["level"]))
        out.append("- 命中偏好密码子: {} / {} 个可比密码子 | CDS GC: {}%".format(
            codon_score["matched"], codon_score["comparable"], codon_score["gc_percent"]))
        if codon_score["nonpreferred_count"]:
            out.append("- 非偏好密码子 {} 个，示例: {}".format(
                codon_score["nonpreferred_count"],
                ", ".join("{}→{}({}{})".format(d["codon"], d["preferred"], d["aa"], d["pos"])
                          for d in codon_score["nonpreferred_detail"][:5])))
        out.append("- ❓ 是否需要在不改变氨基酸序列的前提下，优化为 {} 的最优密码子偏好序列？（确认后执行优化）".format(
            codon_score["host"]))
    return "\n".join(out)


def summarize(items):
    """汇总各级别条目数量"""
    return {
        "error": sum(1 for i in items if i["level"] == "ERROR"),
        "warning": sum(1 for i in items if i["level"] == "WARNING"),
        "info": sum(1 for i in items if i["level"] == "INFO"),
    }


def run_validation(system_name, seq_text, seq_name="", vector="", promoter="", host=None):
    """核心校验入口: 返回结构化 dict。
    供 CLI (main) / HTTP API (api_server.py) / 其他 Agent 程序直接调用。
    """
    if system_name not in ("mrna", "mammalian", "prokaryotic"):
        raise ValueError("系统必须为 mrna / mammalian / prokaryotic")
    rules = load_rules()
    seg = parse_input(seq_text)
    v = Validator(rules)
    items = v.validate(system_name, seg, seq_name=seq_name, vector=vector, promoter=promoter)
    primers = design_primers(seg)
    sop = sop_memo(system_name, rules)
    # 补充输出: 按系统对应宿主的密码子偏好打分 (不自动优化, 交用户确认)
    codon_score = score_codon_preference(system_name, seq_text, host=host, rules=rules)
    return {
        "system": system_name,
        "sequence": seq_name,
        "summary": summarize(items),
        "items": items,
        "primers": primers,
        "sop": sop,
        "codon_preference": codon_score,
        "optimization_prompt": "是否已按 {} 偏好密码子对 CDS 打分（{:.0%}），是否需要在不改变氨基酸序列的前提下优化为最优密码子偏好序列？".format(
            codon_score["host"], codon_score["score"] or 0) if codon_score["score"] is not None else "CDS 无可比密码子，跳过密码子优化。",
    }


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    system_name = sys.argv[1].lower()
    seq_file = Path(sys.argv[2])
    vector = ""
    promoter = ""
    host = None
    do_optimize = "--optimize" in sys.argv
    use_json = "--json" in sys.argv
    # 解析可选参数
    args = sys.argv[3:]
    for i, a in enumerate(args):
        if a == "--vector" and i + 1 < len(args):
            vector = args[i + 1]
        if a == "--promoter" and i + 1 < len(args):
            promoter = args[i + 1]
        if a == "--host" and i + 1 < len(args):
            host = args[i + 1]

    if system_name not in ("mrna", "mammalian", "prokaryotic"):
        sys.exit("[FATAL] 系统必须为 mrna / mammalian / prokaryotic")

    text = seq_file.read_text(encoding="utf-8")
    seq_name = ""
    for ln in text.splitlines():
        if ln.startswith(">"):
            seq_name = ln[1:].strip()
            break

    if do_optimize:
        # 优化模式: 输出最优密码子偏好序列 (仅在用户确认后调用)
        rules = load_rules()
        result = optimize_sequence(system_name, text, host=host, rules=rules)
        result["sequence"] = seq_name
        if use_json:
            import json
            print(json.dumps(result, ensure_ascii=False, indent=2))
        else:
            print("# 密码子偏好优化结果 ({})\n".format(result["host"]))
            print("- 得分: {:.0%} → {:.0%} | GC: {:.1f}% → {:.1f}%".format(
                result["score_before"] or 0, result["score_after"] or 0,
                result["gc_before"], result["gc_after"]))
            print("- 替换密码子: {} 个 | 氨基酸序列不变: {}".format(
                result["changes_count"], result["protein_unchanged"]))
            print("\n## 优化后 CDS\n```\n{}\n```".format(result["optimized_cds"]))
            print("\n## 优化后完整序列\n```\n{}\n```".format(result["optimized_sequence"]))
            print("\n## 变更明细\n" + "\n".join(
                "- {pos}: {from} → {to} ({aa})".format(**c) for c in result["changes"][:30]))
        sys.exit(0)

    result = run_validation(system_name, text, seq_name=seq_name, vector=vector, promoter=promoter, host=host)

    if use_json:
        import json
        out = json.dumps(result, ensure_ascii=False, indent=2)
        out_file = BASE_DIR / "report_{}.json".format(system_name)
        out_file.write_text(out, encoding="utf-8")
        print(out)
        print("\n[JSON已保存] " + str(out_file), file=sys.stderr)
    else:
        report = render_report(system_name, seq_name, result["items"], result["primers"], result["sop"],
                               result.get("codon_preference"))
        out_file = BASE_DIR / "report_{}.md".format(system_name)
        out_file.write_text(report, encoding="utf-8")
        print(report)
        print("\n[报告已保存] " + str(out_file), file=sys.stderr)


if __name__ == "__main__":
    main()
