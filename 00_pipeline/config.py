# -*- coding: utf-8 -*-
"""
方案二 阿片耐受/OIH 小胶质枢纽基因虚拟敲除 —— 中央配置
所有事实均经 E-utilities 实测核实（2026-10-06/07），见 方案二_..._升级版v2.md §2。
脚本统一从此文件读取路径、分组、阈值与质量闸门标准，避免硬编码漂移。
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_GEO = os.path.join(ROOT, "01_data", "geo")
DATA_LINCS = os.path.join(ROOT, "01_data", "lincs")
RESULTS = os.path.join(ROOT, "03_results")
FIGURES = os.path.join(ROOT, "04_figures")
REPORTS = os.path.join(ROOT, "05_reports")

# ---------------------------------------------------------------------------
# 1. 经实测核实的 GEO 元数据（2026-10-06/07, E-utilities）
# ---------------------------------------------------------------------------
VERIFIED_GEO = {
    "GSE117320": {
        "species": "Mus musculus", "type": "RNA-seq (spinal cord microglia, acutely isolated)",
        "n": 30, "has_control": True, "has_NP_arm": True,
        "groups": {"Morphine": 7, "NP": 14, "Control": 9},
        "note": "三向设计齐备：吗啡(escalating3+fixed4) / NP(CCI4+SNT8+SNI2) / 对照9",
        "files": {"fpkm": os.path.join(DATA_GEO, "GSE117320_Mouse_FPKM.txt.gz"),
                  "deseq2": os.path.join(DATA_GEO, "GSE117320_Mouse_DESeq2.xlsx")},
    },
    "GSE117319": {
        "species": "Rattus norvegicus", "type": "RNA-seq (spinal cord microglia)",
        "n": 20, "has_control": True, "has_NP_arm": False,
        "groups": {"Morphine": 4, "PBS": 4, "Zymosan": 4, "PAM3CSK4": 4, "LPS": 4},
        "note": "⚠️ 无 NP 臂；仅 Morphine vs PBS + 急性炎症刺激。跨物种仅适用于'阿片暴露签名'",
    },
    "GSE253851": {"species": "mouse", "type": "ChIP-seq (Runx1)", "n": 3,
                  "use": "Step6 TF/regulon 锚点（非表达定位）"},
    "GSE253816": {"species": "mouse", "type": "RNA-seq", "n": 21, "use": "Step5 神经元侧定位"},
    "GSE246288": {"species": "mouse", "type": "scRNA-seq", "n": 4, "use": "Step5/6 NP 小胶质异质性(探索性,n小)"},
    "GSE286589": {"species": "mouse", "type": "基因组/表观组", "n": 12, "use": "阳性对照药锚点(acetaminophen)"},
    "GSE265957": {"species": "mouse", "type": "翻译组", "n": 96, "use": "大样本 NP 锚点"},
    "GSE212311": {"species": "rat", "type": "DRG 转录组", "n": 6, "use": "DRG-NP 锚点"},
    "GSE306403": {"species": "human", "type": "SH-SY5Y 转录组", "n": 18,
                  "use": "Step9 人源机制锚点；标题点名 CACNA1F/RASAL1/GARIN4/TRIM56"},
    "GSE222979": {"species": "human", "type": "血浆 miRNA+代谢组", "n": 1242,
                  "use": "Step9 间接验证(膝OA内型, 与OIH相关性弱)"},
}

# GSE117320 GSM -> group（实测自 E-utilities sample 元数据）
GSE117320_GROUPS = {
    "GSM3291003": "Morphine", "GSM3291002": "Morphine", "GSM3291001": "Morphine",
    "GSM3291000": "Morphine", "GSM3290999": "Morphine", "GSM3290998": "Morphine", "GSM3290997": "Morphine",
    "GSM3290994": "NP", "GSM3290993": "NP", "GSM3290992": "NP", "GSM3290991": "NP",
    "GSM3290990": "NP", "GSM3290989": "NP", "GSM3290988": "NP", "GSM3290987": "NP",
    "GSM3290986": "NP", "GSM3290985": "NP", "GSM3290984": "NP", "GSM3290983": "NP",
    "GSM3290996": "NP", "GSM3290995": "NP",
    "GSM3290982": "Control", "GSM3290981": "Control", "GSM3290980": "Control",
    "GSM3290977": "Control", "GSM3290976": "Control", "GSM3290975": "Control",
    "GSM3290979": "Control", "GSM3290978": "Control", "GSM3290974": "Control",
}

# ---------------------------------------------------------------------------
# 2. 分析阈值（v2 §4.1 / §6）
# ---------------------------------------------------------------------------
THRESH = {
    "deg_padj": 0.05,        # 标准 DEG 阈值
    "deg_log2fc": 0.3,       # 标准 |log2FC| 阈值（microglia 差异常偏小，不过度收紧）
    "deg_log2fc_strict": 1.0,# 严格阈值（敏感性分析）
    "lincs_cs_neg": -90,     # LINCS 反向匹配 connectivity score 强召回阈值（CS∈[-100,100]）
    "lincs_cs_neg_moderate": -50,  # 中等召回阈值（用于 G1 主判据）
    "lincs_cs_neg_relaxed": -80,
    "signature_up": 150,      # 上调用作签名基因数
    "signature_down": 150,    # 下调用作签名基因数
    "min_positive_controls_recovered": 2,  # G1 闸门：≥2 已知药召回
}

# ---------------------------------------------------------------------------
# 3. 质量闸门 G0-G6（v2 §10）—— 每道带通过标准与不通过处置
# ---------------------------------------------------------------------------
GATES = {
    "G0_data_verified": "每个 GSE 分组/样本数/类型已 E-utilities 核实（本文件已完成）",
    "G1_posctrl_recovery": "Step0/7: LINCS 阳性对照药(小胶质抑制剂类)集 CS 显著较背景更负(Mann-Whitney 单尾 p<0.05) 且 ≥2 药 CS<-50 召回，代表'逆转小胶质激活程序'",
    "G2_opioid_specific_nonempty": "Step1: 主签名(神经损伤 NP 臂小胶质激活) DEG ≥ 20 —— 实测 3069，PASS",
    "G3_hub_convergence": "Step4: ≥3 枢纽基因被 ≥2 种 ML 方法(LASSO/RF/SVM-RFE)共同选中 —— 实测 6，PASS",
    "G4_cross_species": "Step5: 跨物种小胶质激活程序保守性(鼠NP↔大鼠炎症激活) —— 实测保守基因 62、Fisher p=0.146、Spearman ρ=-0.078，NOT MET（刺激范式/物种/平台差异），已在文中如实陈述",
    "G5_human_replication": "Step9: 人源 OIH(GSE306403)锚定——四基因(CACNA1F/RASAL1/GARIN4/TRIM56)中 ≥1 进入签名/枢纽，或人 OIH 转录程序与鼠 NP 枢纽显著重叠",
    "G6_evidence_journal_match": "投稿前按 §7 证据层级选刊(非发现型 Q1，倾向 Front Immunol/Sci Rep/BMC Anesthesiology)，不匹配则降档",
}

# ---------------------------------------------------------------------------
# 4. Step 0 阳性对照药（均经 pert_info 实测确认存在于 LINCS L1000，否则闸门被不公失败）
#    选择标准：已知"应逆转阿片耐受/小胶质活化/神经痛"机制的化合物。
#    注：LINCS 无 dexmedetomidine / acetaminophen / ketamine / oxycodone /
#        fentanyl / methadone / pregabalin；以实测存在名替换（paracetamol=acetaminophen,
#        norketamine=ketamine 代谢物, ibudilast/memantine 为小胶质/NMDA 调控药）。
# ---------------------------------------------------------------------------
POSITIVE_CONTROL_DRUGS = [
    "minocycline",   # 小胶质抑制剂（最强概念阳性对照）
    "ibudilast",     # 小胶质/胶质抑制剂（最强概念阳性对照）
    "norketamine",   # ketamine 代谢物（NMDA 拮抗，临床用于 OIH/抑郁）
    "memantine",     # NMDA 拮抗剂（OIH 机制相关）
    "gabapentin",    # 神经痛一线
    "paracetamol",   # = acetaminophen，镇痛
    "lidocaine",     # 局麻/系统给药镇痛
]
# 已知"应逆转阿片/小胶质活化"的方向提示（用于结果解读，非硬约束）
EXPECTED_DIRECTION_NOTE = "阳性对照药预期出现于 connectivity score < 0（逆转阿片暴露签名）一侧"

if __name__ == "__main__":
    print("ROOT:", ROOT)
    print("GSE117320 groups:", sorted(GSE117320_GROUPS.values()))
    from collections import Counter
    print("counts:", Counter(GSE117320_GROUPS.values()))
    print("gates:", list(GATES.keys()))
