# Task_Convergence: 热容海平面 (TSLA) 收敛性与混合项分析 (Convergence & Cross-Term Analysis)

## 📌 项目背景与目的 (Purpose)

本文件夹下的脚本旨在回答以下三个核心科学问题，以验证我们在计算热容海平面 (Thermosteric Sea Level, TSLA) 时所采用的近似方法的准确性：

1.  **收敛性验证 (Convergence)**：在泰勒展开计算 TSLA 时，保留到 **第 3 阶 (Order 3)** 是否足够？残差 (Residual) 是否已收敛至可忽略的程度？
2.  **缺项补充 (Missing Physics)**：之前的分析主要集中在 **纯温度项 (Pure T)**。本分析补充了被忽略的 **纯盐度项 (Pure S)** 和 **温盐耦合混合项 (Cross Terms)**，以评估其量级与影响。
3.  **近似准确性 (Accuracy)**：通过对比 **"精确解 (Exact)"**（直接由非线性状态方程计算）与 **"近似解 (Approx)"**（由泰勒展开各项求和），定量评估误差。

---

## 📁 文件夹结构 (Directory Structure)

```
Task_Convergence/
├── 01_Data_Preprocessing/      # 数据预处理：计算平均态
├── 02_Exact_Calculation/       # 精确解计算：非线性状态方程
├── 03_Terms_Calculation/       # 泰勒展开项计算：T项、S项、混合项
├── 04_CrossDetail_Calculation/ # 详细混合项计算：各阶混合项细节
├── 05_Convergence_Analysis/    # 收敛性分析：精确解 vs 近似解
├── 06_OrderBreakdown_Analysis/ # 阶数分解分析：各阶贡献
├── 07_CrossTerms_Analysis/     # 混合项专项分析：量级与分布
├── 08_MultiModel_Analysis/     # 多模式验证：EN4/IAP/Ishii
├── 09_Utils/                   # 工具函数
├── 10_Debug_Scripts/           # 调试脚本
└── README.md                   # 本文件
```

---

## 🚀 使用指南与执行顺序 (Workflow)

请严格按照以下步骤顺序运行 MATLAB 脚本，以确保数据依赖关系正确。

### **第一阶段：数据预处理 (Data Preprocessing)**

运行 `01_Data_Preprocessing/` 中的脚本：
- `Calc_Mean_IAP.m` → 计算 IAP 平均态
- `Calc_Mean_Ishii.m` → 计算 Ishii 平均态

### **第二阶段：数据生成 (Data Generation)**
*耗时较长，需计算全非线性解和高阶导数项。*

运行 `02_Exact_Calculation/` 中的脚本：
- `Calc_Exact_EN4_Avg.m` / `Calc_Exact_EN4_Std.m`
- `Calc_Exact_IAP_Avg.m` / `Calc_Exact_IAP_Std.m`
- `Calc_Exact_Ishii_Avg.m` / `Calc_Exact_Ishii_Std.m`

运行 `03_Terms_Calculation/` 中的脚本：
- `Calc_T_Terms_*.m` → 纯温度项
- `Calc_S_Terms_*.m` → 纯盐度项
- `Calc_Cross_and_S_Terms_Avg.m` → 混合项 (平均态)
- `Calc_Std_Terms_HighOrder.m` → 高阶项 (标准态)

运行 `04_CrossDetail_Calculation/` 中的脚本：
- `Calc_CrossDetail_*.m` → 详细混合项

### **第三阶段：分析与可视化 (Analysis & Visualization)**

| 分析类型 | 脚本位置 | 功能 |
|----------|----------|------|
| 收敛性分析 | `05_Convergence_Analysis/` | 精确解 vs 近似解对比 |
| 阶数分解 | `06_OrderBreakdown_Analysis/` | 各阶项贡献分析 |
| 混合项分析 | `07_CrossTerms_Analysis/` | 混合项量级与分布 |
| 多模式验证 | `08_MultiModel_Analysis/` | EN4/IAP/Ishii 对比 |

---

## 📊 各文件夹说明

| 文件夹 | 功能 | 文件数 |
|--------|------|--------|
| `01_Data_Preprocessing` | 计算海洋数据平均态 | 2 |
| `02_Exact_Calculation` | 计算精确解 (无近似) | 6 |
| `03_Terms_Calculation` | 计算泰勒展开各项 | 10 |
| `04_CrossDetail_Calculation` | 计算详细混合项 | 6 |
| `05_Convergence_Analysis` | 收敛性验证分析 | 4 |
| `06_OrderBreakdown_Analysis` | 阶数贡献分解 | 2 |
| `07_CrossTerms_Analysis` | 混合项专项分析 | 3 |
| `08_MultiModel_Analysis` | 多模式验证分析 | 2 |
| `09_Utils` | 工具函数 | 1 |
| `10_Debug_Scripts` | 调试脚本 | 7 |

---

## 📁 数据依赖说明 (Data Dependencies)

*   **输入数据源**:
    *   `D:\work\EN4_analyses_c13_last20years\` (原始 NC 数据)
    *   `D:\work\IAP_05_24\` (IAP 数据)
    *   `D:\work\Ishii_05_24\` (Ishii 数据)
*   **中间数据**:
    *   所有生成的 `.mat` 文件保存在 `D:\work\MAT_Data\` 或 `D:\work\*_TSLA_Terms\`
*   **结果图片**:
    *   保存在 `D:\work\Figures\Convergence\`, `StdRef\`, `Average\` 等子文件夹中。

---

## 🔧 常见问题 (FAQ)

*   **Q: 为什么计算 Exact 解很慢？**
    *   A: `Calc_Exact_*.m` 需要逐月读取 20 年的全球 3D 数据，并对每个格点进行高精度的 TEOS-10 密度计算 (涵盖 42 层深度)，计算量巨大，通常需要 30-60 分钟。

*   **Q: 标准态 (Std) 和 平均态 (Avg) 有什么区别？**
    *   A:
        *   **平均态**: 在每个格点上，以当地的时间平均 T/S 为展开中心。优点是展开项收敛快 (3阶只需计算很少项)，物理意义直观（距平）。
        *   **标准态**: 全球统一使用 $S=35, \theta=0$ 为展开中心。优点是基准统一，但由于某些海区（如地中海、极地）距离展开中心较远，需要更高阶（如 8 阶）才能收敛。

*   **Q: 遇到数据读取错误怎么办？**
    *   A: 运行 `10_Debug_Scripts/` 中的检查脚本，确认变量名和维度是否正确。
