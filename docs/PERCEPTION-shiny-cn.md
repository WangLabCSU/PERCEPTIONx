# PERCEPTION-shiny 使用教程

**PERCEPTION-shiny** 是 **PERCEPTIONx** R 包的交互式网页应用（Shiny）。它把完整的分析流水线——数据加载、模型训练、药物敏感性预测、结果可视化——封装成点选式界面，无需写代码即可完成从"患者单细胞表达矩阵"到"克隆级杀伤分数 + 患者级响应分层"的全流程分析。

> 方法学基础：PERCEPTION（PERsonalized single-Cell Expression-based Planning for Treatments In ONcology），在 DepMap 细胞系数据上训练弹性网络模型，预测患者对药物的响应与耐药。

---

## 目录

- [0. 环境要求与启动](#0-环境要求与启动)
- [1. 界面总览](#1-界面总览)
- [2. Data 页：加载数据](#2-data-页加载数据)
- [3. Train 页：训练模型（可选）](#3-train-页训练模型可选)
- [4. Predict 页：预测杀伤值](#4-predict-页预测杀伤值)
- [5. Visualize 页：可视化](#5-visualize-页可视化)
- [6. Help 页](#6-help-页)
- [7. 常见问题（FAQ）](#7-常见问题faq)
- [8. 引用与联系](#8-引用与联系)

---

## 0. 环境要求与启动

### 0.1 环境要求

- **R** ≥ 4.1.0
- 主要依赖包：`devtools`、`shiny`、`bslib`、`Seurat`、`ggplot2`、`ggiraph`、`glmnet`、`caret`（缺失时应用会在相应功能处提示安装）

### 0.2 启动方式

在包源码根目录打开 R/RStudio，运行：

```r
devtools::load_all()          # 从源码加载 PERCEPTIONx
run_perception_app()          # 启动应用（自动打开浏览器）
```

或直接运行源码目录中的应用：

```r
shiny::runApp("inst/shiny/app")
```

### 0.3 整体流程

```
DepMap 参考数据 ──► 模型训练 ──► 克隆/患者预测 ──► 可视化与验证
      ▲               ▲             ▲                 ▲
   Data 页         Train 页     Predict 页        Visualize 页
（或直接加载     （或直接用    （选模型后
  预训练模型）    44 药预训练   一键预测）
                 模型，跳过）
```

> 提示：想最快看到效果，按 **Load Demo →  Predict → Visualize** 的顺序点击即可，全程不需要训练。

---

## 1. 界面总览

![BQACAgUAAyEGAASHRsPbAAEYrqVqdGKkjlcuARJe5YzI5KXXCeDMvgACfjQAAtsboFcyiuFv_R2_bj0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYrqVqdGKkjlcuARJe5YzI5KXXCeDMvgACfjQAAtsboFcyiuFv_R2_bj0E.png)

顶部导航栏包含 6 个标签页：**Home、Data、Train、Predict、Visualize、Help**。数据从左向右流动：先加载数据（Data），再训练/加载模型（Train），然后预测（Predict），最后可视化与验证（Visualize）。

首页包含：介绍、四步引导（Load Data → Train Model → Predict → Visualize，点击可跳转到对应页面）、数据加载状态一览、功能特性与引用信息。页面上的 **Quick Start** 与 **Load Demo** 按钮可一键加载演示数据。

---

## 2. Data 模块：加载数据

Data 页是分析的起点，负责加载四种数据：**演示数据 / DepMap 参考数据 / 用户表达矩阵 / 临床响应**。每一项加载成功后，左侧状态徽章会变绿。

> **Load Demo**：
>
> 点击 **Load Demo** 一键生成合成演示数据：**49 个基因 × 400 个细胞 × 20 例患者**，并自动完成 Seurat 聚类、秩归一化，随后训练演示模型。用于快速测试整个流程，感受各功能页的交互方式。

### 2.1 加载 DepMap 参考数据（训练必需）

两种方式：

1. **Download & Load**：点击按钮自动从官方镜像下载 DepMap 参考集（约 567 MB，1.5 万+ 基因 × 1,000+ 细胞系），下载完成后自动加载。这是标准训练/参考输入，也是内存和磁盘占用最大的步骤。
2. **本地上传 .RDS**：如果本地已有 DepMap 文件（`DepMap.RDS` 格式），直接 Browse 选中即自动加载。

### 2.2 加载模型

两种方式：

1. **Download & Load**：一键下载 44 种 FDA 批准药物的预训练模型（如 `abemaciclib`、`erlotinib`、`osimertinib`），多选时每项右侧有 × 可单独删除
2. **本地上传 .RDS**：选择训练好的模型文件（`train_models()` 或 Train 页导出的 RDS），选中即自动加载

### 2.3 上传表达矩阵（患者 scRNA-seq）

- **支持的格式**：CSV / TSV / TXT / Excel（.xlsx/.xls）/ RDS
- **支持的 R 对象**（RDS）：数值矩阵、data.frame、Seurat 对象
- **矩阵方向**：基因 × 细胞（行为基因、列为细胞）；若第一列是基因名字符列（Excel/csv 导出常见），会自动转为行名
- **关于归一化**：表达数据需要是**秩归一化**（rank-normalized）后使用。若上传的是原始计数，勾选后由应用在聚类流程中自动归一化。

### 2.4 上传细胞-患者映射（Mapping）

- **支持的格式**：CSV / TSV / TXT / Excel / RDS
- **必须包含两列**：`cell_id`（细胞名，与表达矩阵列名一致）和 `patient_id`（细胞所属患者）。列名大小写均可（`Patient`、`PATIENT` 都能识别）
- **RDS 命名列表**：若 RDS 是"患者 → 细胞名向量"的命名列表（如论文 demo 的 `PRJNA591860_sample_cell_names.RDS`），会自动转换为长格式；无细胞的空样本自动丢弃

### 2.5 上传临床响应（Response，可选但推荐）

- **支持的格式**：CSV / TSV / TXT / Excel / RDS
- **必须包含两列**：`patient`（与 Mapping 的 `patient_id` 完全一致）和 `response`（患者对药物的真实临床结局）
- **响应值自动归一化**：`responder` / `responsive` / `r` / `sensitive` → **Responder**；`non-responder` / `nr` / `resistant` → **Non-responder**
- 上传时会自动检查 patient 是否与 Mapping 对齐；若完全不重叠会发出警告
- **为什么需要**：只有知道真实结局，才能验证预测准不准（画 ROC、响应者/非响应者箱线图）。只想看预测分数可以不传

### 2.6 聚类分群（Seurat）

*Don't Forget to Run Seurat!*

可选择分群方法：UMAP / tSNE。

| 参数 | 说明 | 建议 |
|---|---|---|
| Seurat Resolution | 聚类分辨率，越大分群越细 | 0.5–1.0 |
| Seurat Dims | PCA 降维使用的维度数 | 10–30 |
| Seurat NFeatures | 高变基因数（用于聚类） | 1000–3000，默认 2000 |

> Notes：聚类用于定义clone（细胞亚群），后续预测会以克隆为单位进行。

![BQACAgUAAyEGAASHRsPbAAEYrrpqdGNtTytttAcxi0NLh0PzdHU_SQAClTQAAtsboFfopMIK8x6Grz0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYrrpqdGNtTytttAcxi0NLh0PzdHU_SQAClTQAAtsboFfopMIK8x6Grz0E.png)

---

## 3. Train 模块：训练模型（可选）

**前置条件**：需先在 Data 页加载 DepMap 参考数据（若缺失，页面顶部会提示）。

> 如果只是想用现有模型预测，**完全可以不训练**——直接用 Data 页加载的 44 个预训练模型即可。Train 页适合：换药物、换癌症类型、自定义基因集、想查看模型在三种数据验证集上的表现。

### 3.1 参数说明

| 参数 | 说明 |
|---|---|
| **Drug Name** | 药物多选（支持组合用药），每项带 × 可删除 |
| **Cancer Type (include)** | 训练包含的癌种，默认 PanCan（泛癌） |
| **Cancer Type (exclude)** | 排除的癌种（防止自我验证） |
| **Gene Symbols** | 基因列表，留空 = 使用 DepMap 全部基因（推荐）；可粘贴文本或上传 .txt/.csv |
| **Top k Features** | 特征排名取前 k 个基因进入模型 |
| **Algorithm** | 弹性网络 `glmnet`（推荐）或随机森林 `rf` |
| **CPU Cores** | 并行核数 |

### 3.2 结果解读

- **Model Summary**：每个药物的模型类型与超参数（glmnet 的 alpha/lambda，或 rf 的 ntree/RMSE）
- **Performance Plot**：模型性能曲线（阈值-相关关系）
- **Performance Metrics**：模型在 **Bulk / Pseudo-bulk / Single-cell** 三个层面上的预测-真实相关系数与 p 值。相关系数越高、p 值越小 = 模型越好
- **Download Model (.RDS)**：导出训练好的模型，可日后在 Data 页直接上传复用

---

## 4. Predict 模块：预测杀伤值

选择已加载的模型（Data 页加载的预训练模型或 Train 页训练出的模型），点击预测：

1. **Clone-level（克隆级）**：每个克隆 × 每个药物的活力分数。分数语义：模型输出的是 **viability（存活度）**，**值越高 = 越耐药、越低 = 越敏感**；页面展示的 `viability_scaled` 已归一化到 [0, 1]——**越接近 1 表示越耐药（活力越高）**
2. **Patient-level（患者级）**：按克隆占比把克隆分数聚合到患者（默认 `weighted_max`），得到每例患者的药物敏感性分层

输出包括交互式热图（克隆 × 药物，plotly）与可下载的预测表格。

![BQACAgUAAyEGAASHRsPbAAEYrr1qdGOEsPafeye9qfPgBPkIPFAL_gACmDQAAtsboFdiPjfKVPqPXD0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYrr1qdGOEsPafeye9qfPgBPkIPFAL_gACmDQAAtsboFdiPjfKVPqPXD0E.png)

---

## 5. Visualize 模块：可视化

此模块用于将预测结果进一步可视化，故而必须先进行预测。

所有图表均为**交互式 SVG**（基于 ggiraph）：鼠标悬停任意点/条即可看到详细信息（克隆 ID、杀伤分数、占比、FPR/TPR 等）。图中工具栏支持放大、平移、下载。

### 5.1 Clone Distribution（stackplot）

![BQACAgUAAyEGAASHRsPbAAEYrsRqdGOsyEvZKa-hX56p8XZe9MdfBAACnzQAAtsboFciL8mgr8EQBT0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYrsRqdGOsyEvZKa-hX56p8XZe9MdfBAACnzQAAtsboFciL8mgr8EQBT0E.png)

展示每例患者内部各克隆的占比构成，一个色带 = 一个克隆（≤15 个克隆时使用内置协调色板）。

> 注意：克隆身份取决于数据来源——全局聚类路径下，同一克隆（同色）确实跨患者共享；克隆级输入（如按患者编号的 c1/c2/c3）时，**跨患者颜色相同不代表同一克隆来源**，它只是共享的类别标签。

### 5.2 Clone Viability（lolliplot）

![BQACAgUAAyEGAASHRsPbAAEYrtRqdGQW8zkwfYWQ_VN994k_iZJH8gACsDQAAtsboFe7I7H3sbW6sz0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYrtRqdGQW8zkwfYWQ_VN994k_iZJH8gACsDQAAtsboFe7I7H3sbW6sz0E.png)

- **展示规则**：全样本、全克隆展示，一例患者一个分面，一根棒 = 一个克隆
- **颜色**：蓝-白-红发散——蓝 = 预测敏感（低活力），红 = 预测耐药（高活力）
- **点大小**：克隆占比（占比大的克隆点大）
- **排序**：患者内按占比降序；有响应数据时 Responder 排在前面
- **y 轴**：Predicted Viability（z-score），0 线为中位参考

### 5.3 ROC Curve

![BQACAgUAAyEGAASHRsPbAAEYrtlqdGQvuGR0iVJf_0IC3FrHl2wx2gACtTQAAtsboFdEeKvbpn7yAT0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYrtlqdGQvuGR0iVJf_0IC3FrHl2wx2gACtTQAAtsboFdEeKvbpn7yAT0E.png)

使用真实临床响应（Data 页上传的 Response）与患者级预测分数绘制 ROC 曲线，AUC 越接近 1 说明预测分层能力越强。若为多药模型预测，可选取要查看的特定模型。

### 5.4 Response Boxplot（R vs NR）

![BQACAgUAAyEGAASHRsPbAAEYruFqdGRmVjR7z6AiKJZrmoSvML7-BgACvTQAAtsboFcaROROGDpEdj0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYruFqdGRmVjR7z6AiKJZrmoSvML7-BgACvTQAAtsboFcaROROGDpEdj0E.png)

分组展示 Responder 与 Non-responder 的预测分数分布，附带显著性检验结果。

### 5.5 Model Performance

![BQACAgUAAyEGAASHRsPbAAEYruZqdGSFqhl8FAYozOIa3fcN_FgDUgACwjQAAtsboFd6-2EPyhNE6D0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYruZqdGSFqhl8FAYozOIa3fcN_FgDUgACwjQAAtsboFd6-2EPyhNE6D0E.png)

展示训练模型的验证性能曲线（需 Train 页训练的模型；预训练模型不适用，页面会提示）。

### 5.6 空间图（UMAP / t-SNE）

![BQACAgUAAyEGAASHRsPbAAEYru1qdGSiWp9UYiS0DlsprEKVOud1_AACyTQAAtsboFfwQwHtvh_32z0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYru1qdGSiWp9UYiS0DlsprEKVOud1_AACyTQAAtsboFfwQwHtvh_32z0E.png)

![BQACAgUAAyEGAASHRsPbAAEYrvNqdGS9_DCHZMgsGZFMcIhxBlsf_wACzzQAAtsboFerQBd3_vv2kz0E.png](https://img.remit.ee/api/file/BQACAgUAAyEGAASHRsPbAAEYrvNqdGS9_DCHZMgsGZFMcIhxBlsf_wACzzQAAtsboFerQBd3_vv2kz0E.png)

选择降维方法和着色变量：

- **Gene Expression**：每个细胞的单基因真实表达量（z-score），连续色标——看基因在哪群细胞高表达
- **Drug Viability**：每个细胞回填所属克隆的预测活力分数，呈"拼接色块"——看哪些克隆被预测耐药（亮）或敏感（暗）
- **Clone / Cluster**：按克隆（聚类）着色——看分群结构

> 对照读法：Gene Expression 高表达的细胞群，如果在 Drug Viability 图中也是暖色（高活力），说明该基因高表达与预测耐药正相关；反之（暗色/低活力）则与敏感相关。注意两图色标刻度不同，只看空间分布模式。
>
> 例如，根据上方两个示例图能得出初步的有效信息：
>
> 1. 耐药标志物线索：最右下角的克隆中，SLC2A1 的高表达区与厄洛替尼预测高活力区在空间上存在重叠，提示 SLC2A1 高表达可能是该耐药克隆的标志物。
> 2. 敏感主群体：位于正上方的几个主要细胞群，SLC2A1 表达接近于零，却恰好是预测药物活力最低的区域。
>
> 因此，在PRJNA591860数据集中，SLC2A1 的局部异常高表达显著指向厄洛替尼耐药；而绝大多数不表达 SLC2A1 的细胞，对药物高度敏感。这为“SLC2A1 可能是耐药靶点”提供了初步的单细胞层面证据。

### 5.7 高清下载

| 格式 | 分辨率 |
|---|---|
| PNG | 600 dpi（Cairo 抗锯齿），默认约 6000 × 4200 像素 |
| PDF / SVG | 矢量格式，无限缩放，推荐投稿使用 |

---

## 6. Help 页

内置帮助文档，包含各页功能说明、FAQ 与使用建议。遇到问题时请优先查看。

---

## 7. 常见问题（FAQ）

**Q1：上传文件提示 "Maximum upload size exceeded"？**
应用已把上传上限设为 1 GB。超过该大小的数据请先本地处理（如按基因子集）再上传；DepMap 优先用 "Download & Load" 按钮（服务器直下，不经浏览器）。

**Q2：Model Performance 图提示缺少 performance 字段？**
该图需要训练流程产生的验证指标，预训练模型（`load_model`）不携带。请先用 Train 页训练，或在 Data 页加载演示模型。

**Q3：response 标签怎么构造？**
构造两列 `patient` / `response`。`patient` 必须与 Mapping 的 `patient_id` 完全一致；`response` 取值会自动归一化（responder → Responder 等）。没有现成结局数据时，可用治疗时间点代替（基线 → Responder，耐药进展 → Non-responder）。

**Q4：预测分数和我想的不一样？**
预测分数是模型基于 DepMap 训练得到的相对排名（已归一化到 0–1），"1" 表示在该批次克隆中最耐药（活力最高），不代表绝对杀伤比例。多克隆异质性、耐药通路激活都会拉高分数（活力升高）。

**Q5：关闭应用后数据会留下吗？**
不会。下载的 DepMap 数据、预训练模型均写入 R 会话的临时目录（tempdir），应用关闭自动清理，不会在磁盘残留。

---

## 8. 引用与联系

**若使用本应用/包，请引用原始方法学论文：**

> Sinha, S., Vegesna, R., Mukherjee, S. *et al.* PERCEPTION predicts patient response and resistance to treatment using single-cell transcriptomics of their tumors. *Nature Cancer* 5, 938–952 (2024). DOI: [10.1038/s43018-024-00756-7](https://doi.org/10.1038/s43018-024-00756-7)

**仓库**：[github.com/WangLabCSU/PERCEPTIONx](https://github.com/WangLabCSU/PERCEPTIONx)

**问题反馈**：jiading682@qq.com

---

*PERCEPTION-shiny © PERCEPTIONx authors. 本文档为使用教程，具体界面以实际版本为准。*
