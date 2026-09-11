# Triptolide与16-hydroxytriptolide配对比较对接结果

## 研究目的和证据边界

本分析用于补充“为什么计算筛选命中16-hydroxytriptolide，而后续结构模拟和实验使用triptolide（TPL）”之间的结构桥接证据。它不把TPL实验解释为16-hydroxytriptolide的实验验证，也不包含16-hydroxytriptolide的分子动力学或湿实验。

项目中的原始说明文件记录：16-hydroxytriptolide是L1000CDS2计算优先候选；由于其下游跟进存在限制，而结构类似物TPL更适合后续工作，因此原结构模拟和THP-1实验采用TPL。本次新增的是此前缺失的平行比较对接证据。

## 完成的计算

- 16-hydroxytriptolide来源：PubChem CID 126556，保持PubChem定义的立体化学。
- TPL来源：PubChem CID 107985。
- 两种化合物使用相同的S100A8受体、网格、评分函数和局部复核参数。
- 受体：5HLO的A/C和B/D两个S100A8同源二聚体晶体学拷贝。
- 全表面盲筛：24个重叠28 Å网格，exhaustiveness 32。
- 局部复核：24 Å界面盒，Vina和Vinardo两种评分函数，每个受体/评分组合10个相同固定种子，exhaustiveness 64。
- 完成24次16-hydroxytriptolide盲筛和40次局部对接，无失败；未进行16-hydroxytriptolide的MD或实验。

## 结构相似性

| 指标 | 结果 |
|---|---:|
| Morgan radius 2、2048位、含手性的Tanimoto | 0.831 |
| MACCS Tanimoto | 0.885 |
| 松弛手性限制后的最大公共子结构 | 26个重原子 |
| 公共子结构覆盖TPL | 26/26，100% |
| 公共子结构覆盖16-hydroxytriptolide | 26/27，96.3% |
| TPL分子式、分子量 | C20H24O6，360.406 Da |
| 16-hydroxytriptolide分子式、分子量 | C20H24O7，376.405 Da |
| TPSA | TPL 84.12 Å²；16-hydroxytriptolide 104.35 Å² |

两种化合物共享完整的26重原子TPL骨架，16-hydroxytriptolide多一个氧原子。但增加的羟基同时改变极性、氢键能力和疏水性，因此“结构类似”不能直接推导为功能等效。

## 全表面盲筛

- B/D二聚体中，两种化合物的全局最佳网格均为BD_tile11：TPL −7.757 kcal/mol，16-hydroxytriptolide −7.605 kcal/mol。
- A/C二聚体中，TPL的前两名依次为AC_tile02和AC_tile01；16-hydroxytriptolide的前两名顺序相反，为AC_tile01和AC_tile02。这两个28 Å网格相互重叠并覆盖同一界面区域。
- 因此，两种化合物均由独立全表面搜索定位至同一类S100A8二聚体界面区域，而不是仅通过强制局部盒得到相似结果。

## 局部多种子复核

| 化合物 | 受体 | 评分 | n | 平均分数 ± s.d. (kcal/mol) | 范围 | 姿势RMSD中位数/最大值 (Å) |
|---|---|---|---:|---:|---:|---:|
| 16-hydroxytriptolide | A/C | Vina | 10 | −7.085 ± 0.021 | −7.106至−7.048 | 0.048/0.087 |
| 16-hydroxytriptolide | A/C | Vinardo | 10 | −4.854 ± 0.020 | −4.870至−4.800 | 0.025/0.104 |
| 16-hydroxytriptolide | B/D | Vina | 10 | −7.577 ± 0.013 | −7.592至−7.552 | 0.037/0.054 |
| 16-hydroxytriptolide | B/D | Vinardo | 10 | −5.296 ± 0.011 | −5.309至−5.280 | 0.039/0.094 |

16-hydroxytriptolide在全部四个受体/评分组合中均高度收敛，最大种子间姿势RMSD为0.104 Å。该结果支持局部最优姿势具有计算重复性。

## 两种化合物的姿势和接触重叠

| 受体/评分 | 共同骨架直接RMSD (Å) | 姿势质心距离 (Å) | 接触Jaccard | 共同接触数 |
|---|---:|---:|---:|---:|
| A/C Vina | 0.618 | 0.204 | 0.90 | 9 |
| A/C Vinardo | 0.131 | 0.259 | 1.00 | 10 |
| B/D Vina | 0.345 | 0.316 | 1.00 | 11 |
| B/D Vinardo | 0.190 | 0.096 | 1.00 | 11 |

共同骨架在同一受体坐标系内的直接RMSD为0.13–0.62 Å，接触残基Jaccard为0.90–1.00。B/D代表姿势共同接触B链Pro43、Gln44、Tyr45、Ala82、Lys85、Ser86、Glu89以及D链Leu2、Glu6、Asn10、Ile13。这是本次分析中最强的结构桥接证据。

## 对接分数不能支持等效亲和力

16-hydroxytriptolide相对TPL的代表姿势分数差为：A/C Vina +0.159、A/C Vinardo +0.405、B/D Vina +0.219、B/D Vinardo −0.047 kcal/mol。方向并非四个组合完全一致，而且对接分数不是实验结合自由能。因此不能据此声称两种化合物具有相同亲和力、相同靶点占有率或相同细胞效应。

## 能否替换原文

### 可以替换的内容

可以把过去缺乏证据的“TPL因属于同类化合物而用于后续研究”升级为：

> Matched comparative docking independently placed triptolide and 16-hydroxytriptolide in the same predicted S100A8 homodimer-interface region, with highly similar common-scaffold orientations and 90–100% overlap of representative contact-residue sets.

该结论能够支持“使用TPL进行探索性结构和实验跟进具有结构合理性”。

### 仍不能替换或声称的内容

- TPL实验验证了16-hydroxytriptolide；
- 两种化合物功能等效或具有相同药理效力；
- 两种化合物均已被证明直接结合S100A8；
- TPL的MD结果适用于16-hydroxytriptolide；
- 对接分数等同于实验Kd或结合自由能。

## 最终判断

本次比较对接足以补上“结构类似物替换缺少平行对接证据”这一缺口，并可用于替换原稿中证据不足的class-level结构推断。但它只能建立共享预测口袋的计算桥梁。所有细胞实验、CETSA和分子动力学结论必须继续严格限定为TPL；16-hydroxytriptolide仍应被描述为计算优先候选，而不是被TPL实验验证的化合物。
