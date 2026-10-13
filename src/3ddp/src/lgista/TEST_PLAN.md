# lgista 测试方案 (Local Graph ISTA 正确性验证)

> TaiWei 集成说明：迁移来源保留了测试 driver，但没有包含本文引用的
> `test/chain.v` fixture。`TDP_LGISTA_BUILD_TEST` 因此默认关闭；本文件保留为测试
> 方法和历史结果记录，不能替代当前版本重新编译和实跑。

本方案验证两件事:

- **T1 — virtualReplaceCell 等价性**:在 local graph 内做 `virtualReplaceCell`
  并重算后,其 **delay / slew** 与"对全局图真实 `replaceCell` 同一个 cell 后"
  在 *local graph 覆盖范围内的相同 pin/arc* 上读到的值一致。
- **T2 — write-back 正确性**:构造一个**多级电路**,按**拓扑序**遍历它的全部
  instance,(可选地换 cell)→ 重算 local timing → `commitToGlobal` 写回;最后
  检查写回后的全局 timing 与"对同样改动做一次完整全局重算"的结果一致。

> 关键参照系(oracle)始终是 **OpenSTA 的全局引擎**(`replaceCell` +
> `updateTiming`/`findArrivals`/`findRequireds`),因为 lgista 的语义就是
> "在 local graph 上复现全局 STA 在该局部的结果"。

---

## 0. 通用前置与口径统一(必须先做,否则比较无意义)

这几条直接决定 T1/T2 能不能"对得上",是方案里最关键的部分:

1. **延迟计算引擎一致**。local 侧用 `lgista::makeLocalDelayCalc()`(DMP/
   ElmoreCeff);global 侧必须配同一套算法。
   - 做法:全局 `set_delay_calculator` 选择与 local 相同的 calc(DMP +
     Elmore/Ceff),或在 C++ 里把 `sta->arcDelayCalc()` 设为同一实现。
   - 若两侧引擎不同,delay/slew 会有系统性偏差 → 比较必然失败。这是 #1 风险。
2. **单一 Scene / 单一 analysis point**。固定一个 corner/scene、一个
   `DcalcAnalysisPt`(取 `scene->...->dcalcAnalysisPt()`),min 或 max 二选一,
   两侧用同一个 `ap_index` 读数。
3. **寄生来源一致**。local graph 在 swap 后会重算 Pi/Elmore 寄生
   (`recomputePtParasitics` / `recomputeSinglePtParasitic`);global 侧的寄生
   必须来自同一来源(同一份估算或同一份 SPEF)。
   - 做法:`LocalSta::setParasiticsEst(est)` 用与全局相同的
     `est::EstimateParasitics`;或两侧都读同一份 SPEF。
4. **EquivCells 已建立**。换 cell 需要等价库:`sta->makeEquivCells(...)`
   (dbSta 通常已建)。候选 cell 用 `sta->equivCells(cell)` 取。
5. **输入边界稳定**。换 *单个* instance I 不影响其上游驱动,因此 I 的输入 pin
   slew 在 local 和 global 两侧相同——这是 T1 能逐点对齐的前提。

**通用读数辅助**(两侧都按同一 `(pin, rf, ap_index)` 读):

```cpp
// global
sta::Vertex* gv = sta->graph()->pinDrvrVertex(pin);   // 或 pinLoadVertex
sta::Slew g_slew = sta->graph()->slew(gv, rf, ap_index);
sta::ArcDelay g_arc = sta->graph()->arcDelay(edge, arc, ap_index);

// local
const lgista::PtVertex& pv = pt_graph->pinToPtVertex(pin);
sta::Slew l_slew = pt_graph->slew(pv, rf, ap_index);
// arc: 用 PtVertexOutEdgeIterator 遍历 pv 的出边,pt_graph->arcDelay(pt_edge, arc, ap_index)
//      wire 段用 pt_graph->wireArcDelay(pt_edge, rf, ap_index)
```

---

## 1. 测试电路(fixtures)

用现成 liberty(如 `Nangate45` / `sky130` / `asap7`,挑有多种 drive 的反相器/
缓冲器)。提供 3 个由小到大的电路;统一加 RC(估算或 SPEF)与时钟约束。

### C1 — 单门(T1 最小用例)
```
PI --(wire)--> U1.A   U1(INVx1) --(wire)--> [L1.A, L2.A]   (两个 fanout load)
PI 上 set_input_delay,U1 输出接 set_load / 两个 load pin。
```
目的:local graph 恰好覆盖 {U1 的输入 pin、U1 输出 pin、两个 load pin、各 wire 段}。

### C2 — 多级组合链 + 重汇聚(T1 加强 + T2 主力)
```
PI -> U1(INV) -> n1 -> U2(BUF) -> n2 ┬-> U3(INV) -> n3 -> U5(BUF) -> PO
                                     └-> U4(INV) -> n4 ----------------/ (n4 汇入 U5 另一输入)
```
- 4~5 级,带一个分叉(U2 同时驱动 U3、U4)和一个重汇聚(U5 两输入)。
- 拓扑序:U1 < U2 < {U3,U4} < U5。验证"上游写回后下游 local 图读到更新值"。

### C3 — 时序版(T2 的 required/slack 验证,可选)
```
FF_launch.Q -> [C2 的组合网络] -> FF_capture.D ;  CLK 接两个 FF。
```
- 用于检查 arrival/required/slack 的写回(含 setup check),验证多级传播下
  endpoint slack 与全局一致。FF 模式较复杂,列为可选增强。

> 电路可用 Verilog+LEF/DEF 或程序化建 odb;关键是**固定、可复现**,并带 RC。

---

## 2. T1:virtualReplaceCell 局部等价性

**对象**:C1(以及 C2 里逐个 instance 重复)。对实例 I、候选等价 cell C'。

**步骤**
1. `prepare()` 让全局 timing 最新。
2. **local 测量**:
   - `pt = makePtGraph(I, update_timing_first=true)`
   - `findLocalDelays(pt, calc); findLocalArrivals(pt); findLocalRequireds(pt)`
   - 记录 baseline(可选,用于自检 baseline==global 当前值)。
   - `virtualReplaceCell(pt, C')` → `findLocalDelays/Arrivals/Requireds` 重算。
   - 采集 **local 集合 L**:I 的输入/输出 pin slew、I 内部 timing arc delay、
     I 输出到各 load 的 wire arc delay、各 load pin slew。
3. **global oracle**:
   - `sta->replaceCell(I, C')`(真实换);`sta->updateTiming(true)`。
   - 在**相同的 pin/arc**上采集 **global 集合 G**。
   - **还原**:`sta->replaceCell(I, C_orig); sta->updateTiming(true)`(保证后续测试干净)。
4. **断言**:对 L 与 G 的每一项,`|local - global| ≤ tol`(见 §4)。

**比较点(local graph 覆盖范围内)**
| 量 | local 读法 | global 读法 |
|---|---|---|
| 输入 pin slew | `pt->slew(pinToPtVertex(in_pin), rf, ap)` | `graph->slew(loadVertex(in_pin), rf, ap)` |
| 输出 pin slew | 同上(out_pin) | 同上(drvrVertex) |
| 门内 arc delay | `pt->arcDelay(pt_edge, arc, ap)` | `graph->arcDelay(g_edge, arc, ap)` |
| wire arc delay | `pt->wireArcDelay(pt_edge, rf, ap)` | `graph->arcDelay(wire_edge, rf, ap)` |
| load pin slew | `pt->slew(pinToPtVertex(load_pin), rf, ap)` | `graph->slew(loadVertex(load_pin), rf, ap)` |

**不纳入比较**:load 的**下游**(local graph 在 load pin 处截断,不建模其后级)。

**通过判据**:L 中所有点全部在 tol 内;对 rise/fall 两个 transition 都成立。

---

## 3. T2:write-back 正确性(多级、拓扑序)

**对象**:C2(组合)主验;C3(时序)可选增强。

### 拓扑序的获取
不用 lgista 自带调度(ISTA 不管顺序)。测试自己排序:
- 用 `sta->ensureLevelized()` 后按 `Graph::level(vertex)` 升序;或对实例做一次
  Kahn 拓扑排序(按 fanin→fanout)。得到实例序列 `[U1, U2, U3/U4, U5]`。

### T2-A:恒等写回(no-change identity)—— 先过这关
目的:在**不换 cell** 的前提下,验证"逐实例建图→重算→写回"不破坏全局 timing。

1. `prepare()`;记录**全局基线** B0 = 所有 pin 的 {slew, arrival, required}。
2. 按拓扑序对每个实例 U:
   - `pt = analyzeInstance(U, print=false)`(建图+算 delay/arrival/required)
   - `commitToGlobal(pt)`(写回 slew + arrival/required)
3. 记录写回后的全局 B1。
4. **断言**:`B1 == B0`(逐点 tol 内)。写回本质是"把本地重算结果覆盖回去",
   未改动设计时必须与原全局一致 → 这是写回路径无副作用的纯净性检查。

### T2-B:换 cell + 多级传播(核心)
1. `prepare()`。
2. 选一组等价替换 `{U_i -> C'_i}`(例如全部升一档 drive)。
3. **lgista 路径**(按拓扑序):对每个 U_i:
   - `pt = makePtGraph(U_i, update_timing_first=true)`
   - `localSta()->virtualReplaceCell(pt, C'_i)`
     —— 注意:若希望写回真正反映"已换 cell",此处应在**全局也真实换上 C'_i**
     后再 commit;见下方两种口径。
   - `findLocalDelays/Arrivals/Requireds(pt)`;`commitToGlobal(pt)`
   - 因为按拓扑序,U_{i+1} 的 local 图在其之后建立,会从全局读到 U_i 刚写回的
     更新 slew/arrival → 实现跨级传播。
4. **global oracle**:在干净副本/重置后,对同一组 `{U_i -> C'_i}` 执行真实
   `replaceCell` 全部 + `sta->updateTiming(true); findRequireds()`,得到 G*。
5. **断言**:lgista 写回后的全局 {slew, arrival, required, endpoint slack}
   与 G* 逐点 tol 内一致。重点比较 **PO / FF.D 端点的 arrival/required/slack**,
   以及每级中间 pin 的 slew/arrival。

**两种口径(择一并写清)**
- **(i) 纯 what-if 写回**:全局 netlist 不真换 cell,仅把 local(virtual 换 C'
  后)的 timing 写回全局图。oracle 也必须是"netlist 不变、仅 timing 反映 C'"
  ——这要求 oracle 也走 replaceCell(因为全局 updateTiming 依赖真实 netlist)。
  因此**纯 what-if 写回无法用全局重算直接对照**,只适合做 T1 式逐实例对照。
- **(ii) 真实换 + 写回(推荐用于 T2-B)**:lgista 路径里对每个 U_i **也真实
  `replaceCell`**(改 netlist),再用 local 图算 + 写回;oracle 用"全部真实
  replaceCell + 一次全局 updateTiming"。两者 netlist 相同,可直接对照。
  → T2-B 用口径 (ii)。

**通过判据**:端点与中间节点的 {slew, arrival, required} 全部 tol 内;
拓扑序保证下,逐级值随上游更新而正确变化(可加一项:打乱顺序应导致中间过程
不一致,反证拓扑序的必要性 —— 作为加分项)。

---

## 4. 容差(tolerance)与依据

- 若 §0 的引擎/寄生/scene 完全对齐,期望**接近逐位一致**:
  - 绝对:`|Δdelay|, |Δslew| ≤ 1e-13 ~ 1e-12 s`
  - 相对:`≤ 1e-4`(取 `abs<=1e-12 || rel<=1e-4` 任一满足)
- arrival/required 误差会沿级累积,放宽到相对 `1e-3`。
- 若 local 与 global 引擎/寄生不完全一致(例如 global 用 lumped、local 用
  ElmoreCeff),改为**趋势+宽容差**验证(相对 1~5%),并在报告里标注偏差来源。
  —— 但首选是把口径对齐到能用紧容差。

---

## 5. 已知风险 / 易错点(评审重点)

1. **延迟引擎不一致**(最常见失败原因)→ §0.1 必须落实。
2. **寄生不一致**:swap 后 local 重算寄生,global 若未同步会偏;用同一 est 源。
3. **virtualReplaceCell 仅对等价 cell 有效**:端口必须一致;先建 equivCells。
4. **多 scene/mode**:新 base 引入 Scene/Mode,务必固定单 scene、单 ap_index。
5. **还原**:T1 的 global oracle 会改 netlist,测完必须 replaceCell 还原,
   否则污染后续用例。
6. **截断边界**:T1 不要比较 load 下游;T2 注意 PI/PO 边界约束(input_delay /
   output load)两侧一致。
7. **写回范围**:`updateTimingFromPtGraph` 只写 RefInput/RefOutput/SiblingLoad,
   virtual buffer pin 不写——T2 不要期望 virtual buffer 在全局出现。

---

## 6. 实现与运行

- 形式:扩展 `lgista_test`,或新增一个 gtest(`test/lgista/`),加载固定 fixture。
- 入口:复用 `LgistaMain` 的加载骨架(.odb + .lib + scene),再加：
  - `t1_virtual_replace(sta, inst_name, target_cell)`
  - `t2a_identity_writeback(sta)` / `t2b_replace_writeback(sta, replacements)`
- 断言:用一个 `expectClose(name, a, b, tol)` 收集所有偏差,末尾汇总
  PASS/FAIL 并打印最大偏差点,便于定位。
- 运行:
  ```
  lgista_test <C2.odb> <lib.lib> --test t2b
  ```
  (需在 `LgistaMain` 加一个 `--test` 分发;当前 main 仅做 analyze 演示。)

---

## 7. 验收清单

- [ ] §0 口径(引擎/scene/寄生/equivCells)全部对齐并在日志打印确认。
- [ ] T1:C1 + C2 每个组合 instance,rise/fall、delay/slew 全部 tol 内。
- [ ] T2-A:恒等写回,B1==B0(tol 内)。
- [ ] T2-B:多级 C2,真实换 + 拓扑序写回,端点与中间节点对齐全局重算 G*。
- [ ] (可选) C3:FF 版 endpoint slack 对齐。
- [ ] (加分) 打乱拓扑序 → 中间结果出现可观偏差(反证顺序必要性)。
- [ ] 报告输出最大偏差点 + 是否 PASS。

---

## 8. 实现与实跑结果(已落地)

实现：`src/3ddp/src/lgista/src/LgistaTests.{hh,cc}`，入口
`lgista_test ... --test t1|t2a|t2b|all`。
Fixture：`src/3ddp/src/lgista/test/chain.v`（当前迁移来源未包含）。

运行:
```
lgista_test src/3ddp/src/lgista/test/chain.v test/Nangate45/Nangate45_typ.lib \
    --lef test/Nangate45/Nangate45.lef --top chain --test all
```
结果:**T1 / T2-A / T2-B 全部 PASS(max_abs=0,逐位一致)**。

### 实跑中确认/修正的关键点(评审重点)

1. **延迟引擎确一致**:T2-A 恒等写回 `max_abs=0`(逐位),证明 local 的
   `makeLocalDelayCalc` 与全局 `updateTiming` 在未改动时完全一致。

2. **单一 analysis point**:liberty 必须按 `MinMaxAll::max()` 读 → 单个 dcalc ap。
   若用 `all()`,会有 min/max 两个 ap,而 PtGraph 内部 slew 存储与全局
   DcalcAPIndex 不对齐(用全局 index 读 PtGraph 会取到错位/未初始化值)。
   读数 ap 取 `scene->dcalcAnalysisPtIndex(MinMax::max())`,与
   `LocalSta::makePtGraph` 设的 `setAnalysisPoint(default, max)` 一致。

3. **oracle 不要用 PtGraph 读全局**:`makePtGraph` 的 `copyInfoFromVertex` 不能
   可靠镜像全局图在引擎 ap 上的 slew。**oracle 直接读全局 `sta::Graph`**
   (`graph->slew/arcDelay`)才是 ground truth。

4. **arc delay 按 arc 取值,不要按迭代位置**:local `PtEdge::arcDelays()` 数组按
   arc 自身 index 存,和 `TimingArcSet::arcs()` 迭代顺序不一致。用
   `PtGraph::arcDelay(ptEdge, arc, ap)` 与 `Graph::arcDelay(edge, arc, ap)`
   两侧按同一 `arc` 取值。

5. **T1 的"输入稳定"约束(重要修正)**:换单个 cell 会改变它的*输入 pin cap*,
   经上游驱动反馈回来改变它*自己的输入 slew*(二阶效应)。单实例 local graph
   **不**建模这个上游反馈。因此 T1 只对"换后输入 slew 不变"的实例(自检:换前后
   读全局输入 slew 比较)做严格输出对比 —— 对这些实例,local 的输出 slew + 门
   arc delay 与全局换 cell **逐位一致**(验证 virtualReplaceCell 正确)。受二阶
   反馈影响的实例被跳过并计数(是 local graph 的设计行为,不是 bug)。

6. **T2-B 必须"先换完所有 cell,再按拓扑序重算写回"**:若把 replaceCell 与重算
   交错,上游实例会看到下游负载的*旧*输入 cap → 与全局重算系统性不符。先全部
   replaceCell(负载终态确定)→ 再拓扑序 build+writeback,可与全局完整重算**逐位
   一致**。这条同时说明:基于 local graph 的增量更新,负载 cap 变化必须在重算前
   全部就位(否则需反向序或迭代收敛)。

### 调试
设环境变量 `LGISTA_DEBUG=1` 打印 T1 首个实例的逐 pin / 逐 arc 对比值。
