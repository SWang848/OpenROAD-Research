# lgista — Local Graph ISTA 使用教程

`lgista` 是从 `lrf` 模块中提取出来的**本地图增量时序引擎**(Local Graph
Incremental STA)。它只做两件事:

- **计算(compute)**:围绕一个实例建一张局部时序图(`PtGraph`),在图内传播
  delay / slew / arrival / required;支持"虚拟换 cell"的 what-if 评估。
- **写回(write)**:把局部图算出的结果提交回全局 `sta::Graph`。

它**不**包含调度、sizing、Lagrangian 优化——遍历顺序和优化策略由"在它之上
构建的优化器"决定。本文以原项目示例电路 `test/chain.v` 为主线,一步步教你用,并用
**「期望结论 / 怎么运行 / 实际结论」**的方式给出可复现的验证。

> TaiWei 集成说明：lgista 现在是 `src/3ddp` 的私有依赖，默认只构建 runtime
> library。迁移来源没有包含本文引用的 `test/chain.v` fixture，因此独立测试程序
> 默认关闭；恢复 fixture 后才可按本文方法重新验证。

---

## 0. 构建

`src/3ddp/CMakeLists.txt` 通过私有子目录接入 lgista，产出：

- `tdp_lgista_lib` —— TDP 使用的时序库
- `lgista_test` —— 可选命令行验证程序（需显式启用
  `TDP_LGISTA_BUILD_TEST` 并补齐 fixture）

```bash
cd build
cmake .. -DTDP_LGISTA_BUILD_TEST=ON
cmake --build . --target tdp_lgista_lib lgista_test --parallel 8
```

---

## 1. 示例电路 `test/chain.v`(Nangate45)

```
in ─▶ u1(INV) ─n1─▶ u2(BUF) ─n2─┬─▶ u3(INV) ─n3─┐
                                 └─▶ u4(INV) ─n4─┴─▶ u5(NAND2) ─▶ out
```

5 个实例,拓扑层级 `u1 < u2 < {u3,u4} < u5`,含一个**分叉**(u2→u3,u4)和一个
**重汇聚**(n3,n4→u5)——足够覆盖"多级传播 + 分支 + 汇聚"。

命令行总览:

```
lgista_test <design.v|.odb> <lib...> [--lef <f.lef>]... [--top <name>]
            [--inst <name>] [--test t1|t2a|t2b|all]
```

> `.v` 走 LEF+Verilog 流程(需 `--lef`/`--top`);`.odb` 直接读数据库。无需
> `openroad` 二进制。

---

## 2. 场景一:查看一个实例的本地时序(compute)

**我期望得到什么结论**
给定一个实例,引擎能围绕它建出局部图,并算出图内每个 pin 的 arrival/required
(以及 slew)。这是 ISTA 最基本的"计算"能力。

**我怎么运行**
```bash
cd <openroad-source>
./build/src/3ddp/src/lgista/lgista_test src/3ddp/src/lgista/test/chain.v \
    test/Nangate45/Nangate45_typ.lib \
    --lef test/Nangate45/Nangate45.lef --top chain --inst u5
```

**结论如何**(节选)
```
=== local timing for instance 'u5' ===
LocalSTA::printLocalTiming: Vertex u3/ZN
  Path: u3/ZN ^ default/max 2, arrival = 40.960055, required = -7.631066
  Path: u3/ZN v default/max 3, arrival = 38.564419, required = -7.862778
LocalSTA::printLocalTiming: Vertex u5/A1
  Path: u5/A1 ^ default/max 2, arrival = 40.960055, ...
...
```
局部图正确包含了 u5 的扇入(u3、u4)与其 pin,arrival 沿级递增(单位 ps)。

> 说明:示例没有 SDC(时钟/约束),所以 `min` corner 的 `required` 显示为未约束
> 的哨兵大数,`max` corner 的 required 也只是无约束下的占位值。**要得到有意义的
> required/slack,请在设计上加时钟约束**;delay/slew/arrival 不需要约束即有意义。

---

## 3. 场景二:虚拟换 cell 的局部结果 == 全局真实换 cell(T1)

**我期望得到什么结论**
在局部图里对某实例做 `virtualReplaceCell(C')` 并重算后,该实例**输出 pin 的
slew 和门 arc delay**,应与"对全局图真实 `replaceCell(C')` 后 `updateTiming`"读到
的**完全一致**。这验证 what-if 评估的保真度。

**我怎么运行**
```bash
./build/src/3ddp/src/lgista/lgista_test src/3ddp/src/lgista/test/chain.v \
    test/Nangate45/Nangate45_typ.lib \
    --lef test/Nangate45/Nangate45.lef --top chain --test t1
```

**结论如何**
```
=== T1: virtualReplaceCell vs global replaceCell ===
    tested 1 instances (skipped 4: input not stable under swap = second-order boundary effect)
[T1] 10 checks, 0 fail, max_abs=0, max_rel=0  -> PASS
```
对"换 cell 后输入 slew 不变"的实例(此处为 PI 驱动的 `u1`),local 与 global
**逐位一致(max_abs=0)**。

> 为什么跳过 4 个:换一个 cell 会改变它的*输入 pin 电容*,这会经上游驱动反馈回来
> 改变它*自己的输入 slew*(二阶效应)。**单实例局部图刻意不建模这种上游反馈**,
> 所以只对"换后输入稳定"的实例做严格逐位对比(测试会自动判定并跳过其余)。这是
> 局部图的设计取舍,不是误差。

---

## 4. 场景三:写回不破坏全局时序(T2-A,恒等性)

**我期望得到什么结论**
在**不换 cell** 的前提下,按拓扑序逐个实例「建局部图→重算→写回全局」后,全局
timing 应与写回前**完全相同**——证明写回这条路径本身没有副作用。

**我怎么运行**
```bash
./build/src/3ddp/src/lgista/lgista_test src/3ddp/src/lgista/test/chain.v \
    test/Nangate45/Nangate45_typ.lib \
    --lef test/Nangate45/Nangate45.lef --top chain --test t2a
```

**结论如何**
```
=== T2-A: write-back identity (no change) ===
[T2-A] 22 checks, 0 fail, max_abs=0, max_rel=0  -> PASS
```
写回前后全局 slew **逐位一致**。顺带证明:local 引擎(`makeLocalDelayCalc`)与
全局 `updateTiming` 在未改动时计算结果完全相同(同源)。

---

## 5. 场景四:多级换 cell + 拓扑序写回 == 全局完整重算(T2-B)

**我期望得到什么结论**
把**所有**实例换成等价 cell 后,按拓扑序逐级「建局部图→重算→写回」,最终全局
timing 应与"对同一最终网表做一次完整全局 `updateTiming`"**完全一致**。这验证
多级跨级传播下的写回正确性。

**我怎么运行**
```bash
./build/src/3ddp/src/lgista/lgista_test src/3ddp/src/lgista/test/chain.v \
    test/Nangate45/Nangate45_typ.lib \
    --lef test/Nangate45/Nangate45.lef --top chain --test t2b
```

**结论如何**
```
=== T2-B: replace-all + topo write-back vs global ===
    replaced 5 instances
[T2-B] 22 checks, 0 fail, max_abs=0, max_rel=0  -> PASS
```
增量写回结果与全局完整重算**逐位一致**。

> 方法学要点(很重要):必须**先把所有 cell 换完**(负载电容进入终态),**再**按
> 拓扑序逐级建图+写回。若把"换 cell"与"重算"交错进行,上游实例会看到下游负载的
> *旧*电容,导致与全局重算系统性不符。这条同时说明:基于局部图的增量更新,负载
> 电容的改变必须在重算前全部就位(否则需反向序或迭代收敛)。

一次跑全部:`--test all`(t1 + t2a + t2b),全绿即 `=== test 'all' -> PASS ===`。

---

## 6. C++ API 用法

ISTA 暴露两组原语,**不**管遍历顺序(那是优化器的事)。

**门面 `LocalGraphIsta`**(最省事):
```cpp
#include "lgista/LocalGraphIsta.hh"

lgista::LocalGraphIsta ista(dbSta);          // dbSta:已加载设计+liberty
ista.prepare();                              // 全局 timing 就绪(updateTiming)
lgista::PtGraph* g = ista.analyzeInstance("u5");  // 建局部图 + 算 delay/arrival/required
ista.commitToGlobal(g);                      // 何时/按什么序写回,由你决定
```

**直接用 `LocalSta` / `PtGraph`**(更灵活,做 what-if):
```cpp
lgista::LocalSta L(dbSta);  L.copyState(dbSta);
sta::ArcDelayCalc* calc = lgista::makeLocalDelayCalc(dbSta);

lgista::PtGraph* pt = L.makePtGraph(inst, /*update_timing_first=*/true);
L.findLocalDelays(pt, calc);  L.findLocalArrivals(pt);  L.findLocalRequireds(pt);

L.virtualReplaceCell(pt, new_cell);          // what-if 换 cell
L.findLocalDelays(pt, calc);                 // 重算
// 读结果:pt->slew(ptVertex, rf, ap) / pt->arcDelay(ptEdge, arc, ap)
lgista::updateTimingFromPtGraph(pt);         // 写回全局(写 API)
```

要"在它之上搭优化器",就是:你提供遍历顺序与决策,循环调用上面的
计算 + 写回原语。

---

## 7. 在你自己的设计上跑

- **LEF + Verilog**:`lgista_test design.v lib1.lib [lib2.lib] --lef tech.lef [--lef cells.lef] --top <top>`
- **OpenDB**:`lgista_test design.odb lib1.lib [lib2.lib]`

`src/LgistaMain.cc` 是最小加载脚手架,可按需扩展(LEF/DEF、SDC、寄生、scene)。

---

## 8. 重要前提(否则结论对不上)

1. **延迟引擎一致**:local 用 `makeLocalDelayCalc`(DMP/ElmoreCeff);全局
   `updateTiming` 须为同一套(T2-A 的逐位一致即此证明)。
2. **单一 analysis point**:示例 liberty 按 `MinMaxAll::max()` 读,读数 ap 取
   `scene->dcalcAnalysisPtIndex(MinMax::max())`,与 `makePtGraph` 设的
   `(default scene, max)` 对齐。多 ap 下 PtGraph 内部存储与全局 DcalcAPIndex 会错位。
3. **oracle 读全局图**:对比时直接读 `sta::Graph`(`graph->slew/arcDelay`);用新建
   PtGraph 读全局不可靠(`copyInfoFromVertex` 不镜像引擎 ap 上的 slew)。
4. **arc delay 按 arc 取值**:`PtGraph::arcDelay(edge,arc,ap)` 与
   `Graph::arcDelay(edge,arc,ap)` 两侧按同一 `arc`;不要按数组下标(顺序不同)。
5. **二阶边界效应**:换 cell 改输入电容→反馈改自身输入 slew,单实例局部图不建模
   (见场景二)。增量全图更新需"先换完再重算"(见场景四)。
6. **调试**:`LGISTA_DEBUG=1` 打印 T1 首个实例的逐 pin / 逐 arc 对比值。

---

## 9. 文件清单

| 文件 | 作用 |
|------|------|
| `PtGraph.{hh,cc}` | 本地时序图:`PtVertex`/`PtEdge`、拓扑、路径、写回 |
| `LocalSta.{hh,cc}` | 本地 STA 引擎:建图、delay/arrival/required、slack、ERC、虚拟换 cell |
| `LocalSearch.{hh,cc}`、`LocalPath.cc` | arrival/required 传播、路径回溯 |
| `LocalParasitics.{hh,cc}`、`LocalReduceParasitic.{hh,cc}` | 寄生网络 + Pi 模型规约 |
| `LocalDmpDelayCalc.{hh,cc}`、`LocalElmoreCeffDelayCalc.{hh,cc}` | 延迟计算后端;`makeLocalDelayCalc()` 工厂 |
| `PtPiElmore.{hh,cc}`、`PtElmoreCeff.{hh,cc}` | 每驱动寄生模型 |
| `ViolationCheck.cc` | slew/cap 限值(ERC)检查 |
| `LrfSearchPred.{hh,cc}` | `SearchMEEPred` 搜索谓词(从 TaskArranger.hh 抽出) |
| `IstaWriteback.{hh,cc}` | 写 API:`updateTimingFromPtGraph()`(本地→全局) |
| `include/lgista/LocalGraphIsta.hh`、`src/LocalGraphIsta.cc` | C++ 门面 |
| `src/LgistaMain.cc` | 命令行驱动(加载 + analyze + `--test`) |
| `src/LgistaTests.cc` | T1/T2-A/T2-B 验证实现 |
| `test/chain.v` | 示例电路 fixture |

> 命名空间用独立的 `lgista`(源自 `lrf`,部分注释仍写 lrf):因为 `dbSta` 引用
> `lrf::IncreSta`,链接 lgista 的可执行文件同时会链 `lrf_lib`,独立命名空间避免
> 同名符号冲突。唯一在本目录外的改动:`src/sta/parasitics/ConcreteParasitics.hh`
> 追加 `friend class lgista::LocalParasitics;`(`LocalParasitics` 经基类指针访问其
> 受保护成员所必需)。
