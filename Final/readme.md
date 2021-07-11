# 使用方法

## Quick Start

直接运行文件夹下`./run_all.sh`即可。其中各参数含义如下：

```shell
time bin/FINAL testcases/scene02_office.txt output/scene03_final.bmp 400000
```

`testcases/scene02_office.txt`为输入文件，`output/scene03_final.bmp`为输出文件，`400000`为单次PPM的光子数乘4，即可得到报告中场景。

## 高级调整

### Path Tracing

只需在`src/main.cpp`中作如下改动即可，此时`400000`代表Path Tracing过程中每个点的光线条数。

```c++
25 Image img = pathTrace(parser, samples);
26 // Image img = photonMap(parser, samples);
```

### 参数调整

以下参数均位于`include/utils.hpp`

```c++
const double kNewtonIterEps = 1e-4; // 牛顿迭代时认为函数值绝对值小于该值即得解
const double kPlaneTextureScale = 5; // 平面上5 * 5代表一张纹理图，调大后平面上纹理图将变稀
const int kPPMTurn = 20; // PPM轮数
const int kMaxDepth = 8; // 光线追踪递归深度，调高后若存在B样条曲线折射将显著降低速度，调低后则会影响效果
const int kMaxIteration = 64; // 牛顿迭代单次循环最大次数
const int kDeriveShrink = 32; // 报告中提到“牛顿下山法”所除的系数x
const int kMaxRandom = 128; // 定义域中均匀取点的个数（对于本场景而言32即可达到应有效果）
```

可以通过调整这些参数，以达到更好的效果，以下为报告中参数值设定：

* 图1为：`kPPMTurn = 1`，脚本中数值为`16000000`
* 图2为：`kPPMTurn = 100`，脚本中数值为`400000`
* 图3为：`kPPMTurn = 100`，脚本中数值为`400000`，输入文件改为`testcases/scene02_depth.txt`
* 图4为：`kPPMTurn = 500`，脚本中数值为`80000`，输入文件改为`testcases/scene02_depth.txt`

