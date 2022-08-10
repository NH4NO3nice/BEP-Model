# The Barotropic Primitive Equation Model

## `1. 预报方程`

​		在地图投影坐标下的正压原始方程模式预报方程组

$$
\left\{\begin{array}{l}
\frac{\partial u}{\partial t} = -m\left(u \frac{\partial u}{\partial x}+v \frac{\partial u}{\partial y}\right)+f^{*} v-m g \frac{\partial z}{\partial x} \\
\frac{\partial v}{\partial t} = -m\left(u \frac{\partial v}{\partial x}+v \frac{\partial v}{\partial y}\right)-f^{*} u-m g \frac{\partial z}{\partial y} \\
\frac{\partial z}{\partial t} = -m^{2}\left[u \frac{\partial}{\partial x}\left(\frac{z}{m}\right)+v \frac{\partial}{\partial y}\left(\frac{z}{m}\right)+\frac{z}{m}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)\right] \\
f^{*} = f+m^{2}\left[v \frac{\partial}{\partial x}\left(\frac{1}{m}\right)-u \frac{\partial}{\partial y}\left(\frac{1}{m}\right)\right]
\end{array}\right.
$$

## `2. 二次平流守恒格式`

$$
\begin{array}{l}
\left\{\begin{array}{l}
\frac{\partial u_{i, j}}{\partial t}=-m_{i, j}\left(\overline{\bar{u}^{x} u_{x}}^{x}+\overline{\bar{v}^{y} u_{y}}^{y}+g \bar{z}_{x}^{x}\right)+\widetilde{f_{i j}^{*}} v_{i, j}=E_{i, j} \\
\frac{\partial v_{i, j}}{\partial t}=-m_{i, j}\left(\overline{\bar{u}^{x} v_{x}}^{x}+\overline{\bar{v}^{y} v_{y}}^{y}+g \bar{z}_{y}^{y}\right)-\widetilde{f_{i, j}^{*}} u_{i, j}=G_{i, j} \\
\frac{\partial z_{i, j}}{\partial t}=-m_{i, j}^{2}\left[\overline{\bar{u}^{x}\left(\frac{z}{m}\right)_{x}}^{x}+\overline{\bar{v}^{y}\left(\frac{z}{m}\right)_{y}}^{y}+\frac{z_{i, j}}{m_{i, j}}\left(\bar{u}_{x}^{x}+\bar{v}_{y}^{y}\right)\right]=H_{i, j}
\end{array}\right.\\
\widetilde{f_{i, j}^{*}}=f_{i, j}+u_{i, j} \bar{m}_{y}^{y}-v_{i, j} \bar{m}_{x}^{x}
\end{array}
$$

## `3. 初始条件`

​		将初始位势高度场、初始风场输入原始方程模式。理论分析和预报实践表明，由于观测的风场与高度场之间的不平衡，以及风场、高度场与模式之间的不协调，直接用观测的风场和高度场作为原始方程模式的初始值容易产生高频振荡，使数值积分变为不稳定。

​		因此，在应用原始方程模式制作数值天气预报之前，必须对初始资料加以处理，即所谓的资料初始化。资料初始化的方法有静力初始化、动力初始化和变分初始化。  

​		假定风场、气压场(位势高度场)之间满足某种平衡关系，根据这种平衡关系由其中的一个场来确定另一个场，这种处理初值的方法称为静力初始化。
$$
地转风初值:t=0时，
\begin{array}{l}
z_{i, j}=z_{i, j}^{0} \\
u_{i, j}=u_{i, j}^{0}=-\frac{m_{i, j} g}{f_{i, j}} \frac{\partial z_{i, j}^{0}}{\partial y} \\
v_{i, j}=v_{i, j}^{0}=\frac{m_{i, j} g}{f_{i, j}} \frac{\partial z_{i, j}^{0}}{\partial x}
\end{array}
$$

# `4. 边界条件`

​		对有限区域的预报，在边界上必须人为地给出水平侧边界条件。水平侧边界条件的给定方法有多种，例如固定边界条件、海绵边界条件等。
$$
固定边界条件：
\left.\frac{\partial u_{i, j}}{\partial t}\right|_{\beta}=\left.\frac{\partial v_{i, j}}{\partial t}\right|_{\beta}=\left.\frac{\partial z_{i, j}}{\partial t}\right|_{\beta}=0
$$


# `5. 时间积分方案`

​		设**F**为一矢量函数，它的分量为模式大气的预报变量**u**、**v** 和 **z/m**。于是，差分形式的正压原始方程组可简写为：
$$
\frac{\partial \boldsymbol{F}_{i, j}}{\partial t}=\boldsymbol{A}_{i, j} \boldsymbol{F}_{i, j}
$$
时间积分：可先采用欧拉-后差格式，可有效地抑制高频振荡，使数值积分稳定。
$$
\left\{\begin{array}{l}
\boldsymbol{F}_{i, j}^{* n+1}=\boldsymbol{F}_{i, j}^{n}+\Delta t \boldsymbol{A}_{i, j}^{n} \boldsymbol{F}_{i, j}^{n} \\
\boldsymbol{F}_{i, j}^{n+1}=\boldsymbol{F}_{i, j}^{n}+\Delta t \boldsymbol{A}_{i, j}^{*,{ }^{*}+1} \boldsymbol{F}_{i, j}^{* n+1}
\end{array}\right.
$$
随后。采用三步法起步的时间中央差格式：
$$
\begin{array}{l}
\boldsymbol{F}_{i, j}^{n+1 / 2}=\boldsymbol{F}_{i, j}^{n}+\frac{1}{2} \Delta t \boldsymbol{A}_{i, j}^{n} \boldsymbol{F}_{i, j}^{n} \\
\boldsymbol{F}_{i, j}^{n+1}=\boldsymbol{F}_{i, j}^{n}+\Delta t \boldsymbol{A}_{i, j}^{n+1 / 2} \boldsymbol{F}_{i, j}^{n+\frac{1}{2}} \\
\boldsymbol{F}_{i, j}^{n+2}=\boldsymbol{F}_{i, j}^{n}+2 \Delta t \boldsymbol{A}_{i, j}^{n+1} \boldsymbol{F}_{i, j}^{n+1}
\end{array}
$$

# `6. 数值求解过程中的其他技术问题`

① 在时间积分过程中，为阻尼高频振荡，抑制计算解的增长，可穿插进行时间平滑（3点平滑）：
$$
\widetilde{F_{i, j}^{n}}^{t}=(1-S) F_{i, j}^{n}+\frac{S}{2}\left(F_{i, j}^{n+1}+F_{i, j}^{n-1}\right)
$$
② 为滤除短波扰动，抑制非线性计算不稳定，可穿插进行空间平滑：

边界内第一圈格点 “边界点” 平滑：  9点平滑  

内点平滑：    5点平滑





















