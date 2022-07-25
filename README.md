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





























