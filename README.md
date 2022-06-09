# The-Barotropic-Primitive-Equstion-Model

## `1.` 预报方程
`1.1` 地图投影坐标下的正压原始方程组

$ \frac{\partial u}{\partial t}=-m\left(u \frac{\partial u}{\partial x}+v \frac{\partial u}{\partial y}\right)+f^{*} v-m g \frac{\partial z}{\partial x} \\ $

$ \frac{\partial v}{\partial t}=-m\left(u \frac{\partial v}{\partial x}+v \frac{\partial v}{\partial y}\right)-f^{*} u-m g \frac{\partial z}{\partial y} \\ $

$ \frac{\partial z}{\partial t}=-m^{2}\left[u \frac{\partial}{\partial x}\left(\frac{z}{m}\right)+v \frac{\partial}{\partial y}\left(\frac{z}{m}\right)+\frac{z}{m}\left(\frac{\partial u}{\partial x}+\frac{\partial v}{\partial y}\right)\right] $

$ f^{*}=f+m^{2}\left[v \frac{\partial}{\partial x}\left(\frac{1}{m}\right)-u \frac{\partial y}{\partial}\left(\frac{1}{m}\right)\right] $

`1.2` 二次平流守恒格式有限差分方程  

