# 名字：NH4NO3nice
# 日期：2022年05月05日

import f90nml
import struct
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker


def find(lon, lat, four_side):
    """four_side = [wlon, elon, slat, nlat]
    本函数目的：在lon,lat中找到一个four_side1，其范围包含four_side"""
    lon, lat = np.array(lon), np.array(lat)
    four_side1 = [0, 0, 0, 0]  # 定义一个列表，存放找到的值
    # 找到左边界
    for i in range(len(lon)):
        if lon[i] > four_side[0]:
            four_side1[0] = lon[i - 1]
            break

    # 找到右边界
    for i in range(len(lon)):
        if lon[i] >= four_side[1]:
            four_side1[1] = lon[i]
            break
    # 找到下边界
    for i in range(len(lat)):
        if lat[i] <= four_side[2]:
            four_side1[2] = lat[i]
            break
    # 找到上边界
    for i in range(len(lat)):
        if lat[i] < four_side[3]:
            four_side1[3] = lat[i - 1]
            break
    return four_side1


with open(r'C:\Users\NH4NO3nice\Desktop\PE_module\Postprocess\uvz.ctl', 'r', encoding='utf-8') as f:
    for i in range(12):
        temp = f.readline()
        if i == 3:  # XDEF
            temp = temp.split()  # 以空格分隔，生成列表
            nx, wlon, dx = temp[1], temp[3], temp[4]
            nx, wlon, dx = float(nx), float(wlon), float(dx)
        if i == 4:  # YDEF
            temp = temp.split()  # 以空格分隔，生成列表
            ny, slat, dy = temp[1], temp[3], temp[4]
            ny, slat, dy = float(ny), float(slat), float(dy)

elon, nlat = wlon + (nx - 1) * dx, slat + (ny - 1) * dy
lon = np.arange(wlon, elon + 1, dx)
lat = np.arange(slat, nlat + 1, dy)
# lat = np.arange(nlat, slat - 1, -1*dy)  # 不能这么写，否则grads画图会出错，上下颠倒
# namelist文件地址
namelist_file = r'C:\Users\NH4NO3nice\Desktop\PE_module\namelist'
# 读取该namelist文件内容
namelist = f90nml.read(namelist_file)
# 提取namelist文件里面的control变量
data_path = namelist['pre_processing']['data_path']
data_name = namelist['pre_processing']['data_name']
start_time = namelist['control']['start_time']
temp = ['ua', 'va', 'za']
for i in range(3):
    data = xr.open_dataset(data_path[i])
    xlon, ylat = data['lon'], data['lat']  # 提取物理量经纬度
    four_side = find(xlon, ylat, [wlon, elon, slat, nlat])  # 找到合适的经纬度，要求原网格必须包含指定网格
    ph = data[data_name[i]].loc[start_time, 500, four_side[3]:four_side[2], four_side[0]:four_side[1]]  # 提取初始时间500hPa区域的物理量
    xlon = xlon.loc[four_side[0]:four_side[1]]
    ylat = ylat.loc[four_side[3]:four_side[2]]
    # 插值
    ph_grid25 = ph.interp(lon=lon, lat=lat, method='linear')

    # 存插值后的数据
    input_path = namelist['control']['input_path']
    with open(input_path + temp[i] + '.grd', 'wb') as f:
        for y in lat:  # 由南向北存
            for x in lon:  # 由西向东存
                s = struct.pack('f', float(ph_grid25.loc[y, x]))
                f.write(s)

    # 绘制插值前后的对比图
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    font = {'family': 'Calibri',
            'style': 'italic',  # 修改倾斜程度
            # 'weight': 'normal',  # 修改粗细
            'color': 'black',  # 颜色
            'size': 16,  # 字体大小
            }  # 设置xlabel、title、text等字体设置

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(15, 5), dpi=300)
    for j in range(2):
        ax = fig.add_subplot(1, 2, j + 1, projection=proj)  # 创建画布
        ax.coastlines()  # 绘制海岸线
        extent = [wlon, elon, slat, nlat]
        ax.set_extent(extent, crs=ccrs.PlateCarree())  # 调整地图经纬度范围,crs很重要
        # 设置feature
        ax.add_feature(cfeature.LAND)  # 添加陆地
        ax.add_feature(cfeature.COASTLINE, lw=0.3)  # 添加海岸线
        ax.add_feature(cfeature.RIVERS, lw=0.3)  # 添加河流
        if j == 0:
            ca = ax.contourf(xlon, ylat, ph, 15, transform=ccrs.PlateCarree(), cmap='RdBu_r')
            ax.contour(xlon, ylat, ph, 15, transform=ccrs.PlateCarree(), colors='k', linewidths=0.5)
            ax.set_title('(a) Before the interpolation', fontdict=font)
        if j == 1:
            ca = ax.contourf(lon, lat, ph_grid25, 15, transform=ccrs.PlateCarree(), cmap='RdBu_r')
            ax.contour(lon, lat, ph_grid25, 15, transform=ccrs.PlateCarree(), colors='k', linewidths=0.5)
            ax.set_title('(b) After the interpolation', fontdict=font)
        # 设置坐标
        ax.set_xticks(np.arange((wlon // 10 + 1) * 10, (elon // 10) * 10, 10), crs=proj)
        ax.set_yticks(np.arange((slat // 10 + 1) * 10, (nlat // 10) * 10, 5), crs=proj)
        lon_formatter = cticker.LongitudeFormatter(zero_direction_label=False)
        lat_formatter = cticker.LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        # 设置刻度字体大小
        ax.tick_params(labelsize=13)
        # 设置等值线上的数字标记
        ax.clabel(ca, fmt='%2.0f', fontsize=6, colors="k")
        # 添加色标
        cb = plt.colorbar(ca, orientation='horizontal', format='%d')
        cb.ax.tick_params(labelsize=13)  # 设置色标刻度字体大小
    plt.savefig('C:\\Users\\NH4NO3nice\\Desktop\\' + temp[i] + '.png', bbox_inches='tight', pad_inches=0.0)
