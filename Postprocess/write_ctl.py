# 名字：NH4NO3nice
# 日期：2022年05月05日

import f90nml
import numpy as np
from calendar import month_name


def normalize(x):
    """规范起止经纬度和网格距
       使个位数和小数位为：0.0，2.5，5.0，7.5"""
    list1 = np.arange(0, 10.1, 2.5)
    temp = x - x // 10 * 10  # 计算个位数和小数部分
    for i in range(len(list1)):
        if list1[i] > temp:
            break
    if (temp - list1[i - 1]) > (list1[i] - temp):
        num = x // 10 * 10 + list1[i]
    else:
        num = x // 10 * 10 + list1[i - 1]
    return num


# namelist文件地址
namelist_file = r'C:\Users\NH4NO3nice\Desktop\PE_module\namelist'
# 读取该namelist文件内容
namelist = f90nml.read(namelist_file)
# 提取namelist文件里面的control变量
start_time, all_time = namelist['control']['start_time'], namelist['control']['all_time']
output_interval = namelist['control']['output_interval']
nx, ny = namelist['control']['nx'], namelist['control']['ny']  # 提取区域网格数
d, dt = namelist['control']['d'], namelist['control']['dt']  # 提取区域网格距、积分步长
clat, clon = namelist['control']['clat'], namelist['control']['clon']  # 提取区域中心经纬度

"""确定区域网格经纬度范围"""
a = 6371  # 地球半径，单位km
d = d / 1000  # 将网格距的单位转换成km
# 确定经度左边界范围
sita = d * (nx - 1) / (a * np.cos(clat / 180 * np.pi))  # 计算纬度跨度rad
sita = sita * 180 / np.pi  # 将弧度转化为角度
wlon = clon - sita // 2  # 经度左边界
wlon = normalize(wlon)
list = np.arange(0.5, 10, 0.5)
# 确定经度间隔
min = 999
for i in list:
    if abs(i * (nx - 1) - sita) < min:
        min = abs(i * (nx - 1) - sita)
        dx = i

# 确定纬度下边界
sita = d * (ny - 1) / a  # 纬度跨度rad
sita = sita * 180 / np.pi  # 将弧度转化为角度
slat = clat - sita // 2  # 纬度下边界
slat = normalize(slat)
# 确定纬度间隔
min = 999
for i in list:
    if abs(i * (ny - 1) - sita) < min:
        min = abs(i * (ny - 1) - sita)
        dy = i

# 写ctl文件
output_path = namelist['control']['output_path']  # 读取输出数据地址
ctl_path = namelist['post_processing']['ctl_path']  # 读取ctl存放地址
month = month_name[1:13]  # 12个月
month = [i[:3] for i in month]
with open(ctl_path + 'uvz.ctl', 'w+', encoding='utf-8') as f:
    f.write('DSET ' + output_path + 'uvz.grd\n')
    f.write('TITLE 500hPa wind and height field\n')
    f.write('UNDEF -9999.0\n')
    f.write('XDEF ' + str(nx) + ' LINEAR ' + str(wlon) + ' ' + str(dx) + '\n')
    f.write('YDEF ' + str(ny) + ' LINEAR ' + str(slat) + ' ' + str(dy) + '\n')
    f.write('ZDEF 1 LINEAR 1 1\n')
    f.write('TDEF ' + str(int(all_time / output_interval) + 2) + ' LINEAR ' + str(start_time[11:13]) + 'Z' +
            start_time[8:10] + month[int(start_time[5:7])] + start_time[:4] + ' ' + str(int(output_interval / 3600)) + 'hr\n')
    f.write('VARS 3\n')
    f.write('U 0 99 500\n')
    f.write('V 0 99 500\n')
    f.write('Z 0 99 500\n')
    f.write('ENDVARS')
