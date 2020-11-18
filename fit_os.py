#!/usr/bin/python
# coding=utf-8
import sys
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
def read_data(filename):
    """ 读取能量和体积，文件的形式如下，第一列为晶格常数，第二列为体积，第三列为能量
        3.74 136.2000 -105.833996
        3.75 136.9300 -105.865334
        3.76 137.6600 -105.892136
        3.78 139.1300 -105.928546
        3.79 139.8600 -105.944722
        3.80 140.6000 -105.955402
        3.81 141.3400 -105.960574
        3.82 142.0900 -105.960563
        3.83 142.8300 -105.954437
        3.84 143.5800 -105.949877
    """
    data = np.loadtxt(filename)     
    return data[:,1], data[:,2] 
def eos_murnaghan(vol, E0, B0, BP, V0):
    """ 输入体积和方程参数，返回能量
    Birch-Murnaghan方程：Phys. Rev. B 1983, 28 (10), 5480–5486.
    """
    return E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1)
def fit_murnaghan(volume, energy):
    """ 拟合Murnaghan状态方程，返回最优参数 
    """
    # 用一元二次方程拟合，得到初步猜测参数
    p_coefs = np.polyfit(volume, energy, 2)
    # 抛物线的最低点 dE/dV = 0 ( p_coefs = [c,b,a] ) V(min) = -b/2a
    p_min = - p_coefs[1]/(2.*p_coefs[0])
    # warn if min volume not in result range 
    if (p_min < volume.min() or p_min > volume.max()):
        print("Warning: minimum volume not in range of results")
    # 从抛物线最低点估计基态能量    
    E0 = np.polyval(p_coefs, p_min)
    # 体积模量估计
    B0 = 2.*p_coefs[2]*p_min
    # 初步猜测参数 (BP通常很小，取4)
    init_par = [E0, B0, 4, p_min]
    print("guess parameters:")
    print(" V0     =  {:1.4f} A^3 ".format(init_par[3]))
    print(" E0     =  {:1.4f} eV  ".format(init_par[0]))
    print(" B(V0)  =  {:1.4f} eV/A^3".format(init_par[1]))
    print(" B'(VO) =  {:1.4f} ".format(init_par[2]))
    best_par, cov_matrix = curve_fit(eos_murnaghan, volume, energy, p0 = init_par)
    return best_par
def fit_and_plot(filename):
    """ 读取文件数据，拟合Murnaghan状态方程，返回最优参数和E-V图形
    """
    # 从文件读取数据
    volume, energy = read_data(filename)   
    # 用Murnaghan状态方程拟合数据
    best_par = fit_murnaghan(volume, energy)
    # 输出最优参数   
    print("Fit parameters:")
    print(" V0     =  {:1.4f} A^3 ".format(best_par[3]))
    print(" E0     =  {:1.4f} eV  ".format(best_par[0]))
    print(" B(V0)  =  {:1.4f} eV/A^3".format(best_par[1]))
    print(" B'(VO) =  {:1.4f} ".format(best_par[2]))
    # 用拟合的参数生成Murnaghan模型
    m_volume = np.linspace(volume.min(), volume.max(), 1000) 
    m_energy = eos_murnaghan(m_volume, *best_par) 
    # 画E-V图
    lines = plt.plot(volume, energy, 'ok', m_volume, m_energy, '--r' )
    plt.xlabel(r"Volume [$\rm{A}^3$]")
    plt.ylabel(r"Energy [$\rm{eV}$]")
    return  best_par, lines
fit_and_plot("SUMMARY.dat")
plt.draw()
plt.show()
