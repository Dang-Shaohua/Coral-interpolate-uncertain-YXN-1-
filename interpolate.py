# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:20:09 2023

@author: dsh19
"""

import pycwt
import math
import numpy as np 
import xlrd
import sys
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import pyleoclim as pyleo
plt.style.use('ggplot')
sys.path.append('E:/YX10_60_2/sr/重新建立温度计校准/bandpass.py') # 使用方式：bandpass.butter_bandpass_filter(data,lowcut,highcut,fs),以月分辨率ENSO为例 1/84，1/36 fs一般选择1
import bandpass
sys.path.append('E:/YX10_60_2/sr/重新建立温度计校准/lowpass.py')
import lowpass

#读取EXCEL数据并存入numpy数组 该文件重新定义校准了西沙YXN-1srca温度计，相比过去提升了 R2（0.72-0.77）
DATA = xlrd.open_workbook('E:/YX10_60_2/sr/重新建立温度计校准/modern coral srca.xlsx')
sheet1 = DATA.sheets()[0]
'a = [i for i in a if i != '']' #列表剔除含空值的字符串
depth = np.array([i for i in sheet1.col_values(0)[1:] if i !=''],dtype = np.float64)
srca_original = np.array([i for i in sheet1.col_values(1)[1:] if i !=''],dtype = np.float64)
ctl_points = np.array([i for i in sheet1.col_values(5)[1:] if i !=''],dtype = np.float64)
depth_ctlpoints = np.array([i for i in sheet1.col_values(4)[1:] if i !=''],dtype = np.float64)



#定义珊瑚一次性插值函数
def coral_chronology_model(depth,depth_ctlpoints,ctl_points_age,data,interpolate_rate): #depth为总深度,其实就是数据涵盖的长度;depth_ctlpoints为控制点所在的深度;ctl_points_age控制点的年龄; data为数据;interpolate_rate为想插值的间隔或目标分辨率
    f_depth_ctlpoints_age = interpolate.interp1d(depth_ctlpoints, ctl_points_age)
    depth_age = f_depth_ctlpoints_age(depth) #根据控制点得出原始数据每个深度的年龄
    depth_age_even = [] #制作空矩阵存放插值后的时间序列
    npts = len(depth_age)
    depth_age_even = np.arange(depth_age[0],depth_age[npts-1]+0.01,interpolate_rate) # np.arange()步长为1/12 可能需要 17-18以上的长度才能保证数据的长度 可以考depth_age[npts-1] +0.01
    age = depth_age_even #插值结果对应的年龄
    dataeven = interpolate.interp1d(depth_age,data)
    interpolate_result = dataeven(age) #插值结果
    return age,interpolate_result

age,srca_interpolated = coral_chronology_model(depth,depth_ctlpoints,ctl_points,srca_original,1/12)
age = age[1:]
srca_interpolated = srca_interpolated[1:]

#定义珊瑚srca误差分析函数
def coral_uncertain_montel(data,montel_number,analytical_mean,analytical_sigma,colony_mean,colony_sigma,
                           slope_mean,slope_sigma,intercept_mean,intercept_sigma):

   testdata = data.reshape(len(data),1)
   #analytical_mean = 0
   #analytical_sigma = 0.016
   analytical_error = []
   for i in range(montel_number):
       analytical_error_pseudo = np.random.normal(analytical_mean, analytical_sigma,(len(data),))
       analytical_error.append(analytical_error_pseudo)
   #print(analytical_error)
   analytical_error = np.array(analytical_error,dtype = float).T
   error1 = analytical_error + testdata #测试误差叠加原始数据形成 error1

   #colony_mean = 0
   #colony_sigma = 0.024
   colony_error = []
   for i in range(montel_number):
       colony_error_pseudo = np.random.normal(colony_mean, colony_sigma,(len(data),))
       colony_error.append(colony_error_pseudo)
   #print(colony_error)
   colony_error = np.array(colony_error,dtype = float).T
   error2 = error1 + colony_error #测试误差叠加种间误差 形成error2
   #slope_mean = -0.05233 
   #slope_sigma = 0.001363
   slope_error = []
   for i in range(montel_number):
       slope_error_pseudo = np.random.normal(slope_mean, slope_sigma,(len(data),))
       slope_error.append(slope_error_pseudo)
   #print(slope_error)
   slope_error = np.array(slope_error,dtype = float).T

   #intercept_mean = 10.214
   #intercept_sigma = 0.03732
   intercept_error = []
   for i in range(montel_number):
       intercept_error_pseudo = np.random.normal(intercept_mean,intercept_sigma,(len(data),))
       intercept_error.append(intercept_error_pseudo)
   #print(intercept_error)
   intercept_error = np.array(intercept_error,dtype = float).T

   formular_error = (error2 - intercept_error)/slope_error
   sstdata = (testdata - intercept_mean)/slope_mean
   rmse = np.sqrt((((formular_error - sstdata)**2).sum(axis = 1))/montel_number)
   rmse = rmse.reshape(len(rmse),1)

   sst_up = sstdata + rmse
   sst_down = sstdata - rmse
   sst_uncertain = np.concatenate((sst_up,sstdata,sst_down),axis = 1)
   return rmse,sstdata,sst_uncertain
#输出结果包括 每个点RMSE，srca-sst，sst在每个点的不确定性最值
rmse,sstdata,sst_uncertain = coral_uncertain_montel(srca_interpolated,2000,0,0.016,0,0.024,-0.05233,0.001363,10.214,0.03732)



