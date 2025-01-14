import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pandas as pd
from scipy.integrate import odeint
import scipy.optimize as sco
import matplotlib


lj = 0.55
b = lj / (2 * np.pi)
l_head = 3.41 #龙头长度
d_head = l_head - 0.275 * 2 #龙头凳两孔间距
l_body = 2.2
d_body = l_body - 0.275 * 2
dt = 0.2
flag = 0    #是否碰撞判断指标
step = 1
t = 0


theta_0 = np.linspace(0,16 * 2 * np.pi,1000)
r = b * theta_0
x = r * np.cos(theta_0)
y = r * np.sin(theta_0)

plt.ion()

def dtheta(theta,t):
    return -1 / (b * np.sqrt(1 + theta ** 2))

theta_0 = 57.032076651015522

#求龙身的位置函数
def f_body(x,y,d,theta):
    fun =lambda theta0:(b * theta0 * np.cos(theta0)- x) ** 2+ (b * theta0 * np.sin(theta0) - y) ** 2 - d ** 2
    dthe = 0.01
    theta0 = sco.fsolve(fun,theta + dthe)
    while theta0 <= theta or np.abs(b * theta0 - b * theta) > lj / 2:
        dthe = dthe + 0.1
        theta0 = sco.fsolve(fun,theta0 + dthe)
    return theta0


#求龙身尾的位置
N = 223
X = np.nan * np.zeros((N + 1,3))
Y = np.nan * np.zeros((N + 1,3))
Theta = np.nan * np.zeros((N + 1 ,3))
Theta[0,2] = theta_0


def f_judge(l1,x1_1,x1_2,l2,x2_1,x2_2,n,m):
    k1 = (x1_1[1] - x1_2[1]) / (x1_1[0] -x1_2[0])
    k1_ = -1 / k1
    k2 = (x2_1[1] - x2_2[1]) / (x2_1[0] - x2_2[0])
    k2_ = -1 / k2
    x1_cen = (x1_1 + x1_2) / 2
    x2_cen = (x2_1 + x2_2) / 2
    A = np.array([[k1_[0],-1],[k2_[0],-1]])
    P = np.linalg.lstsq(A,[k1_[0] * x1_cen[0] - x1_cen[1],k2_[0]* x2_cen[0] - x2_cen[1]])
    vec1 = x1_cen - P[0]
    vec2 = x2_cen - P[0]
    theta1 = np.angle(complex(vec1[0],vec1[1]))
    theta2 = np.angle(complex(vec2[0],vec2[1]))
    d_theta = theta1 - theta2
    d_1 = np.linalg.norm(vec1)
    d_2 = np.linalg.norm(vec2)
    x_0,y_0 = np.meshgrid(np.linspace(d_1 - 0.15,d_1 + 0.15,n),np.linspace(-l1 / 2,l1 / 2,m))
    T = np.array([[np.cos(d_theta),-np.sin(d_theta)],[np.sin(d_theta),np.cos(d_theta)]])
    xy_new = T @ np.array([x_0.ravel(order='f'),y_0.ravel(order='f')])
    a1 = (np.abs(xy_new[0,:] - d_2) < 0.15) 
    a2 = (np.abs(xy_new[1,:]-0) < l2/2)
    for jj in range(len(a1)):
        if a1[jj] == True and a2[jj] == True:
            flag = True
        else:
            flag = False
    return flag,xy_new

while flag == 0:

    t = 300 + step * dt
    X[:,0] = X[:,2]
    Y[:,0] = Y[:,2]
    Theta[:,0] = Theta[:,2]
    tspan = np.linspace(0,dt,3)
    theta = odeint(dtheta,Theta[0,0],tspan)
    X1 = b * theta * np.cos(theta)
    Y1 = b * theta * np.sin(theta)
    X[0,:] = np.transpose(X1)
    Y[0,:] = np.transpose(Y1)
    Theta[0,:] = np.transpose(theta)
    for i in range(1,len(tspan)):
        for j in range(1,N+1):
            if j == 1:
                d = d_head
            else:
                d = d_body
            Theta[j,i] = f_body(X[j-1,i],Y[j-1,i],d,Theta[j-1,i])
            X[j,i] = b * Theta[j,i] * np.cos(Theta[j,i])
            Y[j,i] = b * Theta[j,i] * np.sin(Theta[j,i])

    fig1 = plt.figure(1)
    plt.clf()
    plt.plot(x,y)
    plt.plot(X[:,-1],Y[:,-1],linestyle='-', marker='o', color='r', linewidth=1.2, markersize=6)
    plt.title('t = %.1f' % t,loc = 'right',y = 0)
    plt.title('从300s开始板凳龙行进碰撞模拟图', fontproperties="SimHei")
    plt.pause(0.001)

    #判断碰撞  
    for k in range(N):
        x_1 = X[k,1]
        x_2 = X[k+1,1]
        y_1 = Y[k,1]
        y_2 = Y[k+1,1]
        theta_1 = Theta[k,1]
        theta_2 = Theta[k+1,1]
        index1 = np.transpose(np.where((theta_1 + 2 * np.pi - Theta[:,1]) > 0))
        index1 = index1[-3:]
        index2 = np.transpose(np.where(Theta[:,1] - (theta_2 + 2 * np.pi) > 0))
        if len(index2) == 0:
            break
        else:
            index2 = index2[0:min(3,len(index2))]
        index_i =np.arange(index1[0], index2[-1] + 1,1)
        n = 10
        m = 20
        for ii in range(len(index_i)-1):
            x2_1 =np.array( [[X[index_i[ii],1]],[Y[index_i[ii],1]]])
            x2_2 = np.array([[X[index_i[ii+1],1]],[Y[index_i[ii+1],1]]])
            judge,xy = f_judge(l_head * (k<=0) + l_body * (k>0),np.array([[x_1],[y_1]]),np.array([[x_2],[y_2]]),l_body,x2_1,x2_2,n,m)
            if judge != False:
                flag = 1
                break
        if flag == 1:
            break
    step = step + 1

plt.ioff()
#求速度
V = -b * np.sqrt(1 + Theta[:,-1] ** 2) * (Theta[:,-1]-Theta[:,-2]) / (dt/2)
plt.figure(2)
plt.plot(V,linestyle='-',linewidth=1.2,color='r') 
plt.ylim(0,1.1)
plt.title("碰撞停止时各个板凳的速度图",fontproperties="SimHei" )

#最终时刻图
plt.figure(3)
plt.plot(x,y)
plt.plot(X[:,-1],Y[:,-1],linestyle='-', marker='o', color='r', linewidth=1.2, markersize=6)
plt.title("碰撞停止时各个板凳位置图",fontproperties="SimHei" )

#
#  #导出数据
Dataxyv = np.zeros((N+1 ,3))
Dataxyv[:,0]=np.round(X[:,-1],6)
Dataxyv[:,1]=np.round(Y[:,-1],6)
Dataxyv[:,2]=np.round(V[:],6)

pd.DataFrame(Dataxyv).to_excel('dataxyv.xlsx', index=False)
plt.show()




