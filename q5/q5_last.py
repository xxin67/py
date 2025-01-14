import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pandas as pd
from scipy.integrate import odeint
import scipy.optimize as sco
import scipy as sc
import matplotlib

#求龙身的位置函数
def f_body(x,y,d,theta):
    fun =lambda theta0:(b * theta0 * np.cos(theta0)- x) ** 2+ (b * theta0 * np.sin(theta0) - y) ** 2 - d ** 2
    dthe = 0.01
    theta0 = sco.fsolve(fun,theta + dthe)
    while theta0 <= theta or np.abs(b * theta0 - b * theta) > lj / 2:
        dthe = dthe + 0.1
        theta0 = sco.fsolve(fun,theta0 + dthe)
    return theta0

def solve_p21(x1,y1,theta1,d,r1,x_c,y_c,theta_max):
    del_theta = 2 * np.arcsin(d/2/r1)
    if del_theta <= theta_max - theta1:
        flag = 2
        theta = theta1 + del_theta
        x = x_c + r1 * np.cos(theta)
        y = y_c + r1 * np.sin(theta)
    else:
        theta = f_body(x1,y1,d,4.5/b)
        flag = 1
        x = b * theta * np.cos(theta)
        y = b * theta * np.sin(theta)
    return x,y,theta,flag

def solve_p32(x1,y1,theta1,d,r1,x_c1,y_c1,r2,x_c2,y_c2,theta_min):
    del_theta = 2 * np.arcsin(d/2/r2)
    if del_theta <= theta1 - theta_min:
        flag = 3
        theta = theta1 - del_theta
        x = x_c2 + r2 * np.cos(theta)
        y = y_c2 + r2 * np.sin(theta)
    else:
        di = np.sqrt((x1 - x_c1) ** 2 + (y1 - y_c1) ** 2)
        del_theta = np.arccos((di ** 2 + r1 ** 2 - d ** 2) / 2 /di /r1)
        theta_c1_di = np.arctan((y1-y_c1)/(x1-x_c1))
        theta = theta_c1_di + del_theta
        flag = 2
        x = x_c1 + r1 * np.cos(theta)
        y = y_c1 + r1 * np.sin(theta)
    return x,y,theta,flag 

def f_body2(x,y,d,theta):
    fun =lambda theta0:(b * (theta0 + np.pi) * np.cos(theta0)- x) ** 2+ (b * (theta0 + np.pi) * np.sin(theta0) - y) ** 2 - d ** 2
    dthe = -0.1
    theta0 = sco.fsolve(fun,theta + dthe)
    while theta0 > theta or np.abs(b * theta0 - b * theta) > lj / 2:
        dthe = dthe - 0.1
        theta0 = sco.fsolve(fun,theta0 + dthe)
    return theta0



def solve_p43(x1,y1,theta1,d,r2,x_c,y_c,theta_max):
    theta = f_body2(x1,y1,d,theta1)
    if theta >= 4.5 / b - np.pi:
        flag = 4
        x = b * (theta + np.pi) * np.cos(theta)
        y = b * (theta + np.pi) * np.sin(theta)
    else:
        func =lambda t:(x_c + r2 * np.cos(theta_max-t)- x1) ** 2+ (y_c + r2 * np.sin(theta_max - t ) - y1) ** 2 - d ** 2 
        q = -0.1
        del_theta = sco.fsolve(func,theta_max+q)
        theta = theta_max - del_theta
        flag = 3
        x = x_c + r2 * np.cos(theta)
        y = y_c + r2 * np.sin(theta)
    return x,y,theta,flag


lj = 1.7
b = lj / (2 * np.pi)
l_head = 3.41 #龙头长度
d_head = l_head - 0.275 * 2 #龙头凳两孔间距
l_body = 2.2
d_body = l_body - 0.275 * 2
V_0 = np.arange(1,2.1,0.1)
V_m = np.zeros((1,len(V_0)))

for vv in range(len(V_0)):
    v_0 = V_0[vv]
    R = 4.5
    theta_e = R / b
    theta_g = R / b - np.pi
    k = (b * np.sin(theta_e) + R * np.cos(theta_e)) / (b * np.cos(theta_e) - R * np.sin(theta_e))
    theta_m1 = np.arctan(-1 / k) + np.pi
    theta_st = np.arctan(np.tan(theta_e)) + np.pi - theta_m1
    r1_2 = R / np.cos (theta_st)
    r2 = r1_2 / 3
    r1 = r2 *2
    phi = 2 * theta_st
    s1 = r1 * (np.pi - phi)
    s2 = r2 * (np.pi - phi)
    theta_min1 = theta_m1 - s1 / r1
    theta_min2 = theta_min1 - np.pi
    theta_m2 = theta_min2 + s2 / r2
    x1 = R * np.cos(theta_e) + r1 * np.cos(theta_m1 - np.pi)
    y1 = R * np.sin(theta_e) + r1 * np.sin(theta_m1 - np.pi)

    x2 = R * np.cos(theta_g) - r2 * np.cos(theta_m2)
    y2 = R * np.sin(theta_g) - r2 * np.sin(theta_m2)


    def dtheta(theta,t):
        return v_0 / (b * np.sqrt(1 + theta ** 2))
    
    dt = 0.1
    Tto = (s1+s2) / v_0
    T = np.arange(0,Tto,dt)
    xx = np.zeros((224,len(T)))
    yy = np.zeros((224,len(T)))
    th = np.zeros((224,len(T)))
    tt_c1 = np.arange(0,s1/v_0,dt)
    theta_c1 = - v_0 * tt_c1 / r1  + theta_m1
    xx[0,0:len(tt_c1)] = r1 * np.cos(theta_c1) + x1
    yy[0,0:len(tt_c1)] = r1 * np.sin(theta_c1) + y1
    th[0,0:len(tt_c1)] = theta_c1
    tt_c2 = np.arange(tt_c1[-1]+dt,(s1+s2)/v_0,dt)
    theta_c2 = v_0 * (tt_c2 - s1/v_0) / r2 + theta_min1 -np.pi
    th[0,len(theta_c1):(len(theta_c1) + len(tt_c2))] = theta_c2
    xx[0,len(theta_c1):(len(theta_c1) + len(tt_c2))] = r2 * np.cos(theta_c2) + x2
    yy[0,len(theta_c1):(len(theta_c1) + len(tt_c2))] = r2 * np.sin(theta_c2) + y2

    N = 223
    t_total = np.arange(0,(s1+s2)/v_0,dt)
    for ii in range(len(t_total)):
        i = ii
        if t_total[ii] < 0:
            for j in range(1,N+1):
                d = d_head * (j<=1) + d_body * (j> 1)
                th[j,i] = f_body(xx[j-1,i],yy[j-1,i],d,th[j-1,i])
                xx[j,i] = b * th[j,i] * np.cos(th[j,i])
                yy[j,i] = b * th[j,i] * np.sin(th[j,i])
        
        elif t_total[ii] >=0 and t_total[ii] <=s1/v_0:
            flag = 2
            for j in range(1,N+1):
                d = d_head * (j <= 1) + d_body * (j > 1)
                if flag ==2 :
                    xi,yi,thetai,flag = solve_p21(xx[j-1,i],yy[j-1,i],th[j-1,i],d,r1,x1,y1,theta_m1)
                    th[j,i] = thetai
                    xx[j,i] = xi
                    yy[j,i] = yi
                else:
                    th[j,i] = f_body(xx[j-1,i],yy[j-1,i],d,th[j-1,i])
                    xx[j,i] = b * th[j,i] * np.cos(th[j,i])
                    yy[j,i] = b * th[j,i] * np.sin(th[j,i])

        elif t_total[ii] >s1/v_0 and t_total[ii] <=(s1+s2)/v_0:
            flag = 3
            for j in range(1,N+1):
                d = d_head * (j <= 1) + d_body * (j > 1)
                if flag == 3:
                    xi,yi,thetai,flag = solve_p32(xx[j-1,i],yy[j-1,i],th[j-1,i],d,r1,x1,y1,r2,x2,y2,theta_min2)
                    th[j,i] = thetai
                    xx[j,i] = xi
                    yy[j,i] = yi
                elif flag == 2:
                    xi,yi,thetai,flag = solve_p21(xx[j-1,i],yy[j-1,i],th[j-1,i],d,r1,x1,y1,theta_m1)
                    th[j,i] = thetai
                    xx[j,i] = xi
                    yy[j,i] = yi    
                else:
                    th[j,i] = f_body(xx[j-1,i],yy[j-1,i],d,th[j-1,i])
                    xx[j,i] = b * th[j,i] * np.cos(th[j,i])
                    yy[j,i] = b * th[j,i] * np.sin(th[j,i])
        else:
            flag = 4
            for j in range(1,N+1):
                d = d_head * (j <= 1) + d_body * (j > 1)
                if flag == 4:    
                    xi,yi,thetai,flag = solve_p43(xx[j-1,i],yy[j-1,i],th[j-1,i],d,r2,x2,y2,theta_m2)
                    th[j,i] = thetai
                    xx[j,i] = xi
                    yy[j,i] = yi   
                
                elif flag == 3:
                    xi,yi,thetai,flag = solve_p32(xx[j-1,i],yy[j-1,i],th[j-1,i],d,r1,x1,y1,r2,x2,y2,theta_min2)
                    th[j,i] = thetai
                    xx[j,i] = xi
                    yy[j,i] = yi
                elif flag == 2:
                    xi,yi,thetai,flag = solve_p21(xx[j-1,i],yy[j-1,i],th[j-1,i],d,r1,x1,y1,theta_m1)
                    th[j,i] = thetai
                    xx[j,i] = xi
                    yy[j,i] = yi    
                else:
                    th[j,i] = f_body(xx[j-1,i],yy[j-1,i],d,th[j-1,i])
                    xx[j,i] = b * th[j,i] * np.cos(th[j,i])
                    yy[j,i] = b * th[j,i] * np.sin(th[j,i])
    

    Vx = np.zeros((len(xx),len(t_total)))
    Vy = np.zeros((len(yy),len(t_total)))
    Vx[:,0] = (xx[:,1] - xx[:,0]) / dt
    Vx[:,-1] = (xx[:,-1]- xx[:,-2] ) / dt 
    Vx[:,1:-1] = (xx[:,2:] - xx[:,0:-2]) / 2/ dt

    Vy[:,0] = (yy[:,1] - yy[:,0]) / dt
    Vy[:,-1] = (yy[:,-1]- yy[:,-2] ) / dt 
    Vy[:,1:-1] = (yy[:,2:] - yy[:,0:-2]) / 2/ dt

    V = np.sqrt(Vx ** 2 + Vy ** 2)
    V_m [0,vv] = np.max(np.max(V))

print(len(V_0),len(V_m))
p = np.polyfit(V_0,np.transpose(V_m),1)
p1 = plt.plot(V_0,np.transpose(V_m),color = 'r',linewidth = 1)
plt.xlabel('龙头速度',fontproperties="SimHei")
plt.ylabel('队伍把手最高速度',fontproperties="SimHei")
p2 = plt.plot(np.linspace(V_0[0],V_0[-1],100),np.polyval(p,np.linspace(V_0[0],V_0[-1],100)),color =np.random.random((1,3)) )
plt.plot([V_0[0],V_0[-1]],[2,2],color = 'k',linestyle = '--')
plt.legend('队伍最大速度与龙头速度数据点','线性拟合关系')
v0max = (2 - p[1]) / p[0]

print("龙头最大行进速度：", v0max)

plt.show()

