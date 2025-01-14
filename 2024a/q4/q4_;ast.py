import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pandas as pd
from scipy.integrate import odeint
import scipy.optimize as sco
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

def f_body2(x,y,d,theta1):
    func =lambda theta0:(b * (theta0 + np.pi) * np.cos(theta0)- x) ** 2+ (b * (theta0 + np.pi) * np.sin(theta0) - y) ** 2 - d ** 2
    dthe = -0.1
    theta0 = sco.fsolve(func,theta1 + dthe)
    while theta0 >= theta1 or np.abs(b * theta0 - b * theta1) > lj / 2:
        dthe = dthe - 0.1
        theta0 = sco.fsolve(func,theta0 + dthe)
    return theta0



def solve_p43(x1,y1,theta1,d,r2,x_c,y_c,theta_max):
    theta = f_body2(x1,y1,d,theta1)
    if theta >= 4.5 / b - np.pi:
        flag = 4
        x = b * (theta + np.pi) * np.cos(theta)
        y = b * (theta + np.pi) * np.sin(theta)
    else:
        funct =lambda t:(x_c + r2 * np.cos(theta_max-t)- x1) ** 2+ (y_c + r2 * np.sin(theta_max - t ) - y1) ** 2 - d ** 2 
        q = -0.1
        del_theta = sco.fsolve(funct,theta_max+q)
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
v_0 = 0.1

#盘入盘出掉头图像
plt.figure(1)
theta_0 = np.arange(10 * np.pi,0,-0.01)
r = b * theta_0
x = r * np.cos(theta_0)
y = r * np.sin(theta_0)
plt.plot(x,y,color = np.sort(np.random.random((1,3))))

theta_0 = theta_0 - np.pi
r_2 = b * (theta_0 + np.pi)
x2 = r_2 * np.cos(theta_0)
y2 = r_2 * np.sin(theta_0)
x_diao = x2
y_diao = y2
plt.plot(x2,y2,color = np.sort(np.random.random((1,3))),linewidth = 1.3)

R = 4.5
x_d = R * np.cos(theta_0)
y_d = R * np.sin(theta_0)
plt.plot(x_d,y_d,color = 'k',linewidth = 2)

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

plt.plot(x1 + r1 * np.cos(np.linspace(theta_min1,theta_m1,50)),y1 + r1 * np.sin(np.linspace(theta_min1,theta_m1,50)),'r',linewidth = 2)
plt.plot(x1,y1,'r',marker='o')
plt.plot(x2 + r2 * np.cos(np.linspace(theta_min2,theta_m2,50)),y2 + r2 * np.sin(np.linspace(theta_min2,theta_m2,50)),'b',linewidth = 2)
plt.plot(x2,y2,'b',marker='o')
plt.title("调头轨迹图",fontproperties="SimHei")


def dtheta(theta,t):
    return 1 / (b * np.sqrt(1 + theta ** 2))

theta0 = theta_e
dt = 0.1
tspan_s = 0
tspan_e = 100
tspan = np.arange(tspan_s,tspan_e+dt,dt)
theta = odeint(dtheta,theta0,tspan)
X1 = b * theta * np.cos(theta)
Y1 = b * theta * np.sin(theta)
tt = np.transpose(tspan)
tt_e = tt[-1::-1]
xx = np.zeros((224,int(200 / dt + 1)))
yy = np.zeros((224,int(200 / dt + 1)))
th = np.zeros((224,int(200 / dt + 1)))
xx[0,0:len(X1)] = X1[-1::-1,0]
yy[0,0:len(Y1)] = Y1[-1::-1,0]
th[0,0:len(theta)] = theta[-1::-1,0]

tt_1 = np.arange(dt,s1,dt)
theta_1 = - tt_1 / r1 + theta_m1
th[0,len(theta):(len(theta) + len(tt_1))] = theta_1
xx[0,len(X1):(len(X1) + len(tt_1))] = r1 * np.cos(theta_1) + x1
yy[0,len(Y1):(len(Y1) + len(tt_1))] = r1 * np.sin(theta_1) + y1

tt_2 = np.arange(tt_1[-1]+dt,s1+s2,dt)
theta_2 = (tt_2 - s1) / r2 + theta_min1 - np.pi
th[0,(len(theta) + len(theta_1)):(len(theta) + len(theta_1) +len(tt_2))] = theta_2
xx[0,(len(X1) + len(theta_1) ):(len(X1) + len(theta_1) + len(tt_2))] = r2 * np.cos(theta_2) + x2
yy[0,(len(Y1) + len(theta_1) ):(len(Y1) + len(theta_1) + len(tt_2))] = r2 * np.sin(theta_2) + y2

def dtheta_1(theta,t):
    return 1 / (b * np.sqrt(1 + (theta + np.pi) ** 2))

theta0 = theta_g

tspan_s = (tt_2[-1] + dt)
tspan_e = 100
tspan = np.arange(tspan_s,tspan_e,dt)


theta2 = odeint(dtheta_1,theta0,tspan)
X2 = b * (theta2 + np.pi) * np.cos(theta2)
Y2 = b * (theta2 + np.pi) * np.sin(theta2)
th[0,(len(theta)+len(theta_1)+len(tt_2)):] = theta2[0:,0]
xx[0,(len(theta)+len(theta_1)+len(tt_2)):] = X2[0:,0]
yy[0,(len(theta)+len(theta_1)+len(tt_2)):] = Y2[0:,0]

for i in range(0,len(th[0,:]),3):
    plt.plot(xx[0,i],yy[0,i],marker = 'o',markersize = 3,color = 'r')
    plt.pause(0.001)

N = 223
t_all = np.arange(-100,100+dt,dt)
for ii in range(len(t_all)):
    i =  ii
    if t_all[ii] <=0:
        for j in range(1,N+1):
            d = d_head * (j<=1) + d_body * (j> 1)
            th[j,i] = f_body(xx[j-1,i],yy[j-1,i],d,th[j-1,i])
            xx[j,i] = b * th[j,i] * np.cos(th[j,i])
            yy[j,i] = b * th[j,i] * np.sin(th[j,i])
    
    elif t_all[ii] >0 and t_all[ii] <=s1:
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

    elif t_all[ii] >s1 and t_all[ii] <=(s1+s2):
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


plt.figure(4)
for j in range(0,(len(t_all) + int(np.round((t_all[0] + 100) / dt))),2) :
    plt.clf()
    plt.plot(x,y,color = 'g')
    plt.plot(x_diao,y_diao,color = 'k',linewidth = 1.3)
    plt.plot(x_d,y_d,color = 'm',linewidth = 2)
    plt.plot(x1 + r1 * np.cos(np.linspace(theta_min1,theta_m1,50)),y1 + r1 * np.sin(np.linspace(theta_min1,theta_m1,50)),'r',linewidth = 2)
    plt.plot(x1,y1,'r',marker='o')
    plt.plot(x2 + r2 * np.cos(np.linspace(theta_min2,theta_m2,50)),y2 + r2 * np.sin(np.linspace(theta_min2,theta_m2,50)),'b',linewidth = 2)
    plt.plot(x2,y2,'b',marker='o')
    plt.plot(xx[:,j],yy[:,j],linewidth = 1.2,marker = 'o',markersize = 6 ,linestyle='-')
    plt.title("调头模拟图",fontproperties="SimHei")
    plt.title("t = %.1f" % (t_all[j]),loc='right')
    plt.pause(0.01)

Vx = np.zeros((len(xx),len(t_all)))
Vy = np.zeros((len(yy),len(t_all)))
Vx[:,0] = (xx[:,1] - xx[:,0]) / dt
Vx[:,-1] = (xx[:,-1]- xx[:,-2] ) / dt 
Vx[:,1:-1] = (xx[:,2:] - xx[:,0:-2]) / 2/ dt

Vy[:,0] = (yy[:,1] - yy[:,0]) / dt
Vy[:,-1] = (yy[:,-1]- yy[:,-2] ) / dt 
Vy[:,1:-1] = (yy[:,2:] - yy[:,0:-2]) / 2/ dt

V = np.sqrt(Vx ** 2 + Vy ** 2)

nn = 1 / dt 
index =  np.transpose(np.arange(0,len(t_all),nn))
dataxy = np.zeros((2*(N+1),len(index)))
dataxy[0::2,:] = np.round(xx[0::,0::10],6)
dataxy[1::2,:] = np.round(yy[0::,0::10],6)
datav = np.round(V[0::,0::10],6)  

print("最短调头曲线长度 = %.6f"  %(s1 + s2))
pd.DataFrame(dataxy).to_excel('dataxy4.xlsx', index=False)
pd.DataFrame(datav).to_excel('datav4.xlsx', index=False)
plt.pause(0)