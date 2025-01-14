import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import pandas as pd
from scipy.integrate import odeint
import scipy.optimize as sco

lj = 0.55
b = lj / (2 * np.pi)
l_head = 3.41 #龙头长度
d_head = l_head - 0.275 * 2 #龙头凳两孔间距
l_body = 2.200
d_body = l_body - 0.275 * 2

theta = np.linspace(0,16 * 2 * np.pi,1000)
r = b * theta
x = r * np.cos(theta)
y = r * np.sin(theta)
plt.figure(1)
plt.plot(x,y)

def dtheta(theta,t):
    return -1.000 / (b * np.sqrt(1 + theta ** 2))

theta0 = 16 * 2 * np.pi
tspan = np.arange(0,301,0.5)
theta_head = odeint(dtheta,theta0,tspan)
x_head = b * theta_head * np.cos(theta_head) #龙头的x位置
y_head = b * theta_head * np.sin(theta_head) #龙头的y位置
for i in range(len(theta_head)):
    plt.plot(x_head[i],y_head[i],'r',marker='o')
plt.title("300s前龙头在螺线上的位置",fontproperties="SimHei")


#求龙身的位置函数
def f_body(x,y,d,theta):
    fun =lambda theta0:(b * theta0 * np.cos(theta0)- x) ** 2+ (b * theta0 * np.sin(theta0) - y) ** 2 - d ** 2
    dthe = 0.01
    theta0 = sco.fsolve(fun,theta + dthe)
    while theta0 <= theta or np.abs(b * theta0 - b * theta) > lj / 2:
        dthe = dthe + 0.01
        theta0 = sco.fsolve(fun,theta0 + dthe)
    return theta0


#求龙身尾的位置
n = 223
X = np.zeros((n + 1,len(x_head)))
Y = np.zeros((n + 1,len(x_head)))
Theta = np.zeros((n + 1 ,len(x_head)))
X[0,:] = np.transpose(x_head)
Y[0,:] = np.transpose(y_head)
Theta[0,:] = np.transpose(theta_head)
for i in range(len(x_head)):
    for j in range(n + 1):
        if j == 0:
            continue
        elif j == 1:
            d = d_head
        else:
            d = d_body
        Theta[j,i] = f_body(X[j-1,i],Y[j-1,i],d,Theta[j-1,i]) 
        X[j,i] = b * Theta[j,i] * np.cos(Theta[j,i])
        Y[j,i] = b * Theta[j,i] * np.sin(Theta[j,i])

#盘入动态图
fig2 = plt.figure(2)
ax = plt.subplot()
t = np.arange(0,301,0.5)
plt.plot(x,y)
point_ani, = plt.plot(X[:,0],Y[:,0],linestyle='-', marker='o', color='r', linewidth=1.2, markersize=6)
def updata(num):
    ax.set_title('t = %.1f' %(t[num]),loc='right',y=0)
    ax.set_title('300s前板凳龙行进动态图',loc='center',fontproperties="SimHei")
    point_ani.set_data(X[:,num],Y[:,num])
    return [point_ani]
anim = ani.FuncAnimation(fig2,updata,frames=601,interval=100,repeat=False)



#龙身速度
dt = 0.5
V = np.zeros((224,602))
V[:,0]= -b * np.sqrt(1.0 + np.power(Theta[:,0] , 2.0)) * (Theta[:,1]-Theta[:,0]) / dt
V[:,-1]= -b * np.sqrt(1.0 + np.power(Theta[:,-1] , 2.0)) * (Theta[:,-1]-Theta[:,-2]) / dt
V[:,1:-1]= -b * np.sqrt(1.0 + np.power(Theta[:,1:-1] , 2.0)) * (Theta[:,2:]-Theta[:,0:-2]) / 2.0 / dt

#验证龙头速度为1m/s
plt.figure(3)
plt.plot(tspan,V[0,:],linestyle='-',linewidth=1.2,color='r')
plt.ylim(0,1.1)   
plt.title('龙头速度检验图',fontproperties="SimHei")    

#导出数据
index = np.linspace(0,len(tspan),301)
Dataxy = np.zeros((2*(n + 1),len(index)))
Dataxy[0::2,:]=np.round(X[:,0::2],6)
Dataxy[1::2,:]=np.round(Y[:,0::2],6)
Datav = np.round(V[:,0::2],6)


##导出文件
# pd.DataFrame(X).to_excel('x.xlsx', index=False)
# pd.DataFrame(Y).to_excel('y.xlsx', index=False)
# pd.DataFrame(V).to_excel('v.xlsx', index=False)
pd.DataFrame(Dataxy).to_excel('dataxy.xlsx', index=False)
pd.DataFrame(Datav).to_excel('datav.xlsx', index=False)
pd.DataFrame(Theta).to_excel('theta.xlsx', index=False)

plt.show()





