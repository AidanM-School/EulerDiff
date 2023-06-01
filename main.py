import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
def f(x,alpha,beta,delta,gamma):
    return np.array((alpha*x[0]-beta*x[0]*x[1],delta*x[0]*x[1]-gamma*x[1]))

def g(x):
    return np.array((x[1],np.sin(x[0])))

def forward_em(n,dt):
    xs=np.zeros(n)
    alpha=-0.5
    xs[0]=1
    ts=np.linspace(0,n*dt,n)
    for i in range(1,n):
        xs[i]=dt*f(xs[i-1],alpha)+xs[i-1]
    return ts,xs

def backward_em(n,dt):
    xs = np.zeros(n)
    alpha = -0.5
    xs[0] = 1
    ts = np.linspace(0, n * dt, n)
    for i in range(1,n):
        xs[i]=root(lambda x:xs[i-1]+dt*f(x,alpha)-x,xs[i-1])['x']
    return ts,xs

def twod_fem(n,dt):
    thetas=np.zeros(n)
    omegas=np.zeros(n)
    ts = np.linspace(0,n*dt,n)
    thetas[0]=.3
    omegas[0]=0
    for i in range(1,n):
        thetas[i]=dt*omegas[i-1]+thetas[i-1]
        omegas[i]=dt*np.sin(thetas[i-1])+omegas[i-1]
    return ts, thetas, omegas

def twod_bem(n,dt):
    thetas = np.zeros(n)
    omegas = np.zeros(n)
    ts = np.linspace(0, n * dt, n)
    thetas[0] = .3
    omegas[0] = 0
    for i in range(1, n):
        x1 = np.array((thetas[i - 1], omegas[i - 1]))
        x2 = root(lambda x: list(x1 + dt * g(x) - x), x1)['x']
        thetas[i], omegas[i] = list(x2)[0], list(x2)[1]
    return ts, thetas, omegas

def pred_prey(n,dt):
    preys=np.zeros(n)
    preds=np.zeros(n)
    alpha=1.1
    beta=0.4
    delta=.1
    gamma=.4
    ts = np.linspace(0,n*dt,n)
    preys[0]=10
    preds[0]=10
    for i in range(1,n):
        x1=np.array((preys[i-1],preds[i-1]))
        x2=root(lambda x:list(x1+dt*f(x,alpha,beta,delta,gamma)-x),x1)['x']
        preys[i],preds[i]=list(x2)[0],list(x2)[1]
    return ts, preys, preds

def main():
    n=100000
    dt=.0003
    xs,preys,preds=pred_prey(n,dt)
    """xs, ys = forward_em(n, dt)
    plt.plot(xs,ys)
    xs1,ys1=backward_em(n,dt)
    plt.plot(xs1,ys1)"""
    plt.plot(xs,preys)
    plt.plot(xs,preds)
    plt.show()


main()