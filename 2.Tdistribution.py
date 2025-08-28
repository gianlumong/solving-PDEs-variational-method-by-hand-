import numpy as np
import matplotlib.pyplot as plt

#sto considerando un problema per trovare la distribuziome della temperatura
#lungo un'asta lunga L

def GlobalStiffnesA(x,a,k):
    n=len(x)-1
    A=np.zeros((n+1,n+1))
    for i in range(n):
        h=x[i+1]-x[i]
        x_med=(x[i+1]+x[i])/2
        ai=a(x_med)
        A[i,i]=A[i,i]+ai/h
        A[i,i+1]=A[i,i+1]-ai/h
        A[i+1,i]=A[i+1,i]-ai/h
#        A[i+1,i+1]=A[i+1,i+1]+ai/h
        A[i+1,i+1]=ai/h
    A[0,0]+=k[0]
    A[n,n]+=k[1]
    return A


def loadvector(x,f,k,g):
    n=len(x)-1
    b=np.zeros(n+1)
    for i in range(n):
        h=x[i+1]-x[i]
        b[i]=b[i]+f(x[i])*h/2
#        b[i+1]=b[i+1]+f(x[i])*h/2
        b[i+1]=f(x[i+1])*h/2
    b[0]+=k[0]*g[0]
    b[n]+=k[1]*g[1]
    return b 

h=0.1
n_step=int(6/h)
x=np.linspace(2,8,n_step)
k=np.array([1e+6,0])
g=np.array([-1,0])
a = lambda x: 0.1*(5 - 0.6 * x)
f= lambda x: 0.03*(x-6)**4

A=GlobalStiffnesA(x,a,k)
b=loadvector(x,f,k,g)
u=np.linalg.solve(A,b)

plt.plot(x,u,label='u')
plt.legend()
plt.xlabel('x')
plt.ylabel('T')
plt.show()

