import numpy as np
import matplotlib.pyplot as plt


#codice per effettuare il confronto fra f e Pf (tramite L2-projection method)
# con f= x*sin(x)

#definisco la funzione f
f=lambda x: 2*x*np.sin(2*np.pi*x) + 3 

#definisco l'intervallo I in cui la voglio approssimare 
n_nodi=6
n_interval=n_nodi-1
I=[0,1]  #dominio
h=I[-1]/n_interval #intervalli mesh
x=np.linspace(I[0],I[-1],n_nodi)

#voglio andarla ad approssimare tramite la L2-norm Phf
#calcolo la global matrix
def M_global(x):
    n_nodi=int(len(x))
    M=np.zeros((n_nodi,n_nodi))
    #calcolo le local element matrix Mi con i=1,..n
    for i in range(1,n_nodi):
        h=x[i]-x[i-1]
        M_i=np.zeros((n_nodi,n_nodi))
        M_i[i-1,i-1]=h/3
        M_i[i,i-1]=h/6
        M_i[i-1,i]=h/6
        M_i[i,i]=h/3
        M+=M_i    
    return M

def b_load(x,f):
    n_nodi=int(len(x))
    b=np.zeros(n_nodi)
    #calcolo i local load vector bi con i=1,..n
    for i in range(1,n_nodi):
        h=x[i]-x[i-1]
        b_i=np.zeros(n_nodi)
        b_i[i-1]=f(x[i-1])*h/2
        b_i[i]=f(x[i])*h/2
        b+=b_i    
    return b

M=M_global(x)
b=b_load(x,f)
Pf=np.linalg.solve(M,b)

print('shape(M): ', np.shape(M))

print('shape(b): ', np.shape(b))

print('shape(Pf): ', np.shape(Pf))


plt.plot(x,Pf,label='Pf')
x=np.linspace(I[0],I[-1],500)
plt.plot(x,f(x),label='f')
plt.title('costruzione con somma sugli elementi')
plt.legend()
plt.xlabel('x')
plt.ylabel('f')
plt.show()
