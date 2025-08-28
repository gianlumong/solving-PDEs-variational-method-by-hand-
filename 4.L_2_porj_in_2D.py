import numpy as np
from matplotlib import pyplot as plt, cm


#definisco point matrix p e connectivity matrix
p=np.array([
    [0,1,2,2,0],
    [0,0,0,1,1]
])

T=np.array([
    [1,1,2],
    [4,2,3],
    [5,4,4]
])


#definisco la funzione che voglio approssimare 
f=lambda x: x[0]*x[1]

def polyarea(x, y):
    """Calcola l'area di un triangolo dati i suoi vertici (x, y)"""
    return 0.5 * abs((x[1] - x[0]) * (y[2] - y[0]) - 
                     (x[2] - x[0]) * (y[1] - y[0]))

def Global2D(p,T):
    n_p=p.shape[1] #numero di nodi N
    n_r=T.shape[1] #numero di elementi k

    M=np.zeros((n_p,n_p)) #inizializzo la global mass matrix 
    for k in range(n_r):

        loc2glob=T[:,k] -1  #sarebbero r,s,t? (-1 perché non c'è N_0, ma si parte direttamente da N_1 nel conteggio dei nodi)
        #(il -1 è da togliere se gli inidici dei nodi partissero da N_0!!!!)
        coordinate_nodi=p[:,loc2glob] 

        x=coordinate_nodi[0,:]
        y=coordinate_nodi[1,:]
        Area_k=polyarea(x,y)

        #calcolo la local element matrix 
        M_local=np.ones((3,3))
        for i in range(3):
            M_local[i,i]=2 #aggiorno solo gli elementi sulla diagonale
        M_local=M_local * Area_k * 1/12 

        for i in range(3):
            for j in range(3):
                M[loc2glob[i], loc2glob[j]] += M_local[i,j]

    return M, loc2glob

M, loc2 = Global2D(p,T)

#print('M= ', M)


#vado adesso a calcolare il load vector 
def load_b_vector(p,T,f):
    n_p=p.shape[1] #numero di nodi 
    n_r=T.shape[1] #numero di elementi (triangolari)

    b=np.zeros(n_p) #inizializzo b 
    for k in range(n_r):
        loc2glob=T[:,k] -1  #sarebbero r,s,t cioè gli indici dei nodi in P dell'elemento j-esimo 
        coordinate_nodi=p[:,loc2glob] 

        x=coordinate_nodi[0,:]
        y=coordinate_nodi[1,:]
        Area_k=polyarea(x,y)

        b_local=np.ones(3)
        for i in range(3):
            b_local[i]=1/3*f(([p[0,loc2glob[i]],p[1,loc2glob[i]]]))*Area_k
            b[loc2glob[i]]+=b_local[i]
    return b


b=load_b_vector(p,T,f)

print('shape(M): ', np.shape(M))

print('shape(b): ', np.shape(b))

Pf=np.linalg.solve(M,b)

print('shape(Pf): ', np.shape(Pf))




##===== PLOTTIAMO Pf ED f PER FARE UN CONFRONTO =====
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(14, 6))

# === SUBPLOT 1: funzione esatta ===
ax1 = fig.add_subplot(121, projection='3d')
#Ricordo che il dominio è un dominio rettangolare del tipo [0,2]x[0,1]
#decido di discretizzarlo lungo x e lungo y con 50 nodi
nx=50
ny=50
x_vals = np.linspace(0, 2, nx)
y_vals = np.linspace(0, 1, ny)


X, Y = np.meshgrid(x_vals, y_vals)
F = X * Y  # f(x, y)=x[0]*x[1]

ax1.plot_surface(X, Y, F, cmap='viridis', edgecolor='none', alpha=0.9)
ax1.set_title('Funzione esatta f(x, y) = x*y')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('f(x, y)')
ax1.view_init(elev=30, azim=135)

# === SUBPLOT 2: L²-proiezione ===
ax2 = fig.add_subplot(122, projection='3d')
x = p[0, :]
y = p[1, :]
triangles = (T - 1).T

ax2.plot_trisurf(x, y, Pf, triangles=triangles, cmap='viridis', edgecolor='k', linewidth=0.5, alpha=0.9)
ax2.set_title('Proiezione L² di f(x, y)')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('P_f(x, y)')
ax2.view_init(elev=30, azim=135)

plt.tight_layout()
plt.show()
