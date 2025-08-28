import numpy as np
import matplotlib.pyplot as plt


#andiamo a risolvere un problema di Poisson con Dirichelet BC
# -u''= delta (con delta= e^(-c*|x-0.5|^2)
# e con le seguenti Dirichelet BC: u(0)=u(-1)=0
# Scelgo di risolvere il problema prima su di una mesh rada (di 5 nodi)
# e poi di andare a fare una mesh refinement (... volte)!!
n_loop_refinement=10

# Parametri iniziali
dominio = np.array([0, 1])
n_nodi_iniz = 5
x = np.linspace(dominio[0], dominio[1], n_nodi_iniz)
n_el = len(x) - 1  # Numero di elementi

# Source term (delta approssimata)
c = 100
delta = lambda x: np.exp(-c * np.abs(x - 0.5)**2)
f = delta

#dato che ho delle Dirichelet costruisco A e b con solo gli elementi interni!!!!!
def A_global(x):
    n_nodi=int(len(x))
    M=np.zeros((n_nodi,n_nodi))
    #calcolo le local element matrix Mi con i=1,..n
    for i in range(1,n_nodi):
        h=x[i]-x[i-1]
        M_i=np.zeros((n_nodi,n_nodi))
        M_i[i-1,i-1]=1/h
        M_i[i,i-1]=-1/h
        M_i[i-1,i]=-1/h
        M_i[i,i]=1/h
        M+=M_i    
    M=M[1:-1,1:-1] #considero solo gli elementi interni!!
    return M

def b_load(x,f):
    n_nodi=int(len(x))
    b=np.zeros(n_nodi)
    #calcolo i local load vector bi con i=1,..n
    #e poi assemblo il load vector complessivo 
    for i in range(1,n_nodi):
        h=x[i]-x[i-1]
        b_i=np.zeros(n_nodi)
        b_i[i-1]=f(x[i-1])*h/2
        b_i[i]=f(x[i])*h/2
        b+=b_i    
    b= b[1:-1] #considero solo gli elementi interni!! Perché ho delle Dirichelet BC
    return b

def mesh_refinement(x,f,alpha=0.9):
    #adesso vado a fare il loop di mesh refinement per ... volte
    for loop in range(n_loop_refinement):
        #calcolo della soluzione
        A=A_global(x)
        b=b_load(x,f)
        uh=np.linalg.solve(A,b)

        # Aggiungo condizioni di Dirichlet 
        u = np.zeros(len(x))  # Vettore soluzione completo
        u[1:-1] = uh  # Inserisci la soluzione nei nodi interni
        #applico le Dirichelet BC
        u[0] = 0  # u(0) = 0
        u[-1] = 0  # u(1) = 0

        # 1)calcolare i vari element residual
        n_nodi=len(u)
        n_elementi=n_nodi-1
        eta=np.zeros(n_elementi)
        for i in range(1,n_nodi):
            n=i-1 #contatore nodi (i sarà invece il contatore degli elementi)
            z=i-1 #contatore indici di eta
            h=x[i]-x[i-1]
            eta_i=h*np.sqrt((f(x[i-1])+f(x[i]))/2 *h)
            eta[z]=eta_i
        
        # 2)criterio di refinement
        for i in range(1,n_nodi):
            n=i-1 #contatore nodi 
            z=i-1 #contatore indici di eta
            if eta[z] > alpha * np.max(eta):
                xmed=(x[i-1]+x[i])/2
                x=np.append(x,xmed)
                x=np.sort(x) #riordino i nodi 
    #calcolo la soluzione finale alla fine del loop        
    A=A_global(x)
    b=b_load(x,f)
    uh=np.linalg.solve(A,b)

    # Aggiungo condizioni di Dirichlet 
    u = np.zeros(len(x))  # Vettore soluzione completo
    u[1:-1] = uh  # Inserisci la soluzione nei nodi interni
    #applico le Dirichelet BC
    u[0] = 0  # u(0) = 0
    u[-1] = 0  # u(1) = 0
    return A,b,u,x,alpha


alpha=0.9
A,b,u,x,alpha= mesh_refinement(x,f)
print(f'con alpha={alpha} ho len(x)= ',len(x))

# Plot
plt.figure(figsize=(10, 4))

# Primo subplot (1 riga, 2 colonne, primo grafico)
plt.subplot(1, 2, 1)
plt.plot(x, u, label='u')
plt.scatter(x, u, c='red', marker='o', label='Nodi')
plt.title(f'Con alpha= {alpha}')
plt.legend()
plt.xlabel('x')
plt.ylabel('u')

#######################################################
#Vado a confrontare il risultato usando un valore diverso di alpha per il refinement
x = np.linspace(dominio[0], dominio[1], n_nodi_iniz)
alpha=1.2
A,b,u,x,alpha= mesh_refinement(x,f,alpha)
print(f'con alpha={alpha} ho len(x)= ',len(x))
        
plt.subplot(1, 2, 2)
plt.plot(x, u, label='u')
plt.scatter(x, u, c='red', marker='o', label='Nodi')
plt.title(f'Con alpha= {alpha}')
plt.legend()
plt.xlabel('x')
plt.ylabel('u')


plt.tight_layout()
plt.show()
