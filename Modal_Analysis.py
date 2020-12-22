## PROGRAM FOR DYNAMIC ANALYSIS OF 1D BAR SUBJECTED TO UNDAMPED FREE VIBRATIONS ##
# import numpy 


import numpy as np
import pandas as pd


#Enter Input Data like number of elements,Length,E,density(specific weight),Cross-section Area. 
e=int(input("Enter number of elements : "))
n = e+1
print("Number of nodes : ", n)
K=np.zeros((n, n),dtype=float)
M=np.zeros((n,n),dtype=float)
btyp=input("Enter type of bar ('s'for stepped bar,'u' for uniform bar ) : ")
if btyp in ['s','S']:
    for i in range(e):
        l=float(input("Enter length of element %d : "%(i+1)))
        A=float(input("Enter cross-section area for element %d : "%(i+1)))
        E=float(input("Enter Young's modulus for element %d : "%(i+1)))
        rho=float(input("Enter material density for element %d : "%(i+1)))
        C=(A*E)/l
        D=(rho*A*l)/6
        mat1=np.matrix([[1,-1],[-1,1]])
        mat2=np.matrix([[2,1],[1,2]])
        k1=np.multiply(C,mat1)              # Element stiffness matrix for bar element
        print("\n")
        print("Elemental Stiffness matrix for element %d" % (i+1))
        print(k1)
        print("\n")
        K[i:i+2,i:i+2]+=k1
        m1=np.multiply(D,mat2)
        print("Elemental Mass matrix for element %d" % (i+1))      # Element mass matrix for bar element
        print(m1)
        print("\n")
        M[i:i+2, i:i+2] += m1
elif btyp in ['u','U']:
    A=float(input("Enter cross section area of uniform bar : "))
    E=float(input("Enter Young's modulus for uniform bar : "))
    rho=float(input("Enter material density for uniform bar : "))
    for o in range(e):
        l=float(input("Enter length of element %d : "%(o+1)))
        C=(A*E)/l
        D=(rho*A*l)/6
        mat1=np.matrix([[1,-1],[-1,1]])
        mat2=np.matrix([[2,1],[1,2]])
        k1=np.multiply(C,mat1)
        print("\n")
        print("Elemental Stiffness matrix for element %d"%(o+1))
        print(k1)
        print("\n")
        K[o:o+2,o:o+2]+=k1
        m1=np.multiply(D,mat2)
        print("Elemental Mass matrix for element %d" % (o+1))
        print(m1)
        print("\n")
        M[o:o+2,o:o+2]+=m1

# Global Stiffness And Global Mass Matrices 
print("\n")
print("------------------------------Global Stiffness matrix (K)--------------------------------")
Kd=pd.DataFrame(K[:,:],index=range(1,n+1),columns=range(1,n+1))
print(Kd)
print("\n")
print("------------------------------Global Mass matrix (M)-------------------------------------")
Md=pd.DataFrame(M[:,:],index=range(1,n+1),columns=range(1,n+1))
print(Md)
print("\n")
redsmat=[]
redmmat=[]
sup=int(input("Enter '1 ' if only one side is fixed \nEnter '2' if both sides are fixed\n"))
if sup==1:
    kd=Kd.drop([1],axis=0)
    md= Md.drop([1], axis=0)
    Kr=kd.drop([1],axis=1)
    gsd=Kr
    Mr=md.drop([1],axis=1)
    gmd=Mr
    redsmat=gsd.to_numpy()
    redmmat=gmd.to_numpy()
elif sup==2:
    kd=Kd.drop([1,n],axis=0)
    md=Md.drop([1,n], axis=0)
    Kr=kd.drop([1,n],axis=1)
    Mr=md.drop([1,n], axis=1)
    gsd=Kr
    gmd=Mr
    redsmat=gsd.to_numpy()
    redmmat=gmd.to_numpy()

# To Find Natural Frequency 
KINV=np.linalg.inv(redsmat)
EIG=np.matmul(KINV,redmmat)
val=np.linalg.eigvals(EIG)
omega=np.sqrt(1/val)
freq=omega/(2*np.pi)
print("\n-----------------------Natural Frequency (Hz)---------------------")
f=pd.DataFrame(data=freq,index=range(1,len(freq)+1),columns=["Frequency"])
print(f)
