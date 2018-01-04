import math as math
import pandas as pd
import numpy as np
import sys as sys
def LU(A):
    n=len(A)
    L=pd.DataFrame(np.zeros([n,n]))
    U=pd.DataFrame(np.zeros([n,n]))
    P=pd.DataFrame(np.zeros([n,n],dtype=int))
    A[n]=pd.DataFrame([x for x in range(n)])
    for i in range(n):
        L.iloc[i,i]=1
        choose = i;
        for j in range(i+1,n):
            if A.iloc[choose,i]<abs(A.iloc[j,i]):
                choose=j
        if A.iloc[choose,i]==0:
            continue
        tmp=A.iloc[choose,:].copy()
        A.iloc[choose,:]=A.iloc[i,:]
        A.iloc[i,:]=tmp
        for j in range(i+1,n):
            A.iloc[j,i]=A.iloc[j,i]/A.iloc[i,i]
            A.iloc[j,i+1:n]=A.iloc[j,i+1:n]-A.iloc[j,i]*A.iloc[i,i+1:n]
    for i in range(1,n):
        L.iloc[i,0:i]=A.iloc[i,0:i]
    for i in range(n):
        U.iloc[i,i:n]=A.iloc[i,i:n]
    for i in range(n):
        P.iloc[i,int(A.iloc[i,n])]=1
    print('The resulf of LU Factorization:')
    print(' L:\n',L.values)
    print(' U:\n',U.values)
    print(' P:\n',P.values)

def QRGram(A):
    m,n=A.shape
    Q=pd.DataFrame(np.zeros([m,n]))
    R=pd.DataFrame(np.zeros([n,n]))
    for i in range(n):
        Q[i]=A[i]
        for j in range(i):
            r=0
            for k in range(m):
                r=r+Q.iloc[k,j]*A.iloc[k,i]
            R.iloc[j,i]=r
            Q[i]=Q[i]-r*Q[j]
        sum=0
        for j in range(m):
            sum=sum+Q.iloc[j,i]*Q.iloc[j,i]
        R.iloc[i,i]=math.sqrt(sum)
        if R.iloc[i,i] ==0:
            continue
        Q[i]=Q[i]/R.iloc[i,i]

    print('The result of QR with Gram-Schmidt')
    print(' Q:\n',Q.values)
    print(' R:\n',R.values)


def QRHouseholder(A):
    R = A.copy()
    m,n=R.shape
    Q=pd.DataFrame(np.eye(m))
    for i in range(min(m-1,n)):
        u=R.iloc[i:,i]-math.sqrt(sum(R.iloc[i:,i]*R.iloc[i:,i].T))*pd.DataFrame(np.eye(m)).iloc[i:,i]
        r=pd.DataFrame(np.eye(m))
        r.iloc[i:,i:]=pd.DataFrame(np.eye(m)).iloc[i:,i:]-2*pd.DataFrame(u).dot(pd.DataFrame(u).T)/sum(u.T*u)
        R=r.dot(R)
        Q=r.dot(Q)
    Q=Q.T
    print('The result of QR with Householder')
    print(' Q:\n',Q.values)
    print(' R:\n',R.values)
def QRGivens(A):
    R = A.copy()
    m, n = R.shape
    Q = pd.DataFrame(np.eye(m))
    for i in range(min(m - 1, n)):
        #P = pd.DataFrame(np.eye(m))
        for k in range(i+1,m):
            R[i][i]
            c=R[i][i]/math.sqrt(R[i][i]*R[i][i]+R[i][k]*R[i][k])
            s=R[i][k]/math.sqrt(R[i][i]*R[i][i]+R[i][k]*R[i][k])
            Pij=pd.DataFrame(np.eye(m))
            Pij.iloc[[i,k],[i,k]]=[[c,s],[-s,c]]
            R=Pij.dot(R)
            Q=Pij.dot(Q)
    Q = Q.T
    print('The result of QR with Givens')
    print(' Q:\n', Q.values)
    print(' R:\n', R.values)

def Factorization(A,type):
    if type == 0:
        LU(A)
    elif type==1:
        QRGram(A)
    elif type==2:
        QRHouseholder(A)
    elif type==3:
        QRGivens(A)
'''
A=[[1,2,-3,4],[4,8,12,-8],[2,3,2,1],[-3,-1,1,-4]]
B=[[0,-20,-14],[3,27,-4],[4,11,-2]]
LU(A)
QRGram(B)
QRHouseholder(B)
QRGivens(B)
'''
print ("脚本名：", sys.argv[0])
for i in range(1, len(sys.argv)):
    print ("参数", i, sys.argv[i])
A=pd.read_csv(sys.argv[1],header=None)
type=int(sys.argv[2])
print('待分解矩阵 A：\n',A.values)
Factorization(A,type)