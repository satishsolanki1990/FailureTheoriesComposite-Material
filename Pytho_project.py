# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 00:00:14 2019

@author: SATISH
"""

import math
import numpy as np
from numpy.linalg import inv
#import matplotlib.pyplot as plt
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
#import matplotlib

# Material and its Property selection function
# return  elastic Q11,Q22,Q12,Q66
def matQ(mat_key):
    Mat={1:[147e9,10.3e9,0.27,0.02,7e9],2:[169e9,9e9,0.31,0.02,6.5e9],3:[80e9,5.5e9,0.34,0.02,2.2e9],4:[45e9,11e9,0.29,0.06,4.5e9],5:[41e9,10.4e9,0.28,0.06,4.3e9]}
    E=Mat[mat_key]
    E1=E[0]
    E2=E[1]
    v12=E[2]
    v21=E[3]
    G12=E[4]
    
    Q11=E1/(1-v12*v21)
    Q22=E2/(1-v12*v21)
    Q12=v21*E1/(1-v12*v21)
    Q66=G12
    
    A=[Q11,Q12,0,Q12,Q22,0,0,0,Q66]
    Q_12=np.array([[A[0],A[1],A[2]],[A[3],A[4],A[5]],[A[6],A[7],A[8]]])
    return Q_12 

def Q_ij(mat_key,thita):
    Q_12=matQ(mat_key)
    Q11=Q_12[0][0]
    Q12=Q_12[0][1]
    Q22=Q_12[1][1]
    Q66=Q_12[2][2]

    n=math.sin(thita*math.pi/180)
    m=math.cos(thita*math.pi/180)
    #print('m,n:',m,n)
    Qxx=(m**4)*Q11 + (n**4)*Q22 + 2*(m**2)*(n**2)*Q12 + 4*(m**2)*(n**2)*Q66
    Qyy=(n**4)*Q11 + (m**4)*Q22 + 2*(m**2)*(n**2)*Q12+ 4*(m**2)*(n**2)*Q66
    Qxy=(m**2)*(n**2)*Q11 + (m**2)*(n**2)*Q22 + ((m**4)+(n**4))*Q12 - 4*(m**2)*(n**2)*Q66
    Qxs=(m**3)*n*Q11 - m*(n**3)*Q22 - m*n*((m**2)-(n**2))*Q12 - 2*m*n*((m**2)-(n**2))*Q66
    Qys=m*(n**3)*Q11 - (m**3)*n*Q22 + m*n*((m**2)-(n**2))*Q12 + 2*m*n*((m**2)-(n**2))*Q66
    Qss=(m**2)*(n**2)*Q11 + (m**2)*(n**2)*Q22 - 2*(m**2)*(n**2)*Q12 + (((m**2)-(n**2))**2)*Q66
    
    Q=[Qxx,Qxy,Qxs,Qxy,Qyy,Qys,Qxs,Qys,Qss]
    return Q


def symm(L):
    for i in range(len(L)):
        if L[i+1][0]!=L[len(L)-i][0]:
            return False
    return True
    
def bal(L):
    sum=0
    for i in range(len(L)):
        if (L[i+1][0]==0 or L[i+1][0]==90):
            sum+=0
        else:
            sum+=L[i+1][0]
    if sum==0:
        return True
    else:
        return False
    
def cp(L):
    for i in range(len(L)):
        if (L[i+1][0]==0 or L[i+1][0]==90):
            continue
        else:
            return False
    return True
    
def ABD(L,N,t,symm,bal,cp):
    z=[]
    A=[0,0,0,0,0,0,0,0,0]
    B=[0,0,0,0,0,0,0,0,0]
    D=[0,0,0,0,0,0,0,0,0]
    for i in range(N+1):
        z.append(-(N*t/2)+(i*t))
    for i in range(N):
        Q=Q_ij(L[i+1][2],L[i+1][0])
        for j in range(9):
            A[j]+=Q[j]*(z[i+1]-z[i])
            B[j]+=0.5*Q[j]*(z[i+1]**2-z[i]**2)
            D[j]+=(1/3)*Q[j]*(z[i+1]**3-z[i]**3)
    
    #A=[Axx,Axy,Axs,Ayx,Ayy,Ays,Asx,Asy,Ass]
    #B=[Bxx,Bxy,Bxs,Byx,Byy,Bys,Bsx,Bsy,Bss]
    #D=[Dxx,Dxy,Dxs,Dyx,Dyy,Dys,Dsx,Dsy,Dss]
    A=np.array([[A[0],A[1],A[2]],[A[3],A[4],A[5]],[A[6],A[7],A[8]]])
    B=np.array([[B[0],B[1],B[2]],[B[3],B[4],B[5]],[B[6],B[7],B[8]]])
    D=np.array([[D[0],D[1],D[2]],[D[3],D[4],D[5]],[D[6],D[7],D[8]]])
    return A,B,D


def Layout():
    #get input 
    t=float(input('thickness of each ply is(in mm)='))
    t=t*(10**-3)
    #print('please give lay out notation:')
    i=1
    L={}
    flag='y'
    while(flag=='y'):
        print('1. AS4/3501-6 Carbon-Epoxy\n2. IM6G/3501-6 Carbon-Epoxy\n3. Aramid 49-Epoxy\n4. S-Glass-Epoxy\n5. E-Glass-Epoxy')
        mat_key=int(input('The material of ply: '))
        ap=int(input('press 1 if angle of ply is in \u00B1\u03B8: '))
        angle=float(input('Angle of the ply: '))
        N=int(input('no of plies= '))
        if ap==1:
            for k in range(N*2):
                    L[i]=[(-1)**(k)*angle,1,mat_key]
                    i+=1
        else:
            for k in range(N):
                    L[i]=[angle,1,mat_key]
                    i+=1
        flag=input('press y to continue, n to stop:\n ')
    s=int(input('press:\n1. Symmetric with even no of plies\n2.symmetric with odd no of plies\n3.Not symmetric\n '))
    n=len(L)
    if s==1:
        for k in range(n):
            L[n+k+1]=[L[n-k][0],L[n-k][1],L[n-k][2]]
    if s==2:
        for k in range(n-1):
            L[n+k+1]=[L[n-1-k][0],L[n-1-k][1],L[n-1-k][2]]
    return L

# Calcualte Elastic Laminate Constants
def abcd(A,B,D):

    A1=inv(A)
    B_star=-np.matmul(A1,B)
    C_star=np.matmul(B,A1)
    D_star=D-np.matmul(C_star,B)
    
    a=inv(A) - np.matmul(np.matmul(B_star,inv(D_star)),C_star)
    b=np.matmul(B_star,inv(D_star))
    c=-np.matmul(inv(D_star),C_star)
    d=inv(D_star)
    
    return a,b,c,d

def Elastic_constant():
    L,N,t,symm,bal,cp=Layout()
    h=N*t
    A,B,D=ABD(L,N,t,symm,bal,cp)
    a,b,c,d=abcd(A,B,D)
    Ex=1/(h*a[0][0])
    Ey=1/(h*a[1][1])
    Gxy=1/(h*a[2][2])
    vxy=-a[1][0]/a[0][0]
    vyx=-a[0][1]/a[1][1]
    nsx=-a[0][2]/a[2][2]
    nxs=-a[2][0]/a[0][0]
    nsy=-a[1][2]/a[2][2]
    nys=-a[2][1]/a[1][1]

    return Ex,Ey,Gxy,vxy,vyx,nsx,nxs,nsy,nys



