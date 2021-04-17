# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:03:00 2019

@author: SATISH
"""

import numpy as np

def txt2dct():
    data={}
    with open('MatProp.txt','r') as f:
        for lines in f:
            key,val = lines.split(':')
            key = key.strip()
            val=val.strip('[]\n')
            val = val.replace(',''',' ')
            val=val.split()
            for i in range(len(val)-1):
                val[i+1]=float(val[i+1])
            data[int(key)] = val
    return data

def dct2txt(dct):
    with open('MatProp.txt','w') as f:
        for key in dct:
            f.write(str(key).rstrip('\n'))
            f.write(':'.rstrip('\n'))
            f.write(str(dct[key]))
            f.write('\n')
    return None

def Material():
    # Read/converting text file to dictionary
    matDict=txt2dct()
    # Need to Add new Material or not
    print('Material List format: \nKey: (name,Volfrac,E1, E2, G12, v12,v21, F1t,F2t,F6,F1c,F2c)')
    print('0: Add a New Material')
    try:
        for key in matDict:
            print('%d: %s' %(key, matDict[key]))
    except:
        exit()
    a=int(input('Choose Material or Add a material with property: '))
    
    #if yes to add new Material into database
    if a>len(matDict):
        print('incorrect entry, Start Again!!')
        raise ValueError('')
    elif a==0:           
        name=input('Enter Name of Lamina Material (eg. Carbon/Epoxy): ')
        VolumeFraction=float(input('Fiber Volume Fraction (between 0 and 1): '))
        if (VolumeFraction>1 or VolumeFraction<0):
            print('Fiber Volume fraction must be in (0,1)')
            raise ValueError('')

        E1=float(input('Longitudnal Modulus, E1, GPa (eg. 64.5): '))
        E2=float(input('Longitudnal Modulus, E1, GPa (eg. 14.5): '))
        G12=float(input('In-Plane Shear Modulus, G12, GPa (eg. 5.1): '))
        v12=float(input('Major Poison Ratio,v12 (eg. 0.27): '))
        if (v12>0.5 or v12<0):
            print('Poison Ratio is not correct')
            raise ValueError('')
            
        v21=float(input('Minor Poison Ratio,v21 (eg. 0.02): '))
        if (v21>0.5 or v21<0):
            print('Poison Ratio is not correct')
            raise ValueError('')
            
        F1t=float(input('Longitudnal Tensile Strength, F1t, MPa (eg. 3250): '))
        F2t=float(input('Transverse Tensile Strength, F2t, MPa (eg. 62): '))
        F6=float(input('In-Plane Shear Strength, F6, MPa(eg. 75): '))
        F1c=float(input('Longitudnal Compressive Strength, F1c, MPa (eg. 1500): '))
        F2c=float(input('Transverse Compressive Strength, F2c, MPa (eg. 200): '))
        matDict[len(matDict)+1]=[name,VolumeFraction,E1*1e9,E2*1e9,G12*1e9,v12,v21,F1t*1e6,F2t*1e6,F6*1e6,F1c*1e6,F2c*1e6]
        mat=matDict[len(matDict)]
        # save/write back dictionary to text File
        dct2txt(matDict)        
    #Else slected material
    else:
        mat=matDict[a]
    return mat
   
def Compliance(mat,thita):

    Q11=mat[2]/(1-mat[5]*mat[6])
    Q22=mat[3]/(1-mat[5]*mat[6])
    Q12=mat[6]*mat[2]/(1-mat[5]*mat[6])
    Q66=mat[4] 
   
    n=np.sin(thita*np.pi/180)
    m=np.cos(thita*np.pi/180)
    Qxx=(m**4)*Q11 + (n**4)*Q22 + 2*(m**2)*(n**2)*Q12 + 4*(m**2)*(n**2)*Q66
    Qyy=(n**4)*Q11 + (m**4)*Q22 + 2*(m**2)*(n**2)*Q12+ 4*(m**2)*(n**2)*Q66
    Qxy=(m**2)*(n**2)*Q11 + (m**2)*(n**2)*Q22 + ((m**4)+(n**4))*Q12 - 4*(m**2)*(n**2)*Q66
    Qxs=(m**3)*n*Q11 - m*(n**3)*Q22 - m*n*((m**2)-(n**2))*Q12 - 2*m*n*((m**2)-(n**2))*Q66
    Qys=m*(n**3)*Q11 - (m**3)*n*Q22 + m*n*((m**2)-(n**2))*Q12 + 2*m*n*((m**2)-(n**2))*Q66
    Qss=(m**2)*(n**2)*Q11 + (m**2)*(n**2)*Q22 - 2*(m**2)*(n**2)*Q12 + (((m**2)-(n**2))**2)*Q66
    
    Q=np.array([[Qxx,Qxy,Qxs],[Qxy,Qyy,Qys],[Qxs,Qys,Qss]])
    
    return Q

    
def sxy_to_12(stress_xy,thita):
    n=np.sin(thita*np.pi/180)
    m=np.cos(thita*np.pi/180)
    T=np.array([[m**2,n**2,2*m*n],[n**2,m**2,-2*m*n],[-m*n,m*n,(m**2-n**2)]])
    s12=np.matmul(T,stress_xy)
    return s12
    
def strain_to_stress(mat,strain,thita):
    Q=Compliance(mat,thita)
    stress=Q*strain
    return stress

    
def max_stress(stress_12,mat):
    Sf=[]
    # sigma 1
    if stress_12[0]>0:
        Sf.append(mat[7]/stress_12[0])
    if stress_12[0]<0:
        Sf.append(-mat[10]/stress_12[0])
        
    #sigma 2
    if stress_12[1]>0:
        Sf.append(mat[8]/stress_12[1])
    if stress_12[1]<0:
        Sf.append(-mat[11]/stress_12[1])        
    # Tau 6
    if stress_12[2]!=0:
        Sf.append(mat[9]/abs(stress_12[1]))
    return min(Sf)

def max_strain(stress_12,mat):
    Sf=[]
    # sigma 1
    e1=stress_12[0] - (mat[5]*stress_12[1])
    if e1>0:
        Sf.append(mat[7]/e1)
    if e1<0:
        Sf.append(-mat[10]/e1)
        
    #sigma 2
    e2=stress_12[1] - (mat[6]*stress_12[0])
    if e2>0:
        Sf.append(mat[8]/e2)
    if e2<0:
        Sf.append(-mat[11]/e2)        
    
    # Tau 6
    if stress_12[2]!=0:
        Sf.append(mat[9]/abs(stress_12[1]))
    return min(Sf)

def tsai_hill(stress_12,mat):
#Sf is either +1.414 or -1.414
    a=(stress_12[0]/mat[7])**2+(stress_12[1]/mat[8])**2+(stress_12[2]/mat[9])**2-(stress_12[0]*stress_12[1]/mat[7]**2)
    b=(stress_12[0]/mat[10])**2+(stress_12[1]/mat[8])**2+(stress_12[2]/mat[9])**2-(stress_12[0]*stress_12[1]/mat[10]**2)
    c=(stress_12[0]/mat[10])**2+(stress_12[1]/mat[11])**2+(stress_12[2]/mat[9])**2-(stress_12[0]*stress_12[1]/mat[10]**2)
    d=(stress_12[0]/mat[7])**2+(stress_12[1]/mat[11])**2+(stress_12[2]/mat[9])**2-(stress_12[0]*stress_12[1]/mat[7]**2)

    if (stress_12[0]>0 and stress_12[1]>0):
        if a<=1:
            Sf=1.414
        else:
            Sf=-1.414
    elif (stress_12[0]<0 and stress_12[1]>0):
        if b<=1:
            Sf=1.414
        else:
            Sf=-1.414
    elif (stress_12[0]<0 and stress_12[1]<0):
        if c<=1:
            Sf=1.414
        else:
            Sf=-1.414
    else:
        if d<=1:
            Sf=1.414
        else:
            Sf=-1.414
    return Sf

def tsai_wu(stress_12,mat):
    f1=(1/mat[7])-(1/mat[10])
    f2=(1/mat[8]) - (1/mat[11])
    f11=1/(mat[7]*mat[10])
    f22=1/(mat[8]*mat[11])
    f66=1/mat[9]**2
    f12=-0.5*np.sqrt(f11*f22)    
    a=f1*stress_12[0]+f2*stress_12[1]+f11*stress_12[0]**2+f22*stress_12[1]**2+f66*stress_12[2]**2+2*f12*stress_12[0]*stress_12[1]
    if a<=1:
        Sf=1
    else:
        Sf=0
    return Sf

def lamina_failure():
    # get material
    mat=Material()
    # get fiber orintation
    thita=float(input('Enter fiber Orintation in Lamina ( In Degree)'))
    # get strate of stress or strain
    print('What do you have:\n\t1.State of Stress(sigma_x,sigma_y,tau_6)\n\t2.State of Strain(strain in x,y,xy)')
    k=int(input(''))
    if k==1:
        a=float(input('sigma_x(pa)= '))
        b=float(input('sigma_y(pa)= '))
        c=float(input('Tau_xy(pa)= '))
        stress_xy=np.array([[a],[b],[c]])
        stress_12=sxy_to_12(stress_xy,thita)
    elif k==2:
        a=float(input('epsilon_x= '))
        b=float(input('epsilon_y= '))
        c=float(input('gama_xy(pa)= '))
        strain_xy=np.array([[a],[b],[c]])
        stress_xy=strain_to_stress(mat,strain_xy,thita)
        stress_12=sxy_to_12(stress_xy,thita)
    else:
        print('incorrect entry. Start again!!')
        raise ValueError('')

    # applie theories and get Sf from all
    Sf=[]
    Sf.append(max_stress(stress_12,mat))
    Sf.append(max_strain(stress_12,mat))
    Sf.append(tsai_hill(stress_12,mat))
    Sf.append(tsai_wu(stress_12,mat))

    # print massage based onf Sf for individual
    Theory=['Max Stress','Max Strain','Tsai Hill','Tsai Wu']
    a='y'
    while(a=='y'):
        print('\nWhich Theory you want to apply:\n\t1.Max Stress\n\t2.Max Strain\n\t3.Tsai-Hill\n\t4.Tsai-Wu\n\t5. All Together')
        k=input('')
        if k=='5':
            if min(Sf)>=1:
                print('Lamina is safe under according to all theories')
            else:
                for i in range(len(Sf)):
                    if Sf[i]<1:
                        print('Fail by %s Theory'%Theory[i])
                    else:
                        print('Pass by %s Theory'%Theory[i])
        elif k=='4':
            if Sf[-1]>=1:
                print('Lamina is safe under given state of stress by %s theory'%Theory[-1])
            else:
                print('Fail by %s Theory'%Theory[-1])
        elif k=='3':
            if Sf[-2]>=1:
                print('Lamina is safe under given state of stress by %s theory'%Theory[-2])
            else:
                print('Fail by %s Theory'%Theory[-2])
        elif k=='2':
            if Sf[-3]>=1:
                print('Lamina is safe under given state of stress by %s theory'%Theory[-3])
            else:
                print('Fail by %s Theory'%Theory[-3])
        elif k=='1':
            if Sf[-4]>=1:
                print('Lamina is safe under given state of stress by %s theory'%Theory[-4])
            else:
                print('Fail by %s Theory'%Theory[-4])
        else:
            print('Incorrect entry. Start again!!')
            raise ValueError('')
        
        a=input('\nDo you want to try another Theory (y/n): ')
        
    return None

a='y'
while (a=='y'):
    try:
        lamina_failure()
        a=input('\nContinue for another Lamina analysis (y/n):  ')
    except:
        a=input('\nContinue for another Lamina analysis (y/n):  ')
        
input('\n\nThank you for using This Software! Press any key to exit.')

