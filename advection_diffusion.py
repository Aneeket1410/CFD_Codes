#slug flow problem
#In discrete form convection diffusion equation equation is given as
# Unsteady term + Convection term = Diffusion term + Source term
# Use combined Fourier and CFL criteria to find del_t

import numpy as np
import scipy as sp


from numpy import *
from math import pi

#Material Properties
#km = float(input("What is the thermal conductivity of material?(W/mK) \n >"))
#rho = float(input("What is the density of the matierial?(kg/m^3) \n >"))
#Cp = float(input("What is the specific heat of the material?(J/kgK) \n >"))
#qg = float(input("Enter the internal heat generation rate(W/m3) \n >"))
km = 50; rho = 1000 ; Cp = 4187; qg=0
#Geometrical properties
#T_0 = input("What is the reference temperature in deg Cel.? \n >")
#T_0 = float(T_0)

#L = float(input("What is the length of channel? \n >")) #Length is along x-direction
#W = float(input("What is the width of the channel? \n >")) #Width is along y-direction
#t = float(input("What is the thickness of the plate? \n >"))
L = 0.5; W = 0.25

#Inputs required for grid generation
print("For better results select cell dimensions so that you will get around 4000 cells.")

dx = float(input("What is cell length? \n >"))
dy= float(input("What is cell width? \n >"))
dv = dx*dy
cells = int((W*L)/(dx*dy))
print(f"You have {cells} number of cells in grid.")
imax = int(L/dx)
jmax = int(W/dy)
#eps = float(input("Enter the error criteria epsilon = "))
eps = 0.0001
alp = km/(rho*Cp)

#Flow Properties
u = float(input("what is the flow veocity in x-direction u ="))
v = float(input("What is the flow velocity in y-direction v = "))
Ti = 50; Tw = 100
dtd = 0.4/(alp*(1/dx**2 + 1/dy**2))
print(dtd)
dtc = 0.72/(u/dx + v/dy)
print(dtc)
dt = min(dtd,dtc)

print(dt)

Fmw = rho*u ; Fme = Fmw
Fms = rho*v ; Fmn = Fms

#Initial and boundary conditions
phi = np.empty((imax+1,jmax+1));  phi_old = np.empty((imax+1,jmax+1))
phe = np.empty((imax+1,jmax+1)); phw = np.empty((imax+1,jmax+1))
phn = np.empty((imax+1,jmax+1)); phs = np.empty((imax+1,jmax+1))
C = np.empty((imax+1,jmax+1)); D = np.empty((imax+1,jmax+1))
S = np.empty((imax+1,jmax+1))
phi_old[:,:] = 0; phe[:,:] = 0; phw[:,:] = 0; phs[:,:] = 0; phn[:,:] = 0
C[:,:] = 0; D[:,:] = 0 ; S[:,:] = 0;
phi[:,:] = Ti #initial condition
#Boundary conditions
phw[1,:] = Ti
phe[imax,:] = phi[imax,:]
phs[:,1] = Tw
phn[:,jmax] = Tw

#scheme = int(input("Which convection scheme you want?\n1)FOU 2)CD\n3)SOU\n4)QUICK\n>"))
scheme =1
def diffusion(i,j):





    if i == 1 and j == 1: #Southwest corner
        D[i,j] = km*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phw[i,j])/(dx**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phs[i,j])/(dy**2))
    elif i == 1 and (j>1 and j < jmax): # Left Boundary
        D[i,j] = km*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phw[i,j])/(dx**2) + (phi_old[i,j+1]-2*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
    elif i == 1 and j == jmax: #Northwest corner
        D[i,j] = km*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phw[i,j])/(dx**2) + (2*phn[i,j]-3*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
        #print(f"Values {i} {j} {dt*D[i,j]/(rho*Cp)}")
    elif j == 1 and (i>1 and i < imax): #South Boundary
        D[i,j] = km*((phi_old[i+1,j]-2*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phs[i,j])/(dy**2))
    elif j == 1 and i== imax: #Southeast corner
        D[i,j] = km*((2*phe[i,j]-3*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phs[i,j])/(dy**2))
    elif j == jmax and (i>1 and i < imax): #North Boundary
        D[i,j] = km*((phi_old[i+1,j]-2*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (2*phn[i,j]-3*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
    elif i == imax and (j>1 and j < jmax): #Right Boundary
        D[i,j] = km*((2*phe[i,j]-3*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-2*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
    elif i == imax and j == jmax: #Northeast corner
        D[i,j] = km*((2*phe[i,j]-3*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (2*phn[i,j]-3*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))

    if (i>1 and i < imax) and (j>1 and j<jmax): #Interior
        D[i,j] = km*((phi_old[i+1,j]-2*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-2*phi_old[i,j]+phi_old[i,j-1])/(dy**2))
        #print(f"Values {i} {j} {dt*D[i,j]/(rho*Cp)}")
    return(D[i,j])

n=0
time = 0
sum = 0
rmsphi = 0.5
while rmsphi/dt > eps:
    n = n+1
    if scheme == 1:
        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                phi_old[i,j] = phi[i,j]
                phw[1,j] = Ti; phe[1,j] = phi[1,j]
                phw[imax,j] = phi[imax-1,j]; phe[imax,j] = phi[imax,j]
                phs[i,1] = Tw; phn[i,1] = phi[i,1]
                phs[i,jmax] = phi[i,jmax-1]; phn[i,jmax] = Tw

                if i>1 and (i< imax):
                    phw[i,j] = phi[i-1,j]; phe[i,j] = phi[i,j]
                if j>1 and (j<jmax):
                    phs[i,j] = phi[i,j-1]; phn[i,j] = phi[i,j]

                C[i,j] =  ((Fme* phe[i,j]- Fmw* phw[i,j])*dy+(Fmn* phn[i,j]- Fms* phs[i,j])*dx)/dv
                D[i,j] = diffusion(i,j)
                S[i,j] =  qg
                phi[i,j] = phi_old[i,j] + (dt/(rho*Cp))*(D[i,j]-C[i,j] +S[i,j])
                sum = sum + ((phi[i,j]-phi_old[i,j])**2)

    elif scheme == 2:
        for j in range(1,imax+1):
            for i in range(1,jmax+1):
                phe[imax,j] = phi[imax,j]
                phw[1,j] = Ti
                phs[i,1] = Tw
                phn[i,jmax] = Tw
                phi_old[i,j] = phi[i,j]
                if i>1 and i<imax:
                    phe[i,j] = (phi[i+1,j]+ phi[i,j])/2 #East face
                    phw[i,j] = (phi[i-1,j]+ phi[i,j])/2 #West face

                if j>1 and j<jmax:
                    phs[i,j] = (phi[i,j-1]+ phi[i,j])/2 #South face
                    phn[i,j] = (phi[i,j+1]+ phi[i,j])/2 #North face

                C[i,j] =  Cp*((Fme* phe[i,j] * dy)+(Fmn* phn[i,j] *dx)-(Fmw* phw[i,j] * dy)-(Fms* phs[i,j] *dx))/dv
                D[i,j] = diffusion(i,j)
                S[i,j] =  qg
                phi[i,j] = phi_old[i,j] + dt/(rho*Cp)*(D[i,j] -C[i,j] +S[i,j])
                sum = sum + ((phi[i,j]-phi_old[i,j])**2)

    rmsphi = np.sqrt(sum/(imax*jmax))
    time = time + dt
    print(f"Calculating phi after {time}secs.")
    print(f"rms error for {n}th iteration is {rmsphi}")
    if n%100 ==0:
        if scheme == 1:
            target = open(f"conv_diff_FOU_{time}.dat",'w')
            target.write("VARIABLES = X, Y, phi\n")
        elif scheme ==2:
            target = open(f"cov_diff_CD_{time}.dat",'w')
            target.write("VARIABLES = X, Y, phi\n")
        elif scheme ==3:
            target = open(f"conv_diff_SOU_{time}.txt", 'w')
            target.write("VARIABLES = X, Y, phi\n")
        elif scheme ==4:
            target = open(f"conv_diff_QUICK_{time}.txt", 'w')
            target.write("VARIABLES = X Y phi\n")

        for j in range(1,jmax+1):
            target.write(f"0 {(j-1/2)*dy} {phw[1,j]}\n")
            target.write(f"{L} {(j-1/2)*dy} {phe[imax,j]}\n")
        for i in range(1,imax+1):
            target.write(f"{(i-1/2)*dx} 0 {phs[i,1]}\n")
            target.write(f"{(i-1/2)*dx} {W} {phn[i,jmax]}\n")

        for i in range(1,imax+1):
            for j in range(1,jmax+1):
                target.write(f"{(i-1/2)*dx} {(j-1/2)*dy} {phi[i,j]}\n")
