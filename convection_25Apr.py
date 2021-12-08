#slug flow problem
#In discrete form convection diffusion equation equation is given as
# Unsteady term + Convection term = Diffusion term + Source term
# Use combined Fourier and CFL criteria to find del_t

import numpy as np
import scipy as sp
from numpy import *
from math import pi

#Material Properties
km = float(input("What is the thermal conductivity of material?(W/mK) \n >"))
rho = float(input("What is the density of the matierial?(kg/m^3) \n >"))
Cp = float(input("What is the specific heat of the material?(J/kgK) \n >"))
qg = float(input("Enter the internal heat generation rate(W/m3) \n >"))
#km = 69.3; rho = 8050 ; Cp = 350; qg=0
#Geometrical properties

L = float(input("What is the length of channel? \n >")) #Length is along x-direction
W = float(input("What is the width of the channel? \n >")) #Width is along y-direction
#L = 0.5; W = 0.25

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
eps = float(input("Enter convergence criteria."))
alp = km/(rho*Cp)

#Flow Properties
u = float(input("what is the flow veocity in x-direction u ="))
v = float(input("What is the flow velocity in y-direction v = "))
Ti = 1; Tw = 0
dtd = round(0.5/(alp*(1/(dx**2) + 1/(dy**2))),3)
print(f"dtd={dtd} secs")
dtc = round(0.9/(u/dx + v/dy),3)
print(f"dtc={dtc} secs")
dtu = float(input("Your custom time step other than this="))
dt = min(dtd,dtc,dtu)

Fmw = rho*u ; Fme = Fmw
Fms = rho*v ; Fmn = Fms

#Initial and boundary conditions
phi = np.empty((imax+1,jmax+1));  phi_old = np.empty((imax+1,jmax+1))
phx = np.empty((imax+2,jmax+2))
phy = np.empty((imax+2,jmax+2))
C = np.empty((imax+1,jmax+1)); D = np.empty((imax+1,jmax+1))
S = np.empty((imax+1,jmax+1))
#Initializing solution
phi_old[:,:] = 0; phx[:,:] = 0; phy[:,:] = 0
C[:,:] = 0; D[:,:] = 0 ; S[:,:] = 0;
phi[:,:] = 0 #initial condition
#Boundary conditions
phx[1,:] = Ti
for j in range(1,jmax):
    phx[imax+1,j] = phi[imax,j]
phy[:,1] = Tw
phy[:,jmax+1] = phy[:,1]

scheme = int(input("Which convection scheme you want?\n1)FOU 2)CD\n3)SOU\n4)QUICK\n>"))

def diffusion(i,j):
    if i == 1 and j == 1: #Southwest corner
        D[i,j] = km*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phx[i,j])/(dx**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phy[i,j])/(dy**2))
    elif i == imax and j == jmax: #Northeast corner
        D[i,j] = km*((2*phx[i+1,j]-3*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (2*phy[i,j+1]-3*phi_old[i,j]+ phi_old[i,j-1])/(dy**2)) # Was I very stupid to take (i-1) instead j-1!
    elif i == 1 and j == jmax: #Northwest corner
        D[i,j] = km*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phx[i,j])/(dx**2) + (2*phy[i,j+1]-3*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
    elif j == 1 and i== imax: #Southeast corner
        D[i,j] = km*((2*phx[i+1,j]-3*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phy[i,j])/(dy**2))

    elif i == 1 and (j>1 and j < jmax): # Left Boundary
        D[i,j] = km*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phx[i,j])/(dx**2) + (phi_old[i,j+1]-2*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
    elif i == imax and (j>1 and j < jmax): #Right Boundary
        D[i,j] = km*((2*phx[i+1,j]-3*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-2*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
    elif j == 1 and (i>1 and i < imax): #South Boundary
        D[i,j] = km*((phi_old[i+1,j]-2*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phy[i,j])/(dy**2))
    elif j == jmax and (i>1 and i < imax): #North Boundary
        D[i,j] = km*((phi_old[i+1,j]-2*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (2*phy[i,j+1]-3*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))

    if (i>1 and i < imax) and (j>1 and j<jmax): #Interior
        D[i,j] = km*((phi_old[i+1,j]-2*phi_old[i,j]+ phi_old[i-1,j])/(dx**2) + (phi_old[i,j+1]-2*phi_old[i,j]+ phi_old[i,j-1])/(dy**2))
        #print(f"Values {i} {j} {dt*D[i,j]/(rho*Cp)}")
    return(D[i,j])

n=0
time = 0
rmsphi = 0.5
while rmsphi/dt > eps:
    n = n+1
    if scheme == 1:
        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                phi_old[i,j] = phi[i,j]
                phx[1,j] = Ti
                phx[imax+1,j] = phi[imax,j]
                phy[i,jmax+1] = Tw
                phy[i,1] = phy[i,jmax+1] #Boundary conditions update
                if i<imax:
                    phx[i+1,j] = phi[i,j]##FOU Scheme x-direction array
                if j<jmax:
                    phy[i,j+1] = phi[i,j] ##FOU scheme y-direction array

                C[i,j] =  Cp*((Fme* phx[i+1,j]- Fmw* phx[i,j])*dy+(Fmn* phy[i,j+1]- Fms* phy[i,j])*dx)/dv
                D[i,j] = diffusion(i,j)
                S[i,j] =  qg
                phi[i,j] = phi_old[i,j] + (dt/(rho*Cp))*(D[i,j]-C[i,j] +S[i,j])
        sum = 0 #Initial value was carried to next iteration as it is and hence my error was remaining constant after so many iterations
        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                sum = sum + ((phi[i,j]-phi_old[i,j])**2)

    elif scheme == 2:
        for j in range(1,jmax+1):
            for i in range(1,imax+1):
                phi_old[i,j] = phi[i,j]
                phx[1,j] = Ti
                phx[imax+1,j] = phi[imax,j]
                phy[i,jmax+1] = Tw
                phy[i,1] = phy[i,jmax+1] #Boundary conditions update

                if i<imax:
                    phx[i+1,j] = (phi[i+1,j]+ phi[i,j])/2 #East face CD Scheme x-direction
                if j<jmax:
                    phy[i,j+1] = (phi[i,j+1]+ phi[i,j])/2 #North face CD Scheme y-direction

                C[i,j] =  Cp*((Fme* phx[i+1,j]- Fmw* phx[i,j])*dy+(Fmn* phy[i,j+1]- Fms* phy[i,j])*dx)/dv
                D[i,j] = diffusion(i,j)
                S[i,j] =  qg
                phi[i,j] = phi_old[i,j] + dt/(rho*Cp)*(D[i,j] -C[i,j] +S[i,j])
        sum =  0 #Initial value was carried to next iteration as it is and hence my error was remaining constant after so many iterations
        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                sum = sum + ((phi[i,j]- phi_old[i,j])**2)

    elif scheme ==3:
        #East face of left border cells
        #North face of bottom border cells

        for i in range(1,imax+1):
            for j in range(1, jmax+1):
                phi_old[i,j] = phi[i,j]
                phx[1,j] = Ti
                phx[imax+1,j] = phi[imax,j]
                phy[i,jmax+1] = Tw
                phy[i,1] = phy[i,jmax+1] #Boundary conditions update
                #Convection scheme of Second order upwind (SOU)
                phx[2,j] = 2*phi[1,j] - phx[1,j]
                phy[i,2] = 2*phi[i,1] - phy[i,1]
                if i>1 and i< imax:
                    phx[i+1,j] = (3*phi[i,j]- phi[i-1,j])/2 #East face
                if j>1 and j< jmax:
                    phy[i,j+1] = (3*phi[i,j]- phi[i,j-1])/2 #North face
                C[i,j] =  Cp*((Fme* phx[i+1,j]- Fmw* phx[i,j])*dy+(Fmn* phy[i,j+1]- Fms* phy[i,j])*dx)/dv
                D[i,j] = diffusion(i,j)
                S[i,j] =  qg
                phi[i,j] = phi_old[i,j] + dt/(rho*Cp)*(D[i,j] -C[i,j] +S[i,j])
        sum =  0 #Initial value was carried to next iteration as it is and hence my error was remaining constant after so many iterations
        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                sum = sum + ((phi[i,j]- phi_old[i,j])**2)

    elif scheme ==4:
        for i in range(1,imax+1):
            for j in range(1, jmax+1):
                phi_old[i,j] = phi[i,j]
                phx[1,j] = Ti
                phx[imax+1,j] = phi[imax,j]
                phy[i,jmax+1] = Tw
                phy[i,1] = phy[i,jmax+1] #Boundary conditions update
                #convection scheme of QUICK
                phx[2,j] = (1/3*phi[2,j]+phi[1,j]-1/3*phx[1,j]) #East face of left border cells
                phy[i,2] = (1/3*phi[i,2]+phi[i,1]-1/3*phy[i,1]) #North face of bottom border cells
                if i>1 and i<imax:
                    phx[i+1,j] = (3*phi[i+1,j]+6*phi[i,j]-phi[i-1,j])/8 #East face
                if j>1 and j<jmax:
                    phy[i,j+1] = (3*phi[i,j+1]+6*phi[i,j]-phi[i,j-1])/8 #North face
                C[i,j] =  Cp*((Fme* phx[i+1,j]- Fmw* phx[i,j])*dy+(Fmn* phy[i,j+1]- Fms* phy[i,j])*dx)/dv
                D[i,j] = diffusion(i,j)
                S[i,j] =  qg
                phi[i,j] = phi_old[i,j] + dt/(rho*Cp)*(D[i,j] -C[i,j] +S[i,j])
        sum =  0 #Initial value was carried to next iteration as it is and hence my error was remaining constant after so many iterations
        for j in range(1, jmax+1):
            for i in range(1, imax+1):
                sum = sum + ((phi[i,j]- phi_old[i,j])**2)

    rmsphi = np.sqrt(sum/(imax*jmax))
    time = round(time + dt,3)
    print(f"Calculating phi after {time}secs.")
    print(f"rms error for {n}th iteration is {rmsphi}")
    if n%50 ==0:
        if scheme == 1:
            target = open(f"conv_diff_FOU_{time}.dat",'w')
            target.write("VARIABLES = X, Y, phi\n")
            target.write(f"ZONE T = t_{time}\n")
        elif scheme ==2:
            target = open(f"cov_diff_CD_{time}.dat",'w')
            target.write("VARIABLES = X, Y, phi\n")
            target.write(f"ZONE T = t_{time}\n")
        elif scheme ==3:
            target = open(f"conv_diff_SOU_{time}.dat", 'w')
            target.write("VARIABLES = X, Y, phi\n")
            target.write(f"ZONE T = t_{time}\n")
        elif scheme ==4:
            target = open(f"conv_diff_QUICK_{time}.dat", 'w')
            target.write("VARIABLES = X Y phi\n")
            target.write(f"ZONE T = t_{time}\n")

        for j in range(1,jmax+1):
            target.write(f"0 {(j-1/2)*dy} {phx[1,j]}\n")
            target.write(f"{L} {(j-1/2)*dy} {phx[imax+1,j]}\n")
        for i in range(1,imax+1):
            target.write(f"{(i-1/2)*dx} 0 {phy[i,1]}\n")
            target.write(f"{(i-1/2)*dx} {W} {phy[i,jmax+1]}\n")

        for i in range(1,imax+1):
            for j in range(1,jmax+1):
                target.write(f"{(i-1/2)*dx} {(j-1/2)*dy} {phi[i,j]}\n")
