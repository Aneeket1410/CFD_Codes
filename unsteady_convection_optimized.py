import numpy as np
from numpy import *

#Initial conditions
cells =  1600
L,W = 1,1
rho = 1
u, v = 0.7071, 0.7071

#Grid parameters
del_x = np.sqrt(L*W/cells)
del_y = del_x
del_v = del_x * del_y
imax =  int(2*L/del_x +1)
jmax = int(2*W/del_y +1)

#Solution variables matrix
phi = np.empty((imax+1,jmax+1))
phi_old = np.empty((imax+1,jmax+1))
Fm = np.empty((imax+1,jmax+1))
C = np.empty((imax+1,jmax+1))
#Fmr = np.empty((imax+1,jmax+1))

phi[:,:] = 0.5
Fm[:,:] = 0.5
C[:,:] = 0.5
#Fmr[:,:] = 0.0001
#Stability criteria
#del_t = 0.8/((u/del_x) + (v/del_y))
del_t = 0.001


#Defining flux phi as initial condition
for i in range(2,imax):
    for j in range(2,jmax):
        phi[i,j] = 0
#for i in range(2,imax):
    #jc = i
    #if jc%2 == 0:
        #for j in range(1,jc+1): #flux on odd j and north or south face
            #phi[i,j] = 0
    #else:
        #for j in range(2,jc,2): #flux on even j east or west face of the cell
            #phi[i,j] = 0

#for j in range(2,jmax):
    #ic = j
    #if ic%2 == 0:
        #for i in range(1,ic+1): #flux on even j east or west face of the cell
            #phi[i,j] = 1
    #else:
        #for i in range(2,ic,2): #flux on odd j and north or south face
            #phi[i,j] = 1
#Boundary conditions
for j in range(1,jmax+1): #flux on right and left boundary
    phi[imax,j] = 0
    phi[1,j] = 1
for i in range(1,imax+1): #flux on top and bottom boundary
    phi[i,1] = 0
    phi[i,jmax] = 1

#defining Fm on each face of the cell:
for i in range(2, imax,2):
    for j in range(2,jmax,2):
        Fm[i-1,j] = rho*u*del_y
        Fm[i+1,j] = rho*u*del_y
        Fm[i,j-1] = rho*v*del_x
        Fm[i,j+1] = rho*v*del_x


scheme = int(input("Which convection scheme you want?\n0)FOU 1)CD\n2)SOU\n3)QUICK\n>"))
n=0
time = 0
sum = 0
rmsphi = 0.005
while rmsphi > 0.00001:
    n = n+1
    if scheme == 0:
        for i in range(2, imax, 2):
            for j in range(2, jmax, 2):
                phi[i+1,j] = phi[i,j]
                #phi[i-1,j] = phi[i-2,j]
                phi[i,j+1] = phi[i,j]
                #phi[i,j-1] = phi[i,j-2]

    elif scheme == 1:
        for i in range(2,imax+1,2):
            for j in range(2, jmax+1,2):
                phi_old[i,j] = phi[i,j]
                #convection scheme of central difference used.(CD)
                #if ((i>2 and i<imax) and (j>2 and j<jmax)):
                phi[i-1,j] = (phi[i-2,j]+ phi[i,j])/2 #West face
                phi[i,j-1] = (phi[i,j-2]+ phi[i,j])/2 #South face
                phi[i+1,j] = (phi[i+2,j]+ phi[i,j])/2 #East face
                    #phi[i,j] = phi_old[i,j] - (del_t/(rho*del_v))*C[i,j]
                phi[i,j+1] = (phi[i,j+2]+ phi[i,j])/2 #North face

    elif scheme ==2:
        for j in range(4, jmax,2):
            i = 4
            phi_old[i,j] = phi[i,j]
            phi[i-1,j] = (2*phi[i-2,j]- phi[i-3,j]) #East face of left border cells
        for i in range(4, imax,2):
            j = 4
            phi_old[i,j] = phi[i,j]
            phi[i,j-1] = (2*phi[i,j-2]- phi[i,j-3]) #North face of bottom border cells
        for i in range(6,imax,2):
            for j in range(6, jmax,2):
                phi_old[i,j] = phi[i,j]
                #Convection scheme of Second order upwind (SOU)
                phi[i-1,j] = (3*phi[i-2,j]- phi[i-4,j])/2 #West face
                #phi[i+1,j] = (3*phi[i,j]- phi[i-2,j])/2 #East face
                phi[i,j-1] = (3*phi[i,j-2]- phi[i,j-4])/2 #South face
                #phi[i,j+1] = (3*phi[i,j]- phi[i,j-2])/2 #North face
    elif scheme ==3:
        for j in range(4, jmax,2):
            i = 4
            phi_old[i,j] = phi[i,j]
            phi[i-1,j] = (1/3*phi[i,j]+phi[i-2,j]-1/3*phi[i-3,j]) #East face of left border cells
        for i in range(4, imax,2):
            j = 4
            phi_old[i,j] = phi[i,j]
            phi[i,j-1] = (1/3*phi[i,j]+phi[i,j-2]-1/3*phi[i,j-3]) #North face of bottom border cells
        for i in range(6,imax,2):
            for j in range(6, jmax,2):
                phi_old[i,j] = phi[i,j]
                #convection scheme of QUICK
                phi[i-1,j] = (3*phi[i,j]+6*phi[i-2,j]-phi[i-4,j])/8 #West face
                #phi[i+1,j] = (3*phi[i+2,j]+6*phi[i,j]-phi[i-2,j])/8 #East face
                phi[i,j-1] = (3*phi[i,j]+6*phi[i,j-2]-phi[i,j-4])/8 #South face
                #phi[i,j+1] = (3*phi[i+2,j]+6*phi[i,j]-phi[i,j-2])/8 #North face
    #Calculating C[i,j]s for each cell
    for i in range(2,imax,2):
        for j in range(2, jmax,2):
            C[i,j] =  (Fm[i+1,j]*phi[i+1,j])+(Fm[i,j+1]*phi[i,j+1])-(Fm[i-1,j]*phi[i-1,j])-(Fm[i,j-1]*phi[i,j-1])
            phi[i,j] = phi_old[i,j] - (del_t/(rho*del_v))*C[i,j]
            #rms error calculations
            sum = sum + (phi[i,j]-phi_old[i,j])**2
    #rmsphi = np.sqrt(sum/((imax-1)*(jmax-1)/4))
    rmsphi = np.sqrt(sum/(cells))
    time = time + del_t
    print(f"Calculating phi after {time}secs.")
    print(f"rms error for {n}th iteration is {rmsphi}")
    if rmsphi > 0.1:
        if scheme == 0:
            target = open("covection_FOU.txt",'w')
            target.write("X, Y, phi\n")
        elif scheme ==1:
            target = open("covection_CD.txt",'w')
            target.write("X, Y, phi\n")
        elif scheme ==2:
            target = open("convection_SOU.txt", 'w')
            target.write("X, Y, phi\n")
        elif scheme ==3:
            target = open("convetion_QUICK.txt", 'w')
            target.write("VARIABLES = X Y phi\n")

        for j in range(1,jmax):
            target.write(f"0 {(j-1)*del_y/2} {phi[1,j]}\n")
            target.write(f"{W} {(j-1)*del_y/2} {phi[imax,j]}\n")
        for i in range(1,imax):
            target.write(f"{(i-1)*del_y/2} 0 {phi[i,1]}\n")
            target.write(f"{(i-1)*del_y/2} {W} {phi[i,jmax]}\n")

        for i in range(2,imax+1,2):
            for j in range(2,jmax+1,2):
                target.write(f"{(i-1)*del_x/2} {(j-1)*del_y/2} {phi[i,j]}\n")
