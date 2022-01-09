#Inviscid Flow code
import numpy as np
from numpy import *

L, W, L1, W1 = 10, 10, 2, 2  #Lengths are in m
Q_in = 10 #m^2/sec
cells = 2500
del_x = np.sqrt((W*L)/cells)
del_y = del_x

#domain extreames
imax = int(2*L/del_x +1)
jmax = int(2*W/del_y +1)

#internal object grid point extreames
i_left = int((L-L1)/del_x +1)
i_right = int(i_left + 2*L1/del_x)
j_bottom = int((W-W1)/del_y+1)
j_top = int(j_bottom + 2*W1/del_y)
x_cord = np.empty((imax+1))
y_cord = np.empty((jmax+1))
phi = np.empty((imax+1, jmax+1))  #2d_array for variable values
phi_old = np.empty((imax+1, jmax+1)) #2d_array for old variable values

x_cord[:], y_cord[:], phi[:,:], phi_old[:,:] = 0,0,0,0

guess = float(input("What is your initial guess for interior cells? (deg Cel.) \n >"))
for i in range(2,imax,2):
    for j in range(2,jmax,2):
        phi[i,j] = guess

for j in range(1,jmax+1):
    u_in = Q_in/L
    phi[1,j] = phi[1,j-1] + u_in*del_y/2
    phi[imax,j] = phi[imax,j-1] + u_in*del_y/2

for i in range(1,imax+1):
    phi[i,1] = 0
    phi[i,jmax] = 10

for i in range(i_left, i_right+1):
    phi[i,j_bottom] = 5
    phi[i,j_top] = 5
for j in range(j_bottom, j_top+1):
    phi[i_left,j] = 5
    phi[i_right,j] = 5

interior_cells = list()
for i in range(4, imax-2, 2):
    for j in range(4, jmax-2, 2):
        interior_cells.append((i,j))
body_cells = list()
for i in range(i_left+1, i_right,2):
    for j in range(j_bottom+1, j_top,2):
        body_cells.append((i,j))

k1 = (W*L)/(del_x*del_y) - (W1*L1)/(del_x*del_y)
# Actual numbers of cells in the control volume.

n=0
rmsphi = .5
diff = 0
while rmsphi > 0.000001:
    n = n+1
    print(f"Calculating temperatures for {n}th iteration..")
    for i in range(1, imax+1):
        for j in range(1, jmax+1):
            phi_old[i,j] = phi[i,j]
    i, j = 2,2
    phi[i,j] = (phi[i+2,j]+ 2*phi[i-1,j]+phi[i,j+2]+2*phi[i,j-1])/6
    i,j = 2, jmax-1
    phi[i,j] = (phi[i+2,j]+ 2*phi[i-1,j]+2*phi[i,j+1]+ phi[i,j-2])/6
    i,j = imax-1, 2
    phi[i,j] = (2*phi[i+1,j]+ phi[i-2,j]+phi[i,j+2]+2*phi[i,j-1])/6
    i ,j = imax-1,jmax-1
    phi[i,j] = (2*phi[i+1,j]+ phi[i-2,j]+2*phi[i,j+1]+phi[i,j-2])/6

    for i in range(4, imax-1, 2):
        j = jmax-1
        phi[i,j] = (phi[i+2,j]+ phi[i-2,j]+2*phi[i,j+1]+ phi[i,j-2])/5
        j = 2
        phi[i,j] = (phi[i+2,j]+ phi[i-2,j]+phi[i,j+2]+2*phi[i,j-1])/5
    for j in range(4, jmax-1, 2):
        i = 2
        phi[i,j] = (phi[i+2,j]+ 2*phi[i-1,j]+phi[i,j+2]+phi[i,j-2])/5
        i = imax-1
        phi[i,j] = (2*phi[i+1,j]+ phi[i-2,j]+phi[i,j+2]+phi[i,j-2])/5

    for i in range(i_left+1, i_right+1,2):
        j = j_bottom -1
        phi[i,j] = (phi[i+2,j]+ phi[i-2,j]+2*phi[i,j+1]+ phi[i,j-2])/5
        j = j_top + 1
        phi[i,j] = (phi[i+2,j]+ phi[i-2,j]+phi[i,j+2]+2*phi[i,j-1])/5

    for j in range(j_bottom+1, j_top+1,2):
        i = i_left-1
        phi[i,j] = (2*phi[i+1,j]+ phi[i-2,j]+phi[i,j+2]+phi[i,j-2])/5
        i = i_right+1
        phi[i,j] = (phi[i+2,j]+ 2*phi[i-1,j]+phi[i,j+2]+phi[i,j-2])/5


    for i,j in interior_cells:
        if (i,j) not in body_cells:
            phi[i,j] = (phi[i+2,j]+ phi[i-2,j]+phi[i,j+2]+phi[i,j-2])/4
    diff = 0
    for i in range(2, imax):
        for j in range(2, jmax):
            if i not in range(i_left+1, i_right,2):
                if j not in range(j_bottom+1, j_top,2):
                    phi[i,j] = float(phi[i,j])
                    phi_old[i,j] = float(phi_old[i,j])
                    diff = diff + (phi[i,j] - phi_old[i,j])**2
    rmsphi = np.sqrt(diff/k1)
    print(f"RMS Error for {n}th iteration is :{rmsphi}")
    #if n > 1000:
    if rmsphi < 0.000001:
        target_phi = open("inviscid_data.dat", "w")
        target_phi.write("VARIABLES = X Y, phi\nZONE T = sample\n")

        for j in range(2,jmax+1):
            x_cord[1] = 0
            y_cord[j] = y_cord[j-1] + del_y/2
            target_phi.write(f"{0} {y_cord[j]} {phi[1,j]}\n")
            x_cord[imax] = L
            target_phi.write(f"{L} {y_cord[j]} {phi[imax,j]}\n")

        for i in range(2,imax+1):
            x_cord[i] = x_cord[i-1] + del_x/2
            y_cord[jmax] = W
            target_phi.write(f"{x_cord[i]} {W} {phi[i,jmax]}\n")
            y_cord[1] = 0
            target_phi.write(f"{x_cord[i]} {0} {phi[i,1]}\n")

        for i in range(2, imax, 2):
            for j in range(2, jmax, 2):
                target_phi.write(f"{x_cord[i]}  {y_cord[j]} {phi[i,j]}\n")

        quit()
