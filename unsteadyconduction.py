import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

from numpy import *
from math import pi
import scipy.linalg
from datetime import datetime
import xlsxwriter

#Material Properties
km = float(input("What is the thermal conductivity of material?(W/mK) \n >"))
rho = float(input("What is the density of the matierial?(kg/m^3) \n >"))
Cp = float(input("What is the specific heat of the material?(W/m^2K) \n >"))
qg = float(input("Enter the internal heat generation rate(W/m3) \n >"))

#Geometrical properties
#T_0 = input("What is the reference temperature in deg Cel.? \n >")
#T_0 = float(T_0)

L = float(input("What is the length of plate? \n >"))
W = float(input("What is the width of the plate? \n >"))
t = float(input("What is the thickness of the plate? \n >"))

#Inputs required for grid generation
print("For better results select cell dimensions so that you will get around 4000 cells.")

del_x = float(input("What is cell length? \n >"))
del_y= float(input("What is cell width? \n >"))
cells = int((W*L)/(del_x*del_y))
print(f"You have {cells} number of cells in grid.")
alp = km/(rho*Cp)
del_t = 0.4/(alp*(((1/del_x)**2)+((1/del_y)**2)))



# Define coordinates for data points

imax = int(L/del_x + 2)
jmax = int(W/del_y + 2)
#Generate a grid object using numpy
phi = np.empty((imax+1,jmax+1))
phi_old = np.empty((imax+1, jmax+1))
qw = np.empty((imax+1, jmax+1))
#for i in range(2,imax-1):
    #for j in range(2,jmax-1):
        #points.append(((i - 1.5)*del_x,( j - 1.5)*del_y))

a_corner = [(2,2), (2,jmax-1), (imax-1,2), (imax-1,jmax-1)]
corner_points = [(1,1), (1,jmax), (imax,1), (imax,jmax)]

left_bound = list()
for j in range(1,jmax+1):
    left_bound.append((1,j))
    #points.append((0,(j-1.5)*del_y))

right_bound = list()
for j in range(1,jmax+1):
    right_bound.append((imax,j))
    #points.append((L,(j-1.5)*del_y))

top_bound = list()
for f in range(2,imax):
    top_bound.append((f,jmax))
    #points.append(((k-1.5)*del_x,W))

bottom_bound = list()
for f in range(2,imax):
    bottom_bound.append((f,1))
    #points.append(((k-1.5)*del_x,0))

a_boundary = (left_bound + right_bound + top_bound + bottom_bound)

a_border = list()
for j in range(3,jmax-1):
    a_border.append((2,j))
for k in range(3,imax-1):
    a_border.append((k,2))
for j in range(3,jmax-1):
    a_border.append((imax-1,j))
for k in range(3,imax-1):
    a_border.append((k,jmax-1))

a_interior = list()
for i in range(3,imax-1):
    for j in range(3,jmax-1):
        a_interior.append((i,j))

#Setting all phi elements equal to zero
for i in range(1,imax+1):
    for j in range(1,jmax+1):
        phi[i,j] = 0
        phi_old[i,j] = 0



#Boundary conditions
Bound_conds = ['Isothermal(Dirichlett)', 'Heat flux(Neumann)', 'Mixed(Robin)']

T_w_left = 0
T_w_right = 0
T_w_top = 0
T_w_bottom = 0
h_conv_left = 0
h_conv_right = 0
h_conv_top = 0
h_conv_bottom = 0
qw_left = 0
qw_right = 0
qw_top = 0
qw_bottom = 0

#Function defined for constant wall temperature.
def Isothermal(i,j):
    #T_w = float(input("What is the wall temperature? \n >"))
    if i == 1:
        phi[i,j] = T_w_left
    elif i == imax:
        phi[i,j] = T_w_right
    elif j == jmax:
        phi[i,j] = T_w_top
    elif j == 1:
        phi[i,j] = T_w_bottom
    return(phi[i,j])



#Function defined for constant wall flux.
def Heat_flux(i,j):
    #qw = float(input("What is heat flux through the wall? (W/m^2)\n >"))
    if i == 1:
        phi[i,j] = phi[i+1,j] - ((qw_left*del_x)/(2*km))
    elif i == imax:
        phi[i,j] = phi[i-1,j] - ((qw_right*del_x)/(2*km))
    elif j == jmax:
        phi[i,j] = phi[i,j-1] - ((qw_top*del_y)/(2*km))
    elif j == 1:
        phi[i,j] = phi[i,j+1] - ((qw_bottom*del_y)/(2*km))
    return(phi[i,j])


#Function defined for convective heat transfer from the wall.
def Convective(i,j):
    #h_conv = float(input("What is heat convection over the wall? (W/m^2K) \n >"))
    #T_out = float(input("What is the temperature of the convecting fluid? (deg Cel.) \n >"))
    phi_out = T_out
    if i == 1:
        k1 = h_conv_left*del_x/(2*km)
        phi[i,j] = (phi[i+1,j] + k1*phi_out)/(1+ k1)
    elif i == imax:
        k2 = h_conv_right*del_x/(2*km)
        phi[i,j] = (phi[i-1,j] + k2*phi_out)/(1+ k2)
    elif j == jmax:
        k3 = h_conv_top*del_y/(2*km)
        phi[i,j] = (phi[i,j-1] + k3*phi_out)/(1+ k3)
    elif j == 1:
        k4 = h_conv_bottom*del_y/(2*km)
        phi[i,j] = (phi[i, j+1] + k4*phi_out)/(1+ k4)
    return(phi[i,j])



#Using loop to define boundary conditions.
for k in range(1,5):
    if k == 1:
        bound_cond_left = int(input("What is the left boundary condtion? \n 1)Isothermal(Dirichlett) \n 2)Constant Heat flux(Neumann) \n 3)Convective(Robin) \n >"))
        if bound_cond_left == 1:
            T_w_left = float(input("What is the wall temperature? \n >"))
            for p,q in left_bound:
                phi[p,q] = Isothermal(p,q)
        elif bound_cond_left == 2:
            qw_left = float(input("What is heat flux through the wall? (W/m^2)\n >"))
            for p,q in left_bound:
                phi[p,q] = Heat_flux(p,q)
        elif bound_cond_left == 3:
            h_conv_left = float(input("What is heat convection over the wall? (W/m^2K) \n >"))
            T_out = float(input("What is the temperature of the convecting fluid? (deg Cel.) \n >"))
            for p,q in left_bound:
                phi[p,q] = Convective(p,q)

    elif k == 2:
        bound_cond_right = int(input("What is the right boundary condtion? \n 1)Isothermal(Dirichlett) \n 2)Constant Heat flux(Neumann) \n 3)Convective(Robin) \n >"))
        if bound_cond_right == 1:
            T_w_right = float(input("What is the wall temperature? \n >"))
            for p,q in right_bound:
                phi[p,q] = Isothermal(p,q)
        elif bound_cond_right == 2:
            qw_right = float(input("What is heat flux through the wall? (W/m^2)\n >"))
            for p,q in right_bound:
                phi[p,q] = Heat_flux(p,q)
        elif bound_cond_right == 3:
            h_conv_right = float(input("What is heat convection over the wall? (W/m^2K) \n >"))
            T_out = float(input("What is the temperature of the convecting fluid? (deg Cel.) \n >"))
            for p,q in right_bound:
                phi[p,q] = Convective(p,q)


    elif k == 3:
        bound_cond_top = int(input("What is the top boundary condtion? \n 1)Isothermal(Dirichlett) \n 2)Constant Heat flux(Neumann) \n 3)Convective(Robin) \n >"))
        if bound_cond_top == 1:
            T_w_top = float(input("What is the wall temperature? \n >"))
            for p,q in top_bound:
                phi[p,q] = Isothermal(p,q)
        elif bound_cond_top == 2:
            qw_top = float(input("What is heat flux through the wall? (W/m^2)\n >"))
            for p,q in top_bound:
                phi[p,q] = Heat_flux(p,q)
        elif bound_cond_top == 3:
            h_conv_top = float(input("What is heat convection over the wall? (W/m^2K) \n >"))
            T_out = float(input("What is the temperature of the convecting fluid? (deg Cel.) \n >"))
            for p,q in top_bound:
                phi[p,q] = Convective(p,q)

    elif k == 4:
        bound_cond_bottom = int(input("What is the bottom boundary condtion? \n 1)Isothermal(Dirichlett) \n 2)Constant Heat flux(Neumann) \n 3)Convective(Robin) \n >"))
        if bound_cond_bottom == 1:
            T_w_bottom = float(input("What is the wall temperature? \n >"))
            for p,q in bottom_bound:
                phi[p,q] = Isothermal(p,q)
        elif bound_cond_bottom == 2:
            qw_bottom = float(input("What is heat flux through the wall? (W/m^2)\n >"))
            for p,q in bottom_bound:
                phi[p,q] = Heat_flux(p,q)
        elif bound_cond_bottom == 3:
            h_conv_bottom = float(input("What is heat convection over the wall? (W/m^2K) \n >"))
            T_out = float(input("What is the temperature of the convecting fluid? (deg Cel.) \n >"))
            for p,q in bottom_bound:
                phi[p,q] = Convective(p,q)





guess = float(input("What is your initial condition for interior cells? (deg Cel.) \n >"))
for i in range(2,imax):
    for j in range(2,jmax):
        phi[i,j] = guess

def phi_corner(i,j):
    #points.append(((i - 1.5)*del_x,( j - 1.5)*del_y))
    if i == 2 and j == 2:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phi_old[i-1,j])/(del_x**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    elif i == 2 and j == jmax-1:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((phi_old[i+1,j]-3*phi_old[i,j]+2*phi_old[i-1,j])/(del_x**2) + (2*phi_old[i,j+1]-3*phi_old[i,j]+phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    elif i== imax-1 and j == 2:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((2*phi_old[i+1,j]-3*phi_old[i,j]+phi_old[i-1,j])/(del_x**2) + (phi_old[i,j+1]-3*phi_old[i,j]+2*phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    elif i == imax-1 and j == jmax-1:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((2*phi_old[i+1,j]-3*phi_old[i,j]+phi_old[i-1,j])/(del_x**2) + (2*phi_old[i,j+1]-3*phi_old[i,j]+phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    return(phi[i,j])

def phi_border(i,j):
    #points.append(((i - 1.5)*del_x,( j - 1.5)*del_y))
    if i == 2:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((phi_old[i+1,j]-3*phi_old[i,j]+ 2*phi_old[i-1,j])/(del_x**2) + (phi_old[i,j+1]-2*phi_old[i,j]+ phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    elif j == jmax-1:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((phi_old[i+1,j]-2*phi_old[i,j]+phi_old[i-1,j])/(del_x**2) + (2*phi_old[i,j+1]-3*phi_old[i,j]+phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    elif i == imax-1:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((2*phi_old[i+1,j]-3*phi_old[i,j]+phi_old[i-1,j])/(del_x**2) + (phi_old[i,j+1]-2*phi_old[i,j]+phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    elif j == 2:
        phi[i,j] = phi_old[i,j] + del_t*(alp*((phi_old[i+1,j]-2*phi_old[i,j]+ phi_old[i-1,j])/(del_x**2) + (phi_old[i,j+1]-3*phi_old[i,j]+ 2*phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    return(phi[i,j])

def phi_interior(i,j):
    #points.append(((i - 1.5)*del_x,( j - 1.5)*del_y))
    phi[i,j] = phi_old[i,j] + del_t*(alp*((phi_old[i+1,j]-2*phi_old[i,j]+phi_old[i-1,j])/(del_x**2) + (phi_old[i,j+1]-2*phi_old[i,j]+phi_old[i,j-1])/(del_y**2)) + qg/(rho*Cp))
    return(phi[i,j])


workbook = xlsxwriter.Workbook('Temperature_profile_grid.xlsx')
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': 3})

merge_format = workbook.add_format({
    'bold': 1,
    'align': 'center',
    'valign': 'vcenter'})

row_counter = 3
column_counter = 3
#iter = int(input("How many iterations do you want?"))
n = 0
time = 0
rmsphi = 0.5
while rmsphi/del_t > 0.0001:
    n = n+1
    time = time + del_t
    #worksheet.merge_range(4 + (jmax+3)*(n-1), 3, 3+n*(jmax+3), 3, f'Iteration({n})', merge_format)
    print(f"Calculating temperatures for {n}th iteration and after time {time} sec..")

    for i in range(1,5):
        if i == 1:
            if bound_cond_left == 1:
                for p,q in left_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Isothermal(p,q)
            elif bound_cond_left == 2:
                for p,q in left_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Heat_flux(p,q)
            elif bound_cond_left == 3:
                for p,q in left_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Convective(p,q)

        elif i == 2:
            if bound_cond_right == 1:
                for p,q in right_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Isothermal(p,q)
            elif bound_cond_right == 2:
                for p,q in right_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Heat_flux(p,q)
            elif bound_cond_right == 3:
                for p,q in right_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Convective(p,q)


        elif i == 3:
            if bound_cond_top == 1:
                for p,q in top_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Isothermal(p,q)
            elif bound_cond_top == 2:
                for p,q in top_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Heat_flux(p,q)
            elif bound_cond_top == 3:
                for p,q in top_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Convective(p,q)

        elif i == 4:
            if bound_cond_bottom == 1:
                for p,q in bottom_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Isothermal(p,q)
            elif bound_cond_bottom == 2:
                for p,q in bottom_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Heat_flux(p,q)
            elif bound_cond_bottom == 3:
                for p,q in bottom_bound:
                    phi_old[p,q] = phi[p,q]
                    phi[p,q] = Convective(p,q)


    for l,m in a_corner:
        phi_old[l,m] = phi[l,m]
        phi_corner_new = phi_corner(l,m)
        #worksheet.write((row_counter + m), (column_counter + l), phi_corner_new )

    for r,s in a_border:
        phi_old[r,s] = phi[r,s]
        phi_border_new = phi_border(r,s)
        #worksheet.write((row_counter + s), (column_counter + r), phi_border_new )

    for x,y in a_interior:
        phi_old[x,y] = phi[x,y]
        phi_interior_new = phi_interior(x,y)
        #worksheet.write((row_counter + y), (column_counter + x), phi_interior_new)

    diff = 0
    for i in range(2, imax):
        for j in range(2,jmax):
            diff = diff + (phi[i,j] - phi_old[i,j])**2
    rmsphi = np.sqrt(diff/((imax-2)*(jmax-2)))
    print(f"RMS Error for {n}th iteration is :{rmsphi}")

    if n%100 == 0:
        worksheet.merge_range(4 , 3, 3+(jmax+3), 3, f'Iteration({n})', merge_format)
        target = open(f"temperature_{time}sec.dat", "w")
        #target1 = open("Heat_flux_top.txt", "w")
        #target2 = open("Heat_flux_bottom.txt", "w")
        target.write("VARIABLES = X, Y, T \n")
        #target1.write("VARIABLES = X, qt\n")
        #target2.write("VARIABLES = X, qb\n")

        for c,d in a_boundary:
            phi_boundary_new = float(phi[c,d])
            worksheet.write((row_counter + d), (column_counter + c), phi_boundary_new)
            if c==1 and (d !=1 and d != jmax):
                qw_boundary_new = float(km*(phi[c,d]-phi[c+1,d])/(del_x/2))
                target.write(f"{0} {(d-1.5)*del_y} {phi_boundary_new}\n")
                #target2.write(f"{0} {(d-1.5)*del_y} {qw_boundary_new}\n")
            elif c==imax and (d !=1 and d != jmax):
                qw_boundary_new = float(km*(phi[c,d]-phi[c-1,d])/(del_x/2))
                target.write(f"{L} {(d-1.5)*del_y} {phi_boundary_new}\n")
                #target2.write(f"{L} {(d-1.5)*del_y} {qw_boundary_new}\n")
            elif d==1 and (c !=1 and c != imax):
                qw_boundary_new = float(km*(phi[c,d]-phi[c,d+1])/(del_x/2))
                target.write(f"{(c-1.5)*del_x} {0} {phi_boundary_new}\n")
                #target1.write(f"{(c-1.5)*del_x} {qw_boundary_new}\n")
            elif d==jmax and (c !=1 and c != imax):
                qw_boundary_new = float(km*(phi[c,d]-phi[c,d-1])/(del_x/2))
                target.write(f"{(c-1.5)*del_x} {W} {phi_boundary_new}\n")
                #target2.write(f"{(c-1.5)*del_x} {qw_boundary_new}\n")

        for i,j in corner_points:
            if i==1 and j==1:
                phi_boundary_new = float(phi[i,j])
                qw_boundary_new = float(km*(phi[i,j]-phi[i,j+1])/(del_y/2))
                worksheet.write((row_counter + j), (column_counter + i), phi_boundary_new)
                target.write(f"{0} {0} {phi_boundary_new}\n")
                #target2.write(f"{0} {qw_boundary_new}\n")
            elif i==1 and j==jmax:
                phi_boundary_new = float(phi[i,j])
                qw_boundary_new = float(km*(phi[i,j-1]-phi[i,j])/(del_y/2))
                worksheet.write((row_counter + j), (column_counter + i), phi_boundary_new)
                target.write(f"{0} {W} {phi_boundary_new}\n")
                #target1.write(f"{0} {qw_boundary_new}\n")
            elif i==imax and j==1:
                phi_boundary_new = float(phi[i,j])
                qw_boundary_new = float(km*(phi[i,j]-phi[i,j+1])/(del_y/2))
                worksheet.write((row_counter + j), (column_counter + i), phi_boundary_new)
                target.write(f"{0} {L} {phi_boundary_new}\n")
                #target2.write(f"{L} {qw_boundary_new}\n")
            if i==imax and j==jmax:
                phi_boundary_new = float(phi[i,j])
                qw_boundary_new = float(km*(phi[i,j-1]-phi[i,j])/(del_y/2))
                worksheet.write((row_counter + j), (column_counter + i), phi_boundary_new)
                target.write(f"{L} {W} {phi_boundary_new}\n")
                #target1.write(f"{L} {qw_boundary_new}\n")

        for l,m in a_corner:
            phi_corner_new = phi_corner(l,m)
            phi_corner_new = float(phi_corner_new)
            worksheet.write((row_counter + m), (column_counter + l), phi_corner_new)
            target.write(f"{(l - 1.5)*del_x} {( m - 1.5)*del_y} {phi_corner_new}\n")

        for r,s in a_border:
            phi_border_new = phi_border(r,s)
            phi_border_new = float(phi_border_new)
            worksheet.write((row_counter + s), (column_counter + r), phi_border_new)
            target.write(f"{(r - 1.5)*del_x} {( s - 1.5)*del_y} {phi_border_new}\n")

        for x,y in a_interior:
            phi_interior_new = phi_interior(x,y)
            phi_interior_new = float(phi_interior_new)
            worksheet.write((row_counter + y), (column_counter + x), phi_interior_new)
            target.write(f"{(x - 1.5)*del_x} {( y - 1.5)*del_y} {phi_interior_new}\n")


    #row_counter += (jmax+3)


workbook.close()
