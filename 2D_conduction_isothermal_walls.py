import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt


from numpy import *
from math import pi
import scipy.linalg
from datetime import datetime
import xlsxwriter

T_left = 100
T_right = 300
T_top = 400
T_bottom = 200
phi = np.empty((8,8))

a_corner = [(2,2), (2,6), (6,2), (6,6)]

a_boundary = list()
for j in range(1,8):
    a_boundary.append((1,j))
for k in range(1,8):
    a_boundary.append((k,1))
for j in range(1,8):
    a_boundary.append((7,j))
for k in range(1,8):
    a_boundary.append((k,7))
print(a_boundary)

a_border = list()
for j in range(3,6):
    a_border.append((2,j))
for k in range(3,6):
    a_border.append((k,2))
for j in range(3,6):
    a_border.append((6,j))
for k in range(3,6):
    a_border.append((k,6))
print(a_border)

a_interior = list()
for i in range(3,6):
    for j in range(3,6):
        a_interior.append((i,j))
print(a_interior)

for i in range(1,8):
    for j in range(1,8):
        phi[i,j] = 0

for j in range(1,8):
    phi[1,j] = T_left
    phi[7,j] = T_right

for i in range(1,8):
    phi[i,1] = T_bottom
    phi[i,7] = T_top

for p in range(2,7):
    for q in range(2,7):
        phi[p,q] = k

k = float(input("What is your initial guess: "))

def phi_corner(i,j):
    if i == 2 and j == 2:
        phi[i,j] = (phi[i+1,j]+ 2*phi[i-1,j]+phi[i,j+1]+2*phi[i,j-1])/6
    elif i == 2 and j == 6:
        phi[i,j] = (phi[i+1,j]+ 2*phi[i-1,j]+2*phi[i,j+1]+ phi[i,j-1])/6
    elif i== 6 and j == 2:
        phi[i,j] = (2*phi[i+1,j]+ phi[i-1,j]+phi[i,j+1]+2*phi[i,j-1])/6
    elif i == 6 and j == 6:
        phi[i,j] = (2*phi[i+1,j]+ phi[i-1,j]+2*phi[i,j+1]+phi[i,j-1])/6
    return(phi[i,j])

def phi_border(i,j):
    if i == 2:
        phi[i,j] = (phi[i+1,j]+ 2*phi[i-1,j]+phi[i,j+1]+phi[i,j-1])/5
    elif j == 6:
        phi[i,j] = (phi[i+1,j]+ phi[i-1,j]+2*phi[i,j+1]+ phi[i,j-1])/5
    elif i == 6:
        phi[i,j] = (2*phi[i+1,j]+ phi[i-1,j]+phi[i,j+1]+phi[i,j-1])/5
    elif j == 2:
        phi[i,j] = (phi[i+1,j]+ phi[i-1,j]+phi[i,j+1]+2*phi[i,j-1])/5
    return(phi[i,j])
    print(phi[i,j])

def phi_interior(i,j):
    phi[i,j] = (phi[i+1,j]+ phi[i-1,j]+phi[i,j+1]+phi[i,j-1])/4
    return(phi[i,j])


workbook = xlsxwriter.Workbook('Temperature_profile.xlsx')
worksheet = workbook.add_worksheet()
bold = workbook.add_format({'bold': 3})

merge_format = workbook.add_format({
    'bold': 1,
    'align': 'center',
    'valign': 'vcenter'})


row_counter = 3
column_counter = 3

for n in range(1,10):
    worksheet.merge_range(4 + 10*(n-1), 3, 3+n*10, 3, f'Iteration({n})', merge_format)

    for l,m in a_corner:
        phi_corner_new = phi_corner(l,m)
        phi_corner_new = int(phi_corner_new)
        print(f"The {n}th iteration values of phi_corner({l},{m}):", phi_corner_new)
        worksheet.write((row_counter + m), (column_counter + l), phi_corner_new)

    for r,s in a_border:
        phi_border_new = phi_border(r,s)
        phi_border_new = int(phi_border_new)
        print(f"The {n}th iteration values of phi_border({r},{s}):", phi_border_new)
        worksheet.write((row_counter + s), (column_counter + r), phi_border_new)

    for x,y in a_interior:
        phi_interior_new = phi_interior(x,y)
        phi_interior_new = int(phi_interior_new)
        print(f"The {n}th iteration values of phi_interior({x},{y}):", phi_interior_new)
        worksheet.write((row_counter + y), (column_counter + x), phi_interior_new)
        for c,d in a_boundary:
            phi_boundary_new = phi[c,d]
            worksheet.write((row_counter + d), (column_counter + c), phi_boundary_new)
    row_counter += 10

workbook.close()
