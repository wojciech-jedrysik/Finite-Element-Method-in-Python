from __future__ import division

import Structures as S
import Matrices as MS
import GlobalMatrices as GMS
import FileReader as FR
import numpy as np


def calc_equation():
    t0 = []
    for i in range(global_data.nN):
        t0.append(0)
    for i in range(len(grid.elements)):
        for j in range(len(grid.elements[i].nodes)):
            t0[grid.elements[i].nodes[j].id - 1] = grid.elements[i].nodes[j].temp
    new_vector = list(np.array(global_C.gC).dot(np.array(t0)))
    for i in range(global_data.nN):
        new_vector[i] += global_P.gP[i]
    t1_vector = list(np.linalg.solve(np.array(global_H.gH), np.array(new_vector)))
    min_temp = min(t1_vector)
    max_temp = max(t1_vector)
    return t1_vector, min_temp, max_temp, new_vector


global_data = S.GlobalData()
init_temp, sim_time, step_time, amb_temp, alfa, Height, Width, N_Height, N_Width, spec_heat, conduct, dens = FR.read('GlobalData.txt')
global_data.assign_values(init_temp, sim_time, step_time, amb_temp, alfa, Height, Width, N_Height, N_Width, spec_heat, conduct,
                          dens)
grid = S.FEMgrid()
grid.create(global_data)
grid.show()

for i in range(len(grid.elements)):
    points_tab = []
    are_edges = []
    for j in range(len(grid.elements[i].nodes)):
        point = [grid.elements[i].nodes[j].x, grid.elements[i].nodes[j].y]
        points_tab.append(point)
    for j in range(4):
        if grid.elements[i].nodes[j].is_out is True and grid.elements[i].nodes[(j+1)%4].is_out is True:
            are_edges.append(1)
        else:
            are_edges.append(0)
    temp = MS.matrixH()
    temp.calc_H(points_tab, global_data.conduct)
    grid.elements[i].local_H = temp.H
    temp = MS.matrixH_bc()
    temp.calc_H_bc(points_tab, global_data.alfa, are_edges)
    grid.elements[i].local_H_bc = temp.H_bc
    grid.elements[i].add_bc_to_H()
    temp = MS.matrixC()
    temp.calc_C(points_tab, global_data.spec_heat, global_data.dens)
    grid.elements[i].local_C = temp.C
    temp = MS.vectorP()
    temp.calc_P(grid.elements[i].nodes, global_data.amb_temp, global_data.alfa)
    grid.elements[i].local_P = temp.P

global_H = GMS.Global_H()
global_H.create(grid.elements, global_data.nN)
global_C = GMS.Global_C()
global_C.create(grid.elements, global_data.nN)
print("Matrix [C]\n")
for i in global_C.gC:
    print(i)
print()
global_P = GMS.Global_P()
global_P.create(grid.elements, global_data.nN)
for i in range(len(global_C.gC)):
    for j in range(len(global_C.gC)):
        global_C.gC[i][j] /= global_data.step_time
for i in range(len(global_H.gH)):
    for j in range(len(global_H.gH)):
        global_H.gH[i][j] += global_C.gC[i][j]

print("########### SIMULATION ###########\n")
iterations_number = int(global_data.sim_time / global_data.step_time)
print("Minimal temperature:", global_data.init_temp)
print("Maximal temperature:", global_data.init_temp, "\n")
print("\nMatrix [H] = [H] + [C]/dT\n")
for i in global_H.gH:
    print(i)
print("\nMatrix [C] = [C]/dT\n")
for i in global_C.gC:
    print(i)
print("\n{P}\n")
print(global_P.gP, "\n")
for i in range(iterations_number):
    print("Iteration", i)
    temps, min_temp, max_temp, p = calc_equation()
    print("\n{P} = {P} + {[C]/dT}*{T0}")
    print(p, "\n")
    print("Minimal temperature:", round(min_temp, 3))
    print("Maximal temperature:", round(max_temp, 3), "\n")
    for j in range(len(grid.elements)):
        for k in range(4):
            grid.elements[j].nodes[k].temp = temps[grid.elements[j].nodes[k].id - 1]
