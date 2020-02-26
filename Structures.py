import math as m


class GlobalData:
    def __init__(self):
        self.H = 0  # Height
        self.W = 0  # Width
        self.nH = 0  # Nodes per height
        self.nW = 0  # Nodes per width
        self.nN = 0  # All nodes
        self.nE = 0  # All elements
        self.init_temp = 0  # Initial temperature
        self.sim_time = 0  # Simulation time
        self.step_time = 0  # Step time
        self.amb_temp = 0  # Ambient temperature
        self.alfa = 0  # Alfa
        self.spec_heat = 0  # Specific heat
        self.conduct = 0  # Conductivity
        self.dens = 0  # Density

    def assign_values(self, init_temp, sim_time, step_time, amb_temp, alfa, H, W, nH, nW, spec_heat, conduct, dens):
        self.H = H
        self.W = W
        self.nH = nH
        self.nW = nW
        self.nN = self.nH * self.nW
        self.nE = (self.nH - 1) * (self.nW - 1)
        self.init_temp = init_temp
        self.sim_time = sim_time
        self.step_time = step_time
        self.amb_temp = amb_temp
        self.alfa = alfa
        self.spec_heat = spec_heat
        self.conduct = conduct
        self.dens = dens


class Node:
    def __init__(self, id=0, x=0, y=0, is_out=False, temp=0):
        self.id = id
        self.x = x
        self.y = y
        self.temp = temp
        self.is_out = is_out


class Element:
    def __init__(self, node_1, node_2, node_3, node_4):
        self.nodes = [node_1, node_2, node_3, node_4]
        self.local_H = []
        self.local_H_bc = []
        self.local_C = []
        self.local_P = []

    def add_bc_to_H(self):
        for i in range(4):
            for j in range(4):
                self.local_H[i][j] += self.local_H_bc[i][j]


class FEMgrid:
    def __init__(self):
        self.nodes = []
        self.elements = []

    def create(self, global_data):
        delta_x = global_data.W / (global_data.nH - 1)
        delta_y = global_data.H / (global_data.nH - 1)
        id = 1
        for i in range(int(global_data.nW)):
            for j in range(int(global_data.nH)):
                x = i * delta_x
                y = j * delta_y
                if x == 0 or y == 0 or x == global_data.W or y == global_data.H:
                    is_out = True
                else:
                    is_out = False
                node = Node(id, x, y, is_out, global_data.init_temp)
                self.nodes.append(node)
                id += 1
        temp = 0
        for i in range(int(global_data.nE)):
            element = Element(self.nodes[temp], self.nodes[temp + int(global_data.nH)], self.nodes[temp + int(global_data.nH + 1)], self.nodes[temp + 1])
            self.elements.append(element)
            temp += 1
            if (temp + 1) % global_data.nH == 0:
                temp += 1

    def show(self):
        print("Nodes:")
        for node in self.nodes:
            print("ID:", node.id)
            print("x =", node.x)
            print("y =", node.y)
            print("Is out?", node.is_out)
            print("Temperature:", node.temp)
            print()
        print()
        print("Elements:")
        nr = 1
        for element in self.elements:
            print("ELEMENT NR", nr)
            print("Node 1:", element.nodes[0].id)
            print("Node 2:", element.nodes[1].id)
            print("Node 3:", element.nodes[2].id)
            print("Node 4:", element.nodes[3].id)
            print()
            nr += 1
        print()


class UniversalElement:
    def __init__(self):
        self.shapeFunctionsValues = []
        self.byXi = []
        self.byEta = []
        self.local_points = [(-1 / m.sqrt(3), -1 / m.sqrt(3)), (1 / m.sqrt(3), -1 / m.sqrt(3)), (1 / m.sqrt(3), 1 / m.sqrt(3)), (-1/m.sqrt(3), 1/m.sqrt(3))]

        for i in range(4):
            N1234 = []
            N1 = ((1 / 4) * (1 - self.local_points[i][0]) * (1 - self.local_points[i][1]))
            N2 = ((1 / 4) * (1 + self.local_points[i][0]) * (1 - self.local_points[i][1]))
            N3 = ((1 / 4) * (1 + self.local_points[i][0]) * (1 + self.local_points[i][1]))
            N4 = ((1 / 4) * (1 - self.local_points[i][0]) * (1 + self.local_points[i][1]))
            N1234.append(N1)
            N1234.append(N2)
            N1234.append(N3)
            N1234.append(N4)
            self.shapeFunctionsValues.append(N1234)
            N1234 = []

            N1 = -0.25 * (1 - self.local_points[i][1])
            N2 = 0.25 * (1 - self.local_points[i][1])
            N3 = 0.25 * (1 + self.local_points[i][1])
            N4 = -0.25 * (1 + self.local_points[i][1])
            N1234.append(N1)
            N1234.append(N2)
            N1234.append(N3)
            N1234.append(N4)
            self.byXi.append(N1234)
            N1234 = []

            N1 = -0.25 * (1 - self.local_points[i][0])
            N2 = -0.25 * (1 + self.local_points[i][0])
            N3 = 0.25 * (1 + self.local_points[i][0])
            N4 = 0.25 * (1 - self.local_points[i][0])
            N1234.append(N1)
            N1234.append(N2)
            N1234.append(N3)
            N1234.append(N4)
            self.byEta.append(N1234)

    def calc_surface_N(self):
        surface_N = []
        v = 1 / m.sqrt(3)
        ksi_eta = [[-v, -1], [1, -v], [v, 1], [-1, v]]
        for i in range(4):
            temp = []
            N1 = 0.25 * (1 - ksi_eta[i][0]) * (1 - ksi_eta[i][1])
            N2 = 0.25 * (1 + ksi_eta[i][0]) * (1 - ksi_eta[i][1])
            N3 = 0.25 * (1 + ksi_eta[i][0]) * (1 + ksi_eta[i][1])
            N4 = 0.25 * (1 - ksi_eta[i][0]) * (1 + ksi_eta[i][1])
            temp.append(N1)
            temp.append(N2)
            temp.append(N3)
            temp.append(N4)
            surface_N.append(temp)
        return surface_N

    def calc_jacobian(self, points_tab, integer_point_index):
        J = []
        first_line = []
        second_line = []
        results = [0, 0, 0, 0]

        for i in range(4):
            # deX/deXi
            results[0] += self.byXi[integer_point_index][i] * points_tab[i][0]
            # deY/deXi
            results[1] += self.byXi[integer_point_index][i] * points_tab[i][1]
            # deX/deEta
            results[2] += self.byEta[integer_point_index][i] * points_tab[i][0]
            # deY/deEta
            results[3] += self.byEta[integer_point_index][i] * points_tab[i][1]

        first_line.append(results[0])
        first_line.append(results[1])
        second_line.append(results[2])
        second_line.append(results[3])

        J.append(first_line)
        J.append(second_line)

        return J

    def calc_det(self, matrix):
        return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0])

    def calc_J_1(self, matrix, det_m):
        J_1_matrix = []
        w1k1 = matrix[1][1] / det_m
        if w1k1 == -0: w1k1 = 0
        w1k2 = -matrix[0][1] / det_m
        if w1k2 == -0: w1k2 = 0
        w2k1 = -matrix[1][0] / det_m
        if w2k1 == -0: w2k1 = 0
        w2k2 = matrix[0][0] / det_m
        if w2k2 == -0: w2k2 = 0
        J_1_matrix.append([w1k1, w1k2])
        J_1_matrix.append([w2k1, w2k2])

        return J_1_matrix

    def calc_deXY(self, all_J_1):
        byX = []
        for i in range(4):
            temp = []
            for j in range(4):
                N_i = (all_J_1[i][0][0] * self.byXi[i][j]) + (all_J_1[i][0][1] * self.byEta[i][j])
                temp.append(N_i)
            byX.append(temp)

        byY = []
        for i in range(4):
            temp = []
            for j in range(4):
                N_i = all_J_1[i][1][0] * self.byXi[i][j] + all_J_1[i][1][1] * self.byEta[i][j]
                temp.append(N_i)
            byY.append(temp)

        return byX, byY


    def multiply_vector(self, vector):
        result = []
        for j in range(4):
            temp = []
            for k in range(4):
                temp.append(vector[j] * vector[k])
            result.append(temp)

        return result

    def mult_matrix_det(self, matrix, det_m):
        for i in range(4):
            for j in range(4):
                matrix[i][j] *= det_m

        return matrix

    def calc_l(self, points_tab):
        L_matrix = []
        l1 = m.sqrt(pow((points_tab[1][0] - points_tab[0][0]), 2) + pow((points_tab[1][1] - points_tab[0][1]),
                                                                      2))
        l2 = m.sqrt(pow((points_tab[2][0] - points_tab[1][0]), 2) + pow((points_tab[2][1] - points_tab[1][1]),
                                                                      2))
        l3 = m.sqrt(pow((points_tab[3][0] - points_tab[2][0]), 2) + pow((points_tab[3][1] - points_tab[2][1]),
                                                                      2))
        l4 = m.sqrt(pow((points_tab[0][0] - points_tab[3][0]), 2) + pow((points_tab[0][1] - points_tab[3][1]),
                                                                      2))
        L_matrix.append(l1)
        L_matrix.append(l2)
        L_matrix.append(l3)
        L_matrix.append(l4)

        return L_matrix
