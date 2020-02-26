import Structures as S
import math as m

universal_element = S.UniversalElement()
points_tab = [[0, 0], [0.025, 0], [0.025, 0.025], [0, 0.025]]
conduct = 30
spec_heat = 700
dens = 7800
alfa = 25
are_edges = [1, 0, 0, 0]

class matrixH():
    def __init__(self):
        self.H = []

    def calc_H(self, points_tab, K):
        J_matrices = []
        for i in range(4):
            J_matrices.append(universal_element.calc_jacobian(points_tab, i))
        all_det_J = []
        for i in range(4):
            all_det_J.append(universal_element.calc_det(J_matrices[i]))
        all_J_1 = []
        for i in range(4):
            all_J_1.append(universal_element.calc_J_1(J_matrices[i], all_det_J[i]))
        byX, byY = universal_element.calc_deXY(all_J_1)
        dx_dx_T = []
        for i in range(4):
            dx_dx_T.append(universal_element.multiply_vector(byX[i]))
        dy_dy_T = []
        for i in range(4):
            dy_dy_T.append(universal_element.multiply_vector(byY[i]))
        for i in range(4):
            dx_dx_T[i] = universal_element.mult_matrix_det(dx_dx_T[i], all_det_J[i])
        for i in range(4):
            dy_dy_T[i] = universal_element.mult_matrix_det(dy_dy_T[i], all_det_J[i])
        K_matrices = []
        for i in range(4):
            pc_matrix = []
            for j in range(4):
                temp = []
                for k in range(4):
                    temp.append(K * (dx_dx_T[i][j][k] + dy_dy_T[i][j][k]))
                pc_matrix.append(temp)
            K_matrices.append(pc_matrix)
        for i in range(4):
            temp = []
            for j in range(4):
                t = 0
                for k in range(4):
                    t += K_matrices[k][i][j]
                temp.append(t)
            self.H.append(temp)


class matrixC():
    def __init__(self):
        self.C = []

    def calc_C(self, points_tab, c, ro):
        J_matrices = []
        for i in range(4):
            J_matrices.append(universal_element.calc_jacobian(points_tab, i))
        all_det_J = []
        for i in range(4):
            all_det_J.append(universal_element.calc_det(J_matrices[i]))
        C_matrices = []
        for i in range(4):
            matrix = []
            for j in range(4):
                temp = []
                for k in range(4):
                    temp.append(universal_element.shapeFunctionsValues[i][j] * universal_element.shapeFunctionsValues[i][k] * all_det_J[i] * c * ro)
                matrix.append(temp)
            C_matrices.append(matrix)
        for i in range(4):
            temp = []
            for j in range(4):
                t = 0
                for k in range(4):
                    t += C_matrices[k][i][j]
                temp.append(t)
            self.C.append(temp)


class matrixH_bc():
    def __init__(self):
        self.H_bc = []

    def calc_H_bc(self, points_tab, Alfa, are_edges):
        L_matrix = universal_element.calc_l(points_tab)
        pow_matrices = []
        v1, v2 = 1 / m.sqrt(3), 1
        ksi_eta = [[[-v1, -v2], [v1, -v2]], [[v2, -v1], [v2, v1]], [[v1, v2], [-v1, v2]], [[-v2, v1], [-v2, -v1]]]
        for i in range(4):
            matrix = []
            for j in range(2):
                N1 = 0.25 * (1 - ksi_eta[i][j][0]) * (1 - ksi_eta[i][j][1])
                N2 = 0.25 * (1 + ksi_eta[i][j][0]) * (1 - ksi_eta[i][j][1])
                N3 = 0.25 * (1 + ksi_eta[i][j][0]) * (1 + ksi_eta[i][j][1])
                N4 = 0.25 * (1 - ksi_eta[i][j][0]) * (1 + ksi_eta[i][j][1])
                matrix.append([N1, N2, N3, N4])
            pow_matrices.append(matrix)

        H_bc_matrices = []
        matrices1 = []
        matrices2 = []
        for i in range(4):
            matrix1 = []
            matrix2 = []
            for j in range(4):
                temp1 = []
                tmp1 = pow_matrices[i][0][j]
                temp2 = []
                tmp2 = pow_matrices[i][1][j]
                for k in range(4):
                    temp1.append(tmp1 * pow_matrices[i][0][k] * Alfa)
                    temp2.append(tmp2 * pow_matrices[i][1][k] * Alfa)
                matrix1.append(temp1)
                matrix2.append(temp2)
            matrices1.append(matrix1)
            matrices2.append(matrix2)
        for i in range(4):
            Hbc = []
            for j in range(4):
                temp = []
                for k in range(4):
                    temp.append(
                        (matrices1[i][j][k] + matrices2[i][j][k]) * (L_matrix[i] / 2))
                Hbc.append(temp)
            H_bc_matrices.append(Hbc)

        for i in range(4):
            temp = []
            for j in range(4):
                t = 0
                for k in range(4):
                    t += H_bc_matrices[k][i][j] * are_edges[k]
                temp.append(t)
            self.H_bc.append(temp)


class vectorP():
    def __init__(self):
        self.P = [0, 0, 0, 0]

    def calc_P(self, nodes, temp, Alfa):
        surface_N = universal_element.calc_surface_N()
        for i in range(4):
            if nodes[i].is_out is True and nodes[(i+1)%4].is_out is True:
                tmp = 0
                tmp += pow((nodes[i].x - nodes[(i+1)%4].x), 2)
                tmp += pow((nodes[i].y - nodes[(i+1)%4].y), 2)
                tmp = m.sqrt(tmp) * 0.5
                for j in range(4):
                    self.P[i] += surface_N[i][j] * temp * Alfa * tmp
                    self.P[(i+1)%4] += surface_N[i][j] * temp * Alfa * tmp
