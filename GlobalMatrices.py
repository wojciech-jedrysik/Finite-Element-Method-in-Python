class Global_H():
    def __init__(self):
        self.gH = []

    def create(self, elements, nodes_quantity):
        for i in range(nodes_quantity):
            line = []
            for j in range(nodes_quantity):
                line.append(0)
            self.gH.append(line)
        for element in elements:
            ID_list = []
            for node in element.nodes:
                ID_list.append(node.id)
            for i in range(4):
                for j in range(4):
                    self.gH[ID_list[i]-1][ID_list[j]-1] += element.local_H[i][j]


class Global_C():
    def __init__(self):
        self.gC = []

    def create(self, elements, nodes_quantity):
        for i in range(nodes_quantity):
            line = []
            for j in range(nodes_quantity):
                line.append(0)
            self.gC.append(line)
        for element in elements:
            ID_list = []
            for node in element.nodes:
                ID_list.append(node.id)
            for i in range(4):
                for j in range(4):
                    self.gC[ID_list[i]-1][ID_list[j]-1] += element.local_C[i][j]


class Global_P():
    def __init__(self):
        self.gP = []

    def create(self, elements, nodes_quantity):
        for i in range(nodes_quantity):
            self.gP.append(0)
        for element in elements:
            ID_list = []
            for node in element.nodes:
                ID_list.append(node.id)
            for i in range(4):
                self.gP[ID_list[i]-1] += element.local_P[i]



