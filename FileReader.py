def read(filename):
    tab = []
    with open(filename) as file:
        for line in file:
            tab.append(line.strip().split()[0])

    init_temp = float(tab[0])  # initial temperature
    sim_time = float(tab[1])  # simulation time
    step_time = float(tab[2])  # simulation step time
    amb_temp = float(tab[3])  # ambient temperature
    alfa = float(tab[4])  # alfa
    Height = float(tab[5])  # height
    Width = float(tab[6])  # width
    N_Height = int(tab[7])  # nodes per height
    N_Width = int(tab[8])  # nodes per Width
    spec_heat = float(tab[9])  # specific heat
    conduct = float(tab[10])  # conductivity
    dens = float(tab[11])  # density

    return init_temp, sim_time, step_time, amb_temp, alfa, Height, Width, N_Height, N_Width, spec_heat, conduct, dens
