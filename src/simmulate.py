import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from GLOBAL import *
from fdtd import *
from Fourier.functions import ppw_to_sigma
from Math import gaussian

######################
###   PARAMATERS   ###
######################

# Numero de celdas coarse
Ncells_coarse = 20

# Tamaño de la malla coarse
cfl = CFL
Ds_coarse = DELTACOARSE
Dt_coarse = cfl*Ds_coarse/c0

# Numero de instantes temporales
Ntimesteps = NTIMESTEPS
n_verbose =  NVERBOSE

# Numero maximo de grid
max_SG = MAX_SG

# Rango que cubre cada nivel, respecto del inmediatamente superior
# [xi, xf] (lvl(0)), ..., [xi,xf] (lvl(max-1))
if max_SG in [0, 2]:
    fine_range = np.array([
        [4, 6],
        [10, 15],
    ], dtype=int)
elif max_SG == 1:
    fine_range = np.array([
        [10, 15],
    ], dtype=int)
elif max_SG == 3:
    fine_range = np.array([
        [3, 5],
        [4, 8],
        [10, 16],
    ], dtype=int)


##################
###   PROBES   ###
##################

Nprobes = 1
probe_name = [
    PROBE_NAME
]
probe_index = np.array([
    2,
], dtype=int)
probe_Dt = np.array([
    0.99*Dt_coarse,
])


############################
###   GAUSSIAN SOURCES   ###
############################

Ngaussians = 1
gaussian_index = np.array([
    4,
], dtype=int)
gaussian_direction = np.array([
    +1,
], dtype=int)
gaussian_amplitude = np.array([
    1,
], dtype=float)
gaussian_mean = np.array([
    5*ppw_to_sigma(PPW, delta=Ds_coarse, c0=c0, dbDecay=3.),
], dtype=float)
gaussian_spread = np.array([
    ppw_to_sigma(PPW, delta=Ds_coarse, c0=c0, dbDecay=3.),
], dtype=float)


########################
###   FIRST CHECKS   ###
########################

Ncells = np.zeros((max_SG+1,), dtype=int)
Ncells[max_SG] = Ncells_coarse
for sg in range(max_SG-1, -1, -1):
    if fine_range[sg][0] < 0:
        raise ValueError("Grid lvl-{} cannot have lower bound less than 0.".format(sg))
    if fine_range[sg][1] > Ncells[sg+1]:
        raise ValueError("Grid lvl-{} cannot have upper bound ({}) greater than grid lvl-{} cells number ({})".\
            format(sg, fine_range[sg][1], sg+1, Ncells[sg+1]))
    Ncells[sg] = 2 * (fine_range[sg][1] - fine_range[sg][0])


#########################################
###   ALLOCATE FIELDS AND CONSTANTS   ###
#########################################

E = np.zeros((max_SG+1,), dtype=object)
H = np.zeros((max_SG+1,), dtype=object)
Ds = np.zeros((max_SG+1,))
Dt = np.zeros((max_SG+1,))
CE = np.zeros((max_SG+1,))
CH = np.zeros((max_SG+1,))
CE_border = np.zeros((max_SG+1,))
H_border  = np.zeros((max_SG+1, 2))
n_E = np.zeros((max_SG+1,), dtype=int)
for sg in range(max_SG+1):
    E[sg] = np.zeros((Ncells[sg]+1,))
    H[sg] = np.zeros((Ncells[sg]  ,))
    Ds[sg] = Ds_coarse / 2.**sg
    Dt[sg] = Dt_coarse / 2.**sg
    if LTS:
        CE[sg]        = Dt[sg]/(eps0*Ds[sg])
        CH[sg]        = Dt[sg]/(mu0 *Ds[sg])
        if BORDER_LTS:
            CE_border[sg] = Dt[sg]/(eps0*Ds[sg])
        else:
            CE_border[sg] = Dt[sg]/(eps0*1.5*Ds[sg])
    else:
        CE[sg]        = Dt_coarse/(eps0    *Ds[sg])
        CH[sg]        = Dt_coarse/(mu0     *Ds[sg])
        CE_border[sg] = Dt_coarse/(eps0*1.5*Ds[sg])
    H_border [sg, :] = [0., 0.]

E1_mur = 0.
E2_mur = 0.
CMur = (c0*Dt_coarse-Ds_coarse)/(c0*Dt_coarse+Ds_coarse)


##########################
###   PREPARE PROBES   ###
##########################

probe = [[] for i in range(Nprobes)]
probe_t = [[] for i in range(Nprobes)]
probe_next_t = np.zeros((Nprobes,))
probe_n = np.zeros((Nprobes,), dtype=int)


#########################
###   TIME STEPPING   ###
#########################

if PLOT:
    # Plot stuff
    fig, ax = plt.subplots()

    xdata_ini = [None for _ in range(max_SG+1)]
    xdata_fin = [None for _ in range(max_SG+1)]
    ydata_ini = [None for _ in range(max_SG+1)]
    ydata_fin = [None for _ in range(max_SG+1)]
    ln_ini    = [None for _ in range(max_SG+1)]
    ln_fin    = [None for _ in range(max_SG+1)]
    begin = 0.
    for sg in range(max_SG, 0, -1):
        xdata_ini[sg] = np.arange(
            begin,
            begin + fine_range[sg-1][0]/2**(max_SG-sg),
            1/2**(max_SG-sg)
        )
        xdata_fin[sg] = np.arange(
            begin + fine_range[sg-1][1]/2**(max_SG-sg) + 1/2**(max_SG-sg),
            begin +          Ncells[sg]/2**(max_SG-sg) + 1/2**(max_SG-sg),
            1/2**(max_SG-sg)
        )
        ydata_ini[sg] = np.zeros(xdata_ini[sg].shape)
        ydata_fin[sg] = np.zeros(xdata_fin[sg].shape)

        begin += fine_range[sg-1][0] / 2**(max_SG-sg)

        ln_ini[sg] = plt.plot(xdata_ini[sg], ydata_ini[sg], ".-", color=COLOR[max_SG-sg])[0]
        ln_fin[sg] = plt.plot(xdata_fin[sg], ydata_fin[sg], ".-", color=COLOR[max_SG-sg])[0]

    sg = 0
    xdata_ini[0] = np.arange(
        begin,
        begin + Ncells[0]/2**max_SG + 1/2**max_SG,
        1/2**max_SG
    )
    ydata_ini[0] = np.zeros(xdata_ini[0].shape)
    ln_ini[0] = plt.plot(xdata_ini[0], ydata_ini[0], ".-", color=COLOR[max_SG-sg])[0]

    text = plt.text(.95, .95, "N = 0", bbox={"facecolor":'w', "pad":5}, ha="right", va="top", transform=ax.transAxes)
    ax.set_xlim(0, Ncells_coarse)
    ax.set_ylim(-1.2, 1.2)
    plt.grid()

t = 0.
#for n_coarse_time in range(Ntimesteps):
def advance(n_coarse):
    global E1_mur, E2_mur, t, max_SG
    if PLOT:
        global xdata_fin, xdata_ini, ln_fin, ln_ini, text

    if n_coarse % n_verbose == 0:
        print(" N = {} - Max E value: {:1.5E}".format(n_coarse, np.abs(E[max_SG][:]).max()))
        if any(np.abs(E[max_SG][:]) > 10):
            print ("!!!! UNSTABILITIES !!!!")

    # Advance E
    advanceE(E, H, CE, CH, CE_border, H_border, fine_range, n_E, max_SG, max_SG, from_E=False)
    E1_mur, E2_mur = advanceMur(E, E1_mur, E2_mur, CMur, max_SG)

    # Source E
    for i in range(Ngaussians):
        E[max_SG][gaussian_index[i]] += \
            gaussian(t, gaussian_mean[i], gaussian_spread[i], gaussian_amplitude[i])

    # Advance t
    t += 0.5*Dt_coarse

    # Advance H
    advanceH(E, H, CE, CH, CE_border, H_border, fine_range, n_E, max_SG, max_SG)

    # Source H
    for i in range(Ngaussians):
        sum_index = int(-(gaussian_direction[i]+1)/2)
        H[max_SG][gaussian_index[i]+sum_index] += gaussian_direction[i] * \
            gaussian(t, gaussian_mean[i]+0.5*Dt_coarse, gaussian_spread[i], gaussian_amplitude[i]/CE[max_SG]*cfl)

    # Advance t
    t += 0.5 * Dt_coarse

    # Probes
    for i in range(Nprobes):
        if t > probe_next_t[i]:
            probe_n[i] += 1
            probe_next_t[i] += probe_Dt[i]
            probe_t[i].append(t)
            probe  [i].append(E[max_SG][probe_index[i]])

    if PLOT:
        # Plot stuff
        for sg in range(max_SG, 0, -1):
            ydata_ini[sg] = E[sg][0:fine_range[sg-1][0]]
            ydata_fin[sg] = E[sg][fine_range[sg-1][1]+1:]
            ln_ini[sg].set_data(xdata_ini[sg], ydata_ini[sg])
            ln_fin[sg].set_data(xdata_fin[sg], ydata_fin[sg])
        ydata_ini[0] = E[0][:]
        ln_ini[0].set_data(xdata_ini[0], ydata_ini[0])
        text.set_text("N = {}".format(n_coarse))

if PLOT:
    ani = FuncAnimation(fig, advance, frames=Ntimesteps, interval=INTERVAL, repeat=False)
    plt.show()
else:
    for i in range(Ntimesteps):
        advance(i)

# Export probes
for i in range(Nprobes):
    p = np.append(
        np.array(probe_t[i]).reshape((probe_n[i],1)),
        np.array(probe  [i]).reshape((probe_n[i],1)),
        axis=1
    )
    np.savetxt(probe_name[i]+".txt", p, fmt="%20.10E")
