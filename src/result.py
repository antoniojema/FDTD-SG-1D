import numpy as np
from matplotlib import ticker
import matplotlib.pyplot as plt
import re
from Fourier.functions import *
from GLOBAL import *

max_index = 20000
casos = [
    {"coletilla":"_cfl0.700_LTS_BORDERLTS.txt", "cfl":0.700, "label":"Border LTS"},
    {"coletilla":"_cfl0.700_LTS.txt"          , "cfl":0.700, "label":"LTS"},
    {"coletilla":"_cfl0.700.txt"              , "cfl":0.700, "label":"Non-LTS"},
]

print("Begin reading")
data_sg0 = np.zeros((len(casos),), dtype=object)
data_sg1 = np.zeros((len(casos),), dtype=object)
data_sg2 = np.zeros((len(casos),), dtype=object)
i = 0
for caso in casos:
    # COARSE LTS
    with open("..\\probe_sg0"+caso["coletilla"], "r") as fin:
        data_sg0[i] = np.array(list(map(lambda x: list(map(float, re.split("\s+", x)[1:3])), fin.readlines()[:max_index])))

    # SG1 LTS
    with open("..\\probe_sg1"+caso["coletilla"], "r") as fin:
        data_sg1[i] = np.array(list(map(lambda x: list(map(float, re.split("\s+", x)[1:3])), fin.readlines()[:max_index])))
        data_sg1[i][:, 1] -= data_sg0[i][:, 1]

    # SG2 LTS
    with open("..\\probe_sg2"+caso["coletilla"], "r") as fin:
        data_sg2[i] = np.array(list(map(lambda x: list(map(float, re.split("\s+", x)[1:3])), fin.readlines()[:max_index])))
        data_sg2[i][:, 1] -= data_sg0[i][:, 1]

    i += 1

del data_sg0,

# Gaussian source
print("Calculate source", flush=True)
sigma_source = ppw_to_sigma(PPW, delta=DELTACOARSE, c0=c0, dbDecay=3.)
Dt_source = 0.01*sigma_source
times = np.arange(0, 1000*Dt_source, Dt_source)
data_source = gaussian(times, amplitude=1., mean=5*sigma_source, sigma=sigma_source)
data_source = np.append(times.reshape((len(times),1)), data_source.reshape(len(times),1), axis=1)
#np.savetxt("gaussian.txt", data_source, fmt="%20.10E")

# Frecs
max_frec = ppw_to_freq(PPW, delta=DELTACOARSE, c=c0)
max_frec_log = np.log10(max_frec)
frecs = np.logspace(max_frec_log-8, max_frec_log+2, num=1000, base=10)
ppws = freq_to_ppw(frecs, delta=DELTACOARSE, c=c0)

# Fourier transform
print("Begin Fourier Transforms")

print("- Source", flush=True)
FT_source = np.absolute(FT_inhomogeneous(data_source[:, 0], data_source[:, 1], frecs))

FT_sg1 = np.zeros((len(casos),), dtype=object)
FT_sg2 = np.zeros((len(casos),), dtype=object)
for i in range(len(casos)):
    print("- "+casos[i]["label"], flush=True)
    print("-- SG1", flush=True)
    FT_sg1[i]= np.absolute(FT_inhomogeneous(data_sg1[i][:, 0], data_sg1[i][:, 1], frecs))
    print("-- SG2", flush=True)
    FT_sg2[i]= np.absolute(FT_inhomogeneous(data_sg2[i][:, 0], data_sg2[i][:, 1], frecs))

# Plot
print("Begin plotting", flush=True)

# Plot
#plt.rcParams.update({
#    "text.usetex": True,
#    "font.family": "sans-serif",
#    "font.serif": ["Helvetica"],
#    "font.weight": "bold"
#})
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

print("- Time", flush=True)
plt.figure("Time")
#plt.plot(data_source [:,0], data_source [:,1], label="Source")
for i in range(len(casos)):
    plt.plot(data_sg1[i][:,0], data_sg1[i][:,1], label="Max subgrid lvl 1 ("+casos[i]["label"]+") cfl-{:1.3F}".format(casos[i]["cfl"]))
    plt.plot(data_sg2[i][:,0], data_sg2[i][:,1], label="Max subgrid lvl 2 ("+casos[i]["label"]+") cfl-{:1.3F}".format(casos[i]["cfl"]))
plt.legend()
plt.grid()

print("- Fourier", flush=True)
fig, ax = plt.subplots(figsize=(10, 6))
#plt.semilogx(ppws, 20*np.log10(FT_source)       , label="source"           )
for i in range(len(casos)):
    plt.semilogx(ppws, 20*np.log10(FT_sg1[i]/FT_source), label=r"\textbf{Max subgrid lvl 1 ("+casos[i]["label"]+") cfl{:1.3F}".format(casos[i]["cfl"])+"}")
    plt.semilogx(ppws, 20*np.log10(FT_sg2[i]/FT_source), label=r"\textbf{Max subgrid lvl 2 ("+casos[i]["label"]+") cfl{:1.3F}".format(casos[i]["cfl"])+"}")

plt.legend()
plt.grid(which='major', ls="-")
plt.grid(which='minor', ls="--")
ax.axvspan(0.1, 5, alpha=0.5, color="gray")

plt.xlabel(r"\textbf{PPW}")
plt.xlim(1, 1e+6)
locmaj = ticker.LogLocator(base=10,numticks=12)
ax.xaxis.set_major_locator(locmaj)

plt.ylabel(r"\textbf{Reflection coef. (dB)}")
plt.yticks(np.arange(-400, +400, 20))
plt.ylim(-160, 40)

plt.tight_layout()
plt.show()
