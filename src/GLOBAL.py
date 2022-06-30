import numpy as np

eps0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
mu0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
c0 = 1./np.sqrt(eps0*mu0)

COLOR = ["blue", "green", "red", "orange", "purple", "yellow"]

LTS           : bool = True
INTERPOLATING : bool = False
BORDER_LTS    : bool = False
PLOT          : bool = True

NTIMESTEPS  : int = 1000000
NVERBOSE    : int = 10000
DELTACOARSE : float = 0.5
MAX_SG      : int = 1
CFL         : float = 0.9
PPW         : float = 20
INTERVAL    : int = 10
INTERP_OLD  : str = "1/4"
INTERP_NEW  : str = "1-interp_old"
BORDER_LTS_LOOP : type[np.inf] = np.inf

BOUNDARY_INI = "PEC"
BOUNDARY_END = "PEC"

PROBE_NAME = "..\\probe_sg{}_cfl{:1.3F}".format(MAX_SG, CFL)
if LTS:
    PROBE_NAME += "_LTS"
if INTERPOLATING:
    PROBE_NAME += "_INTERP"
if BORDER_LTS:
    PROBE_NAME += "_BORDERLTS"

if not LTS:
    CFL /= 2**MAX_SG
