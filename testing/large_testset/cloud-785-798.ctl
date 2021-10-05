
# Emitters...
NG = 13
EMITTER[0] = CCl4
EMITTER[1] = ClO
EMITTER[2] = ClONO2
EMITTER[3] = CO2
EMITTER[4] = F11
EMITTER[5] = F113
EMITTER[6] = F22
EMITTER[7] = H2O
EMITTER[8] = HNO3
EMITTER[9] = HNO4
EMITTER[10] = NO2
EMITTER[11] = O3
EMITTER[12] = PAN

# Channels...
ND = 33
NU[0] = 785.00
NU[1] = 785.42
NU[2] = 785.84
NU[3] = 786.26
NU[4] = 786.68
NU[5] = 787.10
NU[6] = 787.52
NU[7] = 787.94
NU[8] = 788.36
NU[9] = 788.78
NU[10] = 789.20
NU[11] = 789.62
NU[12] = 790.04
NU[13] = 790.46
NU[14] = 790.88
NU[15] = 791.30
NU[16] = 791.72
NU[17] = 792.14
NU[18] = 792.56
NU[19] = 792.98
NU[20] = 793.40
NU[21] = 793.82
NU[22] = 794.24
NU[23] = 794.66
NU[24] = 795.08
NU[25] = 795.50
NU[26] = 795.92
NU[27] = 796.34
NU[28] = 796.76
NU[29] = 797.18
NU[30] = 797.60
NU[31] = 798.02
NU[32] = 798.44

TBLBASE = /p/fastdata/slmet/slmet111/model_data/jurassic/tab/crista-nf/ascii_785_964/crista_nf 

# Continua...
CTM_CO2 = 1
CTM_H2O = 1
CTM_N2 = 0
CTM_O2 = 0

# Aerosol/Clouds...
SCA_N = 1
SCA_MULT = 1
SCA_EXT = beta_e

MAX_QUEUE = 1e6

# Raytracing...
RAYDS = 10.
RAYDZ = 0.1

# Atmosphere/Climatology...
ZMIN = 0
ZMAX = 80
DZ = 1
CLIMPATH = /p/project/chwu36/hwu361/crista-nf/setup/atmo/clim_pscs.tab

# use the GPU: 0:never, 1:always, -1:if possible
USEGPU = 1 
