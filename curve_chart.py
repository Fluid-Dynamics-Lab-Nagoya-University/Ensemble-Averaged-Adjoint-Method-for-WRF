"""
Plots curves that compare deterministic
and ensemble-averaged control results.

by Shan Jiang, FDL, Nagoya University

"""

import os
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'Helvetica',
    'font.size': 14,
    'axes.labelsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
})

# sg and ig value lists (must match those used when generating the npy files)
s_values = ["0", "0.001", "0.01", "0.05", "0.1", "0.15", "0.2", "0.3", "0.4", "0.5"]
ig_list  = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"]

idx_sg0  = s_values.index("0")
idx_sg03 = s_values.index("0.3")

# Load deltaP matrices
mat_ng0   = np.load("./dats/deltaP_mat_absG_i0M_NOV0.02_ng0_sg_G3x3.npy")
mat_ngEQsg = np.load("./dats/deltaP_mat_absG_i0M_NOV0.02_ngEQsg_sg_G3x3.npy")

x    = [float(v) for v in ig_list]
y_b  = mat_ng0[idx_sg0,  :].tolist()   # sg=0,   ng=0
y_a2 = mat_ng0[idx_sg03, :].tolist()   # sg=0.3, ng=0
y_a  = mat_ngEQsg[idx_sg03, :].tolist() # sg=0.3, ng=0.3

plt.plot(x, y_b,  marker='^', color='#0000FF', label='det. sens., det. valid.')
plt.plot(x, y_a2, marker='s', color='#FFA500', label='ens. sens., det. valid.')
plt.plot(x, y_a,  marker='o', color='#FF0000', label='ens. sens., ens. valid.')

plt.xlabel(r'$C_{\mathrm{act}}$')
plt.ylabel(r'$-\Delta \mathbf{P}^f$ (mm)')

plt.xticks(np.arange(0, 1.0, 0.1))
plt.legend(fontsize=12)
plt.grid(True)

os.makedirs("./figs", exist_ok=True)
plt.savefig('./figs/curve_chart.png', dpi=600, bbox_inches='tight')