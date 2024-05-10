#!/usr/bin/env python3
import sp8py
import numpy as np
import matplotlib.pyplot as plt

def get_max_slope(mu_vec,gap):
    lumo_vec = [mu-gap/2 for mu in mu_vec]
    homo_vec = [mu+gap/2 for mu in mu_vec]
    slope_vec = np.empty(len(mu_vec))
    for idx,(mu,lumo,homo) in enumerate(zip(mu_vec,lumo_vec,homo_vec)):
        v = np.empty(9)
        sp8py.get_sp8_params_max_slope(lumo,homo,v)
        slope_vec[idx] = sp8py.sp8_prim(v,mu)
    return slope_vec

npoints = 500
gap = 1e-5

fig_slope, ax_slope = plt.subplots()

mu_vec = np.linspace(0.001,0.999,npoints)
slope_vec = get_max_slope(mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec)

ax_slope.grid('on')
ax_slope.set_xlim((0, 1))
ax_slope.set_ylim((0, 15))
#ax_slope.legend()
fig_slope.savefig("max_slope_at_mu.pdf", bbox_inches='tight')
