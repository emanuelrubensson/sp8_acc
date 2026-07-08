#!/usr/bin/env python3
import sp8py
import sp8py.sp7 as sp7
import sp8py.sp5 as sp5
import numpy as np
import matplotlib.pyplot as plt

def get_slope(left,right,mu_vec,gap):
    acc_left = left != 0
    acc_right = right != 0
    sp8_spec = (left,right,acc_left,acc_right)
    lumo_vec = [mu-gap/2 for mu in mu_vec]
    homo_vec = [mu+gap/2 for mu in mu_vec]
    slope_vec = np.empty(len(mu_vec))
    map_vec = np.empty(len(mu_vec))
    for idx,(mu,lumo,homo) in enumerate(zip(mu_vec,lumo_vec,homo_vec)):
        v = np.empty(9)
        sp8py.get_sp8_params(lumo,homo,sp8_spec,v)
        slope_vec[idx] = sp8py.sp8_prim(v,mu)
        map_vec[idx] = sp8py.sp8_eval(v,mu)
    return slope_vec, map_vec

def get_slope_sp7(lumo):
    mc,_ = sp7.sp7_polynomial_coeffs(lumo)
    mcder = np.polyder(mc)
    return np.polyval(mcder,0.5)

def get_slope_sp5(lumo):
    mc = sp5.sp5_polynomial_coeffs(lumo)
    mcder = np.polyder(mc)
    return np.polyval(mcder,0.5)

npoints = 150
gap = 1e-5

plt.rc('font', size=16)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig_slope, ax_slope = plt.subplots()
fig_map, ax_map = plt.subplots()

mu_vec = np.linspace(0.001,0.4,npoints)
slope_vec, map_vec = get_slope(0, 7, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(0,7)')
ax_map.plot(mu_vec,map_vec,label='(0,7)')

mu_vec = np.linspace(0.001,0.6,npoints)
slope_vec, map_vec = get_slope(1, 6, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(1,6)')
ax_map.plot(mu_vec,map_vec,label='(1,6)')

mu_vec = np.linspace(0.001,0.8,npoints)
slope_vec, map_vec = get_slope(2, 5, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(2,5)')
ax_map.plot(mu_vec,map_vec,label='(2,5)')

mu_vec = np.linspace(0.05,1-0.05,npoints)
slope_vec, map_vec = get_slope(3, 4, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(3,4)')
ax_map.plot(mu_vec,map_vec,label='(3,4)')

mu_vec = np.linspace(0.05,1-0.05,npoints)
slope_vec, map_vec = get_slope(4, 3, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(4,3)')
ax_map.plot(mu_vec,map_vec,label='(4,3)')

mu_vec = np.linspace(0.2,1-0.001,npoints)
slope_vec, map_vec = get_slope(5, 2, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(5,2)')
ax_map.plot(mu_vec,map_vec,label='(5,2)')

mu_vec = np.linspace(0.4,1-0.001,npoints)
slope_vec, map_vec = get_slope(6, 1, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(6,1)')
ax_map.plot(mu_vec,map_vec,label='(6,1)')

mu_vec = np.linspace(0.6,1-0.001,npoints)
slope_vec, map_vec = get_slope(7, 0, mu_vec, gap)
ax_slope.plot(mu_vec,slope_vec,label='(7,0)')
ax_map.plot(mu_vec,map_vec,label='(7,0)')

slope_sp5 = get_slope_sp5(0.5-gap/2)
ax_slope.plot(0.5,slope_sp5,'o', label='Su5')

slope_sp7 = get_slope_sp7(0.5-gap/2)
ax_slope.plot(0.5,slope_sp7,'x', label='Su7')

ax_slope.grid('on')
ax_slope.set_xlim((0, 1))
ax_slope.set_ylim((0, 14))
ax_slope.legend(ncol=2,fontsize='small')
#ax_slope.set_xlabel('$\mu = (\lambda_{\mathrm{lumo}}+\lambda_{\mathrm{homo}})/2$')
#ax_slope.set_ylabel("$p_{L,R,\lambda_{\mathrm{lumo}},\lambda_{\mathrm{homo}}}'(\mu)$")
ax_slope.set_xlabel('$\mu$')
ax_slope.set_ylabel("$p_{(L,R,\mu-\\varepsilon,\mu+\\varepsilon)}'(\mu)$")
fig_slope.savefig("slope_at_mu.pdf", bbox_inches='tight')

ax_map.grid('on')
ax_map.set_xlim((0, 1))
ax_map.set_ylim((0, 1))
ax_map.legend(ncol=2,fontsize='small',loc='upper center')
#ax_map.set_xlabel('$\mu = (\lambda_{\mathrm{lumo}}+\lambda_{\mathrm{homo}})/2$')
#ax_map.set_ylabel("$p_{L,R,\lambda_{\mathrm{lumo}},\lambda_{\mathrm{homo}}}(\mu)$")
ax_map.set_xlabel('$\mu$')
ax_map.set_ylabel("$p_{(L,R,\mu-\\varepsilon,\mu+\\varepsilon)}(\mu)$")
fig_map.savefig("map_at_mu.pdf", bbox_inches='tight')

#plt.show()
