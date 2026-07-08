#!/usr/bin/env python3
import sp8py
import numpy as np
import matplotlib.pyplot as plt

######### <4,3>
lumo = 0.5; homo = 0.65; left = 4; right = 3
#########
######### <5,2>
# lumo = 0.65; homo = 0.8; left = 5; right = 2
#########
######### <3,4>
# lumo = 0.35; homo = 0.5; left = 3; right = 4
#########


lw = 0.5
mydashes = [6, 4]
fs_small=9
euler_color = 'tab:blue'

def plot_xtick(pos,mc,text,offset=0, ha='center'):
    plt.plot([pos,pos],[-0.01,0.01],'k-')
    plt.plot([pos,pos],[0, np.polyval(mc, pos)],'--',color='black', linewidth=lw, dashes=mydashes)
    plt.text(pos+offset,-0.02,text,ha=ha,va='top')


x = np.linspace(-0.1,1.1,200)
plt.figure()
plt.rc('font', size=9)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.quiver(-0.1,  0, 1.2, 0, angles='xy', scale_units='xy', scale=1, width=0.005)
plt.quiver( 0, 0,  0, 1.1, angles='xy', scale_units='xy', scale=1, width=0.005)


### SP8 STUFF
v = np.empty(9)
acc_left = True
acc_right = True
sp8_spec = (left,right,acc_left,acc_right)
info = sp8py.get_sp8_params(lumo,homo,sp8_spec,v)
if info != 0:
    print(f"get_sp8_params failed with info = {info}")
    exit(1)
print(f'Slope <{left},{right}>: {sp8py.sp8_prim(v,homo)}')
#    print(f'v: {[ "{:0.16f}".format(x) for x in v]}')
mc = np.empty(9)
sp8py.get_sp8_monomial_coefficients(v,mc)
print(mc)
xv = np.linspace(0,1,500)
yv = [np.polyval(mc, x) for x in xv]
plt.plot(xv,yv,zorder=3)
plt.text((homo+lumo)/2, np.polyval(mc, (homo+lumo)/2),f'$p(x)$',ha='right',va='bottom')
#plt.text((homo+lumo)/2, np.polyval(mc, (homo+lumo)/2),f'$p_{{({left},{right},\lambda_{{\mathrm{{lumo}}}},\lambda_{{\mathrm{{homo}}}})}}(x)$',ha='right',va='bottom')
## Plot s1...s7
#for ind,s in enumerate(v[2:]):
#    plt.plot([s,s],[-0.01,0.01],'k-')
#    plt.plot([s,s],[0, np.polyval(mc, s)],'--',color='black', linewidth=lw, dashes=mydashes)
#    plt.text(s,-0.01,f'$s_{ind+1}$',ha='center',va='top')
##
## Plot homo and lumo
plt.plot([lumo,lumo],[-0.01,0.01],'k-')
plt.plot([lumo,lumo],[0, np.polyval(mc, lumo)],'--',color='black', linewidth=lw, dashes=mydashes)
plt.text(lumo+0.01,-0.02,'$\lambda_{\mathrm{lumo}}$',ha='center',va='top')
#
plt.plot([homo,homo],[-0.01,0.01],'k-')
plt.plot([homo,homo],[0, np.polyval(mc, homo)],'--',color='black', linewidth=lw, dashes=mydashes)
plt.text(homo-0.015,-0.02,'$\lambda_{\mathrm{homo}}$',ha='center',va='top')
##
## Plot l1...l_{left}
lvec = v[1+left:1:-1]
print(lvec)
offset_vec = np.zeros(left)
offset_vec[left-1] = 0.005
for ind,l in enumerate(lvec):
    plot_xtick(l,mc,f'$s_{ind+1}$',offset=offset_vec[ind])
# plot_xtick(0.0,mc,f'$0\!=\!l_{left+1}$',offset=0.01, ha='right')
plot_xtick(0.0,mc,f'$0$',offset=0.0, ha='center')
#plt.plot([0, 0],[-0.065,0],'--',color='black', linewidth=lw, dashes=mydashes)
#plt.text(0.0,-0.065,f'$l_{left+1}=0$',ha='center',va='top')

#    plt.plot([l,l],[-0.01,0.01],'k-')
#    plt.plot([l,l],[0, np.polyval(mc, l)],'--',color='black', linewidth=lw, dashes=mydashes)
#    plt.text(l,-0.02,f'$l_{ind+1}$',ha='center',va='top')
##
## Plot l1...l_{right}
rvec = v[2+left:]
print(rvec)
for ind,r in enumerate(rvec):
    plot_xtick(r,mc,f'$s_{left+ind+1}$')
#plot_xtick(1.0,mc,f'$r_{right+1}\!=\!1$',offset=-0.01,ha='left')
plot_xtick(1.0,mc,f'$1$',offset=-0.0,ha='center')
#plt.plot([1, 1],[-0.065,0],'--',color='black', linewidth=lw, dashes=mydashes)
#plt.text(1.0,-0.065,f'$r_{right+1}=1$',ha='center',va='top')
##
## plot p(lumo) and p(homo)
plumo = np.polyval(mc, lumo)
plt.plot([-0.01,0.01],[plumo,plumo],'k-')
plt.text(-0.02,plumo-0.0,'$p(\lambda_{\mathrm{lumo}})$',ha='right',va='center')
plt.plot([0,1.0],[plumo, plumo],'--',color='tab:orange', linewidth=2*lw, dashes=mydashes)
phomo = np.polyval(mc, homo)
plt.plot([-0.01,0.01],[phomo,phomo],'k-')
plt.text(-0.02,phomo-0.0,'$p(\lambda_{\mathrm{homo}})$',ha='right',va='center')
plt.plot([0,1.0],[phomo, phomo],'--',color='tab:orange', linewidth=2*lw, dashes=mydashes)
plt.plot([-0.01,0.01],[1.0,1.0],'k-')
plt.text(-0.02,1.0-0.0,'$1$',ha='right',va='center')
plt.plot([0,1.0],[1.0,1.0],'--',color='black', linewidth=lw, dashes=mydashes)


###

plt.axis('equal')
plt.xlim([-0.2, 1.2])
plt.ylim([-0.1, 1.2])
plt.axis('off')

#plt.show()
plt.savefig(f"sp8_acc_schematic_{left}_{right}.pdf", bbox_inches='tight')

