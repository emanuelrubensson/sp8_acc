#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

lumo = 0.4
homo = 1-lumo
beta = lumo    
alpha = 3 / np.sqrt(12*beta**2 - 18*beta + 9)
lw = 0.5
mydashes = [6, 4]

def sp3acc(x):
    xs = alpha*(x-0.5) + 0.5
    return 3*xs**2 - 2*xs**3

def plot_xtick(pos,text,offset=0, ha='center'):
    plt.plot([pos,pos],[-0.01,0.01],'k-')
    plt.text(pos+offset,-0.02,text,ha=ha,va='top')

plt.figure()
plt.rc('font', size=14)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.quiver(-0.1,  0, 1.2, 0, angles='xy', scale_units='xy', scale=1, width=0.005)
plt.quiver( 0, 0,  0, 1.1, angles='xy', scale_units='xy', scale=1, width=0.005)

xv = np.linspace(0,1,500)

yv = [sp3acc(x) for x in xv]
plt.plot(xv,yv,zorder=3)#, label='$2x-x^2$')

# Plot dotted unscaled McWeeny curve
xv_scaled = np.linspace(alpha*(0-0.5) + 0.5,alpha*(1.0-0.5) + 0.5,500)
plt.plot(xv_scaled,3*xv_scaled**2-2*xv_scaled**3,':',zorder=3,color='tab:blue',linewidth=lw)#, label='$2x-x^2$')


plot_xtick(0.0,f'$0$',offset=0.0, ha='center')
plot_xtick(1.0,f'$1$',offset=-0.0,ha='center')

plt.plot([-0.01,0.01],[1.0,1.0],'k-')
plt.text(-0.02,1.0-0.0,'$1$',ha='right',va='center')

#plt.legend(bbox_to_anchor=(0.22,0.9), loc="upper left")

## Plot lumo and homo
plumo = sp3acc(lumo)
phomo = sp3acc(homo)
## x-axis
plt.plot([lumo,lumo],[-0.01,0.01],'k-')
plt.plot([lumo,lumo],[0, plumo],'--',color='black', linewidth=lw, dashes=mydashes)
plt.text(lumo+0.0,-0.02,'$\lambda_{\mathrm{lumo}}$',ha='center',va='top')
#
plt.plot([homo,homo],[-0.01,0.01],'k-')
plt.plot([homo,homo],[0, phomo],'--',color='black', linewidth=lw, dashes=mydashes)
plt.text(homo-0.0,-0.02,'$\lambda_{\mathrm{homo}}$',ha='center',va='top')
## y-axis
plt.plot([-0.01,0.01],[plumo,plumo],'k-')
plt.text(-0.02,plumo-0.0,'$p(\lambda_{\mathrm{lumo}})$',ha='right',va='center')
plt.plot([0,1.0],[plumo, plumo],'--',color='tab:orange', linewidth=2*lw, dashes=mydashes)
#
plt.plot([-0.01,0.01],[phomo,phomo],'k-')
plt.text(-0.02,phomo-0.0,'$p(\lambda_{\mathrm{homo}})$',ha='right',va='center')
plt.plot([0,1.0],[phomo, phomo],'--',color='tab:orange', linewidth=2*lw, dashes=mydashes)

##
plt.plot([-0.01,0.01],[1.0,1.0],'k-')
plt.text(-0.02,1.0-0.0,'$1$',ha='right',va='center')
plt.plot([0,1.0],[1.0,1.0],'--',color='black', linewidth=lw, dashes=mydashes)

##
x = 0.66
plt.text(x, sp3acc(x),f'$p(x)$',ha='right',va='bottom')

plt.axis('equal')
plt.xlim([-0.2, 1.2])
plt.ylim([-0.1, 1.2])
plt.axis('off')

#plt.show()
plt.savefig(f"sp3_acc_schematic.eps", bbox_inches='tight')
