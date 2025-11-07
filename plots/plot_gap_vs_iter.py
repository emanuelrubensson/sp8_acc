#!/usr/bin/env python3
import sp8py
import numpy as np
import matplotlib.pyplot as plt

#input: a gap
#output: sp8 polynomial for maximizing gap
def get_greedy_sp8_mc(gapvec, acc_left=False, acc_right=False):
    v = np.empty(9)
    mc_opt = np.zeros(9)
    # loop for greedy choice of <left,right>
    for i in range(8):
        mc = np.empty(9)
        left = i
        right = 7-i
        sp8_spec = (left,right,acc_left,acc_right)
        #
        info = sp8py.get_sp8_params(gapvec[0],gapvec[1],sp8_spec,    v) #compute best poly for current root setup
        if info != 0:
            #print(f"get_sp8_params failed with info = {info}")
            sp8py.get_sp8_params(gapvec[0],gapvec[1],(left,right,False,False),    v)
            #print(acc_left)
            #mc = get_greedy_sp8_mc(gapvec, acc_left=False, acc_right=False)
            #exit(1)
        sp8py.get_sp8_monomial_coefficients(v, mc)

        if np.abs(np.polyval(mc, gapvec[1]) - np.polyval(mc, gapvec[0])) > np.abs(np.polyval(mc_opt, gapvec[1]) - np.polyval(mc_opt, gapvec[0])):
            #print(f'Slope <{left},{right}>: {sp8py.sp8_prim(v,gapvec[1])}')
            mc_opt = mc
    return mc_opt


def traceadjusted_p2(X, K=100):
    tr = 1.0
    n = len(X)
    hist = np.zeros((K+1, n))
    hist[0,:] = X
    for k in range(K):
        if sum(hist[k,:]) >= tr:
            hist[k+1,:] = hist[k,:]**2
        else:
            hist[k+1,:] = 2*hist[k,:] - hist[k,:]**2
    return hist

def mcweeny(gapvec, K=100):
    n = len(gapvec)
    hist = np.zeros((K+1, n))
    hist[0,:] = gapvec
    if not(hist[0,0] <= 0.5 and hist[0,1] >= 0.5):
        mu = (hist[0,1] + hist[0,0])/2
        scaling = 0.5/mu
        hist[0,:] = scaling * hist[0,:]
    
    for k in range(K):
        hist[k+1,:] = 3*hist[k,:]**2 - 2*hist[k,:]**3
    return hist

def sp8_gapmax(gapvec, K=100, acc_left=False, acc_right=False):
    puretol = 0.5
    n = len(gapvec)
    hist = np.zeros((K+1, n))
    hist[0,:] = gapvec
    for k in range(K):
        if (hist[k,1] - hist[k,0]) < puretol:
            mc_opt = get_greedy_sp8_mc(hist[k,:], acc_left=acc_left, acc_right=acc_right)
            hist[k+1,:] = np.polyval(mc_opt, hist[k,:])
        else:
            #make three p2 iterations if we are past purification phase in this case
            hist[k+1,:] = traceadjusted_p2(hist[k,:], K=3)[3]
    return hist

#########################
# Fixed gap data
#########################

# no. mults
K = 6*12

#initial gap
gapvec = [0.499, 0.501]
tmpgap = round(gapvec[1]-gapvec[0], 6)

hist_sp2 = traceadjusted_p2(gapvec, K=K)
hist_mw = mcweeny(gapvec, K=int(K/2))
hist_sp8 = sp8_gapmax(gapvec, K=int(K/3))
hist_sp8acc = sp8_gapmax(gapvec, K=int(K/3), acc_left=True, acc_right=True)

gap_sp2 = np.abs(hist_sp2[:,1] - hist_sp2[:,0] - np.ones(len(hist_sp2[:,0])))
gap_mw = np.abs(hist_mw[:,1] - hist_mw[:,0] - np.ones(len(hist_mw[:,0])))
gap_sp8 = np.abs(hist_sp8[:,1] - hist_sp8[:,0] - np.ones(len(hist_sp8[:,0])))
gap_sp8acc = np.abs(hist_sp8acc[:,1] - hist_sp8acc[:,0] - np.ones(len(hist_sp8acc[:,0])))


gap_mw = np.repeat(gap_mw, 2)[:-1]
gap_sp8 = np.repeat(gap_sp8, 3)[:-2]
gap_sp8acc = np.repeat(gap_sp8acc, 3)[:-2]

multsax = np.array([i for i in range(K+1)])
data_full_vary_mu = np.column_stack((multsax, gap_sp2, gap_mw, gap_sp8, gap_sp8acc))
np.savetxt('data/idotemponcy_evolution.txt', data_full_vary_mu, fmt='%.18f', delimiter=' ', header='x y1 y2 y3 y4', comments='')


############
# Plotting
############

#plotting gap
plt.plot(hist_sp2, label="sp2", color="green")
plt.plot([i*2 for i in range(int(K/2+1))], hist_mw, label="McWeeny", color="red")
plt.plot([i*3 for i in range(int(K/3+1))], hist_sp8, label = "sp8", color="blue")
plt.plot([i*3 for i in range(int(K/3+1))], hist_sp8acc, label = "sp8acc", color="cyan")

plt.legend()
plt.xlabel(xlabel="matmuls")
plt.ylabel(ylabel="lumo/homo")

plt.title(f"lumo/homo Evolution, gap={tmpgap}")
plt.savefig("lumo-homo_vs_matmults.pdf", bbox_inches='tight')

plt.close()

#plotting gapsize
plt.plot(gap_sp2, label="sp2")
plt.plot(multsax, gap_mw, label="McWeeny")
plt.plot(multsax, gap_sp8, label="sp8")
plt.plot(multsax, gap_sp8acc, label = "sp8acc", color="cyan")
plt.yscale('log')
plt.legend()
plt.title(f"Gap Evolution, gap={tmpgap}")
plt.xlabel(xlabel="matmuls")
plt.ylabel(ylabel="gap")
plt.savefig("gap_vs_matmults.pdf", bbox_inches='tight')
plt.close()




############################
# Iterations to convergence
############################
K=6*25
Ng = 100 #number of gaps
tol = 1e-15
gapsize = 1e-10
vary_mu = np.linspace(0.5, 1-gapsize, num=Ng) #center of gaps
no_mults_sp2 = np.zeros(Ng)
no_mults_mw = np.zeros(Ng)
no_mults_sp8 = np.zeros(Ng)
no_mults_sp8acc = np.zeros(Ng)

for i in range(Ng):
    gapveci = np.array([vary_mu[i]-gapsize/2, vary_mu[i]+gapsize/2])
    hist_sp2 = traceadjusted_p2(gapveci, K=K)
    hist_mw = mcweeny(gapveci, K=int(K/2))
    hist_sp8 = sp8_gapmax(gapveci, K=int(K/3))
    hist_sp8acc = sp8_gapmax(gapveci, K=int(K/3), acc_left=True, acc_right=True)

    gap_sp2 = hist_sp2[:,1] - hist_sp2[:,0]
    gap_mw = hist_mw[:,1] - hist_mw[:,0]
    gap_sp8 = hist_sp8[:,1] - hist_sp8[:,0]
    gap_sp8acc = hist_sp8acc[:,1] - hist_sp8acc[:,0]

    # get number of mults to convergence:
    no_mults_sp2[i]    = np.where(gap_sp2>(1-tol))[0][0]
    no_mults_mw[i]     = 2*np.where(gap_mw>(1-tol))[0][0]
    no_mults_sp8[i]    = 3*np.where(gap_sp8>(1-tol))[0][0]
    no_mults_sp8acc[i] = 3*np.where(gap_sp8acc>(1-tol))[0][0]

############
# Plotting
############

# plotting number of mults to convergence:
plt.plot(vary_mu, no_mults_sp2,     label="sp2")
plt.plot(vary_mu, no_mults_mw,      label="mw")
plt.plot(vary_mu, no_mults_sp8,     label="sp8")
plt.plot(vary_mu, no_mults_sp8acc,  label="sp8acc")

data_full_vary_mu = np.column_stack((vary_mu, no_mults_sp2, no_mults_mw, no_mults_sp8, no_mults_sp8acc))
np.savetxt('data/data_full_vary_mu.txt', data_full_vary_mu, fmt='%.18f', delimiter=' ', header='x y1 y2 y3 y4', comments='')

plt.legend()
plt.xlim(left=0.49)
plt.ylim(bottom=-0.01)
plt.title(f"no. mults to convergence, gapsize={gapsize}")
plt.xlabel(xlabel="gap location")
plt.ylabel(ylabel="matmuls")
plt.savefig("no_mults_to_convergence.pdf", bbox_inches='tight')
