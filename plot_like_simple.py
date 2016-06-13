############
# Plot loglikelihood as a function of one MCMC parameter, while keeping other parameters fixed
#############

from pylab import *
import numpy as np
import random
import asciitable
import matplotlib.pyplot as plt
from scipy.stats import mode
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import cPickle
from scipy.interpolate import griddata
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
# Load pickled calculated Q_ex file for all wavelengths and reff&veff
mie = cPickle.load(open('pickles/hansenmieff.pkl'))
ab = 152
w = np.asarray(mie['wave'][0:184])
a = np.asarray(mie['a'])[np.arange(0,184*ab,184)]
b = np.asarray(mie['b'])[np.arange(0,184*ab,184)]

data = asciitable.read('Files/fieldred/fieldred.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']
color = data['J-Ks']
name = data['name']
Q_grid = []
Qchi = []
for i in range(ab):
#    loc = np.arange(i,184*ab,184) #rather than loop, just pick out the right locations                                                             
    Q_grid.append(mie['Qex_averaged'][i*184:i*184+184])
        
    Qa = []
    for k in range(len(Q_grid[i])):
        Qa.append(Q_grid[i][k]*np.pi*a[i]**2)
    Qchi.append(Qa)

# Define log probablity, prior, and likelihood functios
def logprior(theta):
#    r, v, tol, n, c = theta
    r = theta
    if 0.05 < r < 0.4 and 0.1 < v < 1.0 and n > 0.0:
        return 0.0
    return -np.inf

#print type(a), type(b)
def loglike(theta,data):
#    r, v, tol, n, c = theta
    v = theta
    s = np.exp(tol)
    #interpolate Q_ext*area and calculate likelihood
    Qint=griddata((a,b), Q_grid, (r,v), method='linear')
    Qint = Qint[0]
#    print 'Qint is', type(Qint), shape(Qint)
    Qint = np.interp(wave, w, Qint)
    Qint = Qint * n * np.pi * r**2 + c # Q*scattering area + constant
    return -0.5 * np.sum(((frat-Qint)**2/((error_frat**2+s**2)))+np.log(2*np.pi*(error_frat**2+s**2)))

def logprob(theta, data):
    r = theta
    if not np.isfinite(logprior(theta)):
        return -np.inf
#    if not np.isfinite(loglike(tol,data)):
#        return -np.inf
    return logprior(theta) + loglike(theta, data)

red = 'u11538_050323'#L4 / 1615+4953       
data = asciitable.read('inputs/extinction_%s.unnorm'%(red))
wave = data['wave']
frat = data['normalized flux ratio']
error_frat = data['flux ratio error']
frat = np.delete(frat,np.where((wave>1.35) & (wave<1.45)))
error_frat = np.delete(error_frat,np.where((wave>1.35) & (wave<1.45)))
wave = np.delete(wave,np.where((wave>1.35) & (wave<1.45)))
frat = np.delete(frat,np.where((wave>1.8) & (wave<2.0)))
error_frat = np.delete(error_frat,np.where((wave>1.8) & (wave<2.0)))
wave = np.delete(wave,np.where((wave>1.8) & (wave<2.0)))

aa = [0.3]
t = [-1.0, -3.0]
nn = [1.0, 3.0] #1e8 cm^-2 in um^-2
#c = 1.7
c = 3.0
bb = np.arange(0.1, 1.0, 1.e-3)
color = ['m', 'c']
style = ['-', '--']
alpha = [1.0, 0.5]

for k in range(0, len(aa)):
    r = aa[k]
    col = color[k]
    for j in range(0, len(nn)):
        n = nn[j]
        sty = style[j]
        for h in range(0, len(t)):
            tol = t[h]
            alp = alpha[h]
#            like = loglike(v, data)
            s = np.exp(tol)
            for v in bb:
                Qchain=griddata((a,b), Q_grid, (r,v), method='linear')
#        Qint = Qint[0]                                                                                                                                 
    #    print 'Qint is', type(Qint), shape(Qint)                                                                                                       
                Qchain = np.interp(wave, w, Qchain)
                Qchain = Qchain * n * np.pi * r**2 + c # Q*scattering area + constant                                                                           
                llike = -0.5 * np.sum(((frat-Qchain)**2/((error_frat**2+s**2)))+np.log(2*np.pi*(error_frat**2+s**2)))
                plt.plot(v, llike, 'b.')
#                plt.plot(v, llike, color=col, linestyle=sty, alpha=alp, label='a = %s[$\mu$m], N = %s [$10^8 cm^{-2}$], log(s) = %s'%(r,n,tol))
#plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
#leg = plt.legend(fancybox=True, loc=1, prop={'size':8})
#leg.get_frame().set_alpha(0.5)
            plt.title('a=%s, N=%s, log(s)=%s, c=%s'%(r,n,tol,c))
            plt.xlabel('b')
            plt.ylabel('loglikelihood')
            plt.savefig('Plots/like_%s_a%s_n%s_t%s_c%s.pdf'%(red,r,n,tol,c))
            plt.clf()
#    plt.show()
quit()
