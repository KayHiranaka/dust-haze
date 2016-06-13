##############
# Integrate Qext over gaussian particle size distribution
##############

import numpy as np
import asciitable
from scipy.integrate import quad
from scipy.integrate import trapz
import matplotlib.pyplot as plt
from scipy.stats import norm
import cPickle
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
mie = cPickle.load(open('pickles/miefftable.pkl', 'rb'))
wave = np.asarray(mie['wave'])
rad = np.asarray(mie['rad'])
Qex = np.asarray(mie['Qex'])


# Define Hansen distribution
#def hansendist(r, a, b):
#    return r**((1-3*b)/b)*np.exp(-r/a/b)

# n(r)*scattering area
def pgauss(r, mu, sigma):
    nr = norm.pdf(r, mu, sigma)
    return np.pi*r**2*nr

# n(r)*extinction coefficient*area
def qgauss(r, mu, sigma, w):
#    Q = Qex[np.where((np.abs(wave-ww)<1e-4) & (np.abs(rad-r)<1e-4))]
    Qw = Qex[np.where(wave==w)]
    rw = rad[np.where(wave==w)]
    Q = np.interp(r, rw, Qw)
    nr = norm.pdf(r, mu, sigma)
    return np.pi*r**2*Q*nr

#aa = np.arange(0.05, 0.45, 0.05)
#bb = np.arange(0.1, 1.05, 0.2)
#aa = [0.2, 0.5]
#bb = [0.4, 0.4]
col = ['m', 'c', 'g', 'm']
sty = ['--', '-']
rr = np.arange(0.01, 10, 0.01)
www = []
ww = wave[np.arange(0,len(wave),1000)]
ww = ww[np.where(ww>0.8)]
ww = ww[np.where(ww<6)]
#ww = wave[np.arange(0,len(wave),1000)]
#pp = [3, 3.5]
m = [0.5, 0.3]
x = np.logspace(-2, 1, 200)
mm = []
qnorm = []
for ppp in range(0, len(m)):
    # Gaussian distributions                                                                                                                                
    mu = m[ppp]
    sigma = 0.1 * mu / np.sqrt(2.)
#snr = trapz(norm.pdf(x, mu, sigma), x)
#    p = pp[ppp]
    print mu
    for w in range(0, len(ww)):
        mm.append(mu)
        print ww[w]
        www.append(ww[w])
#    print a, b, ww[w]
        sq = quad(qgauss, 0.01, 10, args = (mu,sigma, ww[w]))
        sp = quad(pgauss, 0.01, 10, args = (mu,sigma))
        qnorm.append(sq[0]/sp[0])
#    plt.plot(ww, qnorm, color = 'm', linestyle = sty[ppp],label ='$\mu$ = %s $\mu$m'%(mu))
#        c = col[u]
#        s = sty[v]
#        plt.plot(ww, qnorm, color = c, ls = s, label = 'a=%s [$\mu m$], b=%s'%(a,b))

ext_gauss = {'mu': mm, 'wave':www, 'Q_gauss': qnorm}
cPickle.dump(ext_gauss, open('Files/gaussmieff.pkl','wb'))
quit()

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
leg = plt.legend(fancybox=True, loc=1)
leg.get_frame().set_alpha(0.5)
plt.xlim(0.8, 2.5)
plt.ylim(0, 3.5)
plt.title('Gaussian')                                                                                                                                
plt.xlabel('wavelength [$\mu m$]')                                                                                                                      
plt.ylabel('forsterite extinction coefficient')        
plt.savefig('Plots/Qex_gauss.pdf')
#plt.savefig('Plots/Qex_small.pdf')
plt.clf()
        
quit()
obj2 = {'a':aaa, 'b':bbb, 'wave':www, 'Qex_averaged': qnorm}
#cPickle.dump(obj2, open('../../Documents/Kaystuff/pickles/hansenmieff.pkl', 'wb'))
#cPickle.dump(obj2, open('pickles/hansenmieff.pkl', 'wb'))
asciitable.write(obj2, 'Files/hansensmall.txt')
