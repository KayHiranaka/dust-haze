##############
# Integrate Qext over power law particle size distribution
##############

import numpy as np
import asciitable
from scipy.integrate import quad
from scipy.integrate import trapz
import matplotlib.pyplot as plt
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

# Define power law distribution
def power(r, p):
    return r**(-p)

# n(r)*scattering area
def ppower(r, p):
    return np.pi*r**2*r**(-p)

# n(r)*extinction coefficient*area
def qpower(r, p, w):
#    Q = Qex[np.where((np.abs(wave-ww)<1e-4) & (np.abs(rad-r)<1e-4))]
    Qw = Qex[np.where(wave==w)]
    rw = rad[np.where(wave==w)]
    Q = np.interp(r, rw, Qw)
    return np.pi*r**2*Q*r**(-p)

def meanpow(r, p):
#    return np.pi*r**2*r**(-p+1)
    return r**(-p+2)
#aa = np.arange(0.05, 0.45, 0.05)
#bb = np.arange(0.1, 1.05, 0.2)
#aa = [0.2, 0.5]
#bb = [0.4, 0.4]
col = ['b', 'c', 'g', 'm']
sty = ['--', '-']
rr = np.arange(0.01, 10, 0.01)
www = []
ww = wave[np.arange(0,len(wave),1000)]
ww = ww[np.where(ww>0.8)]
ww = ww[np.where(ww<6)]
#ww = wave[np.arange(0,len(wave),1000)]
pp = [3, 3.5]
pppp = []
#meanr = []
qnorm = []
for ppp in range(0, len(pp)):
    p = pp[ppp]
    meanr = quad(meanpow, 0.01, 10, args = (p))
#    sp = quad(ppower, 0.01, 10, args = (p))
    sp = quad(power, 0.01, 10, args = (p))
    mr = meanr[0] / sp[0]
    print np.sqrt(mr)
    if ppp == len(pp)-1:
        break
    else:
        continue
    
    for w in range(0, len(ww)):
        www.append(ww[w])
        pppp.append(p)
        #    print a, b, ww[w]
    #    sq = quad(qpower, 0.01, 10, args = (p, ww[w]))
    #    sp = quad(ppower, 0.01, 10, args = (p))
    #    qnorm.append(sq[0]/sp[0])

#    plt.plot(ww, qnorm, color = col[ppp],label ='p = -%s'%(p))
#        c = col[u]
#        s = sty[v]
#        plt.plot(ww, qnorm, color = c, ls = s, label = 'a=%s [$\mu m$], b=%s'%(a,b))
#print len(pppp), len(www), len(qnorm)
#quit()
#power = {'index': pppp, 'wavelength': www, 'Qpower': qnorm}

#cPickle.dump(power, open('Files/powermieff03.pkl', 'wb'))
quit()

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
leg = plt.legend(fancybox=True, loc=1)
leg.get_frame().set_alpha(0.5)
plt.xlim(0.8, 2.5)
plt.ylim(0, 3.5)
plt.title('$n(r) = r^p$')                                                                                                                                
plt.xlabel('wavelength [$\mu m$]')                                                                                                                      
plt.ylabel('forsterite extinction coefficient')        
plt.savefig('Plots/Qex_power.pdf')
#plt.savefig('Plots/Qex_small.pdf')
plt.clf()
        
quit()
obj2 = {'a':aaa, 'b':bbb, 'wave':www, 'Qex_averaged': qnorm}
#cPickle.dump(obj2, open('../../Documents/Kaystuff/pickles/hansenmieff.pkl', 'wb'))
#cPickle.dump(obj2, open('pickles/hansenmieff.pkl', 'wb'))
asciitable.write(obj2, 'Files/hansensmall.txt')
