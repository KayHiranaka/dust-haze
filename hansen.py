##############
# Integrate Qext over hansen particle size distribution
##############

import numpy as np
import asciitable
from scipy.integrate import quad
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import cPickle
import time
timestr = time.strftime("%Y%m%d-%H%M%S")

#mie = ast.read('Files/Mg2SiO4_fine.miefftable', data_start = 1)
#wave = mie['col1']
#rad = mie['col2']
#Qex = mie['col4']
#obj = {'wave':wave, 'rad':rad, 'Qex':Qex}
mie = cPickle.load(open('pickles/miefftable.pkl', 'rb'))
wave = np.asarray(mie['wave'])
rad = np.asarray(mie['rad'])
Qex = np.asarray(mie['Qex'])
#print max(rad)
#quit()
#aa = [2.0, 1.0, 2.0]
#ww = np.arange(0.8, 5.0, 0.01)

#for j in range(0, len(aa)):
#    r = aa[j]
#    print r
#    Qe = Qex[np.where(np.abs(rad-r)==np.amin(np.abs(rad-r)))]
#    w = wave[np.where(np.abs(rad-r)==np.amin(np.abs(rad-r)))]
#    ra = rad[np.where(np.abs(rad-r)==np.amin(np.abs(rad-r)))]
#    print ra
#    plt.plot(w, Qe, 'g')                                                                                                                                    
#    plt.xlim(0.9, 5)                                                                                                                                   
#    plt.title('%sum'%r)
#    plt.xlabel('wavelength [um]')
#    plt.ylabel('forsterite extinction coefficient')
#    plt.show()       
#    plt.savefig('Plots/%sum.pdf'%r)
#    plt.clf()

# Define Hansen distribution
def hansendist(r, a, b):
    return r**((1-3*b)/b)*np.exp(-r/a/b)

# Hansen*scattering area
def phansen(r, a, b):
    return np.pi*r**2*r**((1-3*b)/b)*np.exp(-r/a/b)

# Hansen*extinction coefficient
def qhansen(r, a, b, w):
#    Q = Qex[np.where((np.abs(wave-ww)<1e-4) & (np.abs(rad-r)<1e-4))]
    Qw = Qex[np.where(wave==w)]
    rw = rad[np.where(wave==w)]
    Q = np.interp(r, rw, Qw)
    return np.pi*r**2*Q*r**((1-3*b)/b)*np.exp(-r/a/b)

#aa = np.arange(0.05, 0.45, 0.05)
#bb = np.arange(0.1, 1.05, 0.2)
aa = [0.25]
bb = [0.5, 0.1, 1.0]
col = ['b', 'c', 'g', 'm', '0.5']
sty = ['-', '--', '-.', ':']
rr = np.arange(0.01, 10, 0.01)
aaa = []
bbb = []
www = []
ww = wave[np.arange(0,len(wave),1000)]
ww = ww[np.where(ww>1)]
ww = ww[np.where(ww<12)]
#ww = wave[np.arange(0,len(wave),1000)]

for u in range(0, len(aa)):
    a = aa[u]
    for v in range(0, len(bb)):
        b = bb[v]
        qnorm = []
        
        for w in range(0, len(ww)):
            aaa.append(a)
            bbb.append(b)
#            wwww.append(ww[w])
#            www.append(ww[w])
            print a, b, ww[w]
            sq = quad(qhansen, 0.01, 10, args = (a,b,ww[w]))
            sp = quad(phansen, 0.01, 10, args = (a,b))
            qnorm.append(sq[0]/sp[0])
#        plt.plot(ww, qnorm, 'g')
        c = col[u]
        s = sty[v]
        plt.plot(ww, qnorm, color = c, ls = s, label = 'a=%s [$\mu m$], b=%s'%(a,b))
plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
leg = plt.legend(fancybox=True, loc=1, prop={'size':7})
leg.get_frame().set_alpha(0.5)
#plt.xlim(1, 12)
#plt.ylim(1, 4)
#plt.title('extinction at long wavelengths')                                                                                                                                
#plt.title('Model Grid')                                          
plt.xlabel('wavelength [$\mu m$]')                                                                                                                      
plt.ylabel('forsterite extinction coefficient')        
plt.savefig('Plots/Qex_long_%s.pdf'%(timestr))
#plt.savefig('Plots/Qex_small.pdf')
plt.clf()

quit()
obj2 = {'a':aaa, 'b':bbb, 'wave':www, 'Qex_averaged': qnorm}
#cPickle.dump(obj2, open('../../Documents/Kaystuff/pickles/hansenmieff.pkl', 'wb'))
#cPickle.dump(obj2, open('pickles/hansenmieff.pkl', 'wb'))
asciitable.write(obj2, 'Files/hansensmall.txt')
