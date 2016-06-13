############
# MCMC parameter N (column denisity) and constant offset (C) with s
# Start with minimum chi squared & restart with maximum probability 
# Add red field L dwarfs (alphas&betas)
# Compare with Gaussian and power law grain size distributions
#############

from pylab import *
import numpy as np
import random
import asciitable
import matplotlib.pyplot as plt
from scipy.stats import mode
import triangle
import matplotlib.mlab as mlab
from matplotlib.backends.backend_pdf import PdfPages
import emcee
import cPickle
from scipy.interpolate import griddata
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
# Load pickled calculated Q_ex file for all wavelengths and reff&veff
mie = cPickle.load(open('pickles/hansenmieff.pkl'))
#ab = 152
#w = np.asarray(mie['wave'][0:184])
#a = np.asarray(mie['a'])[np.arange(0,184*ab,184)]
#b = np.asarray(mie['b'])[np.arange(0,184*ab,184)]
lenp = 2
lenm = 2
power = cPickle.load(open('Files/powermieff.pkl'))
wp = np.asarray(power['wavelength'][0:121]) # len 121
p = np.asarray(power['index'])[np.arange(0,121*lenp,121)]
#print len(wp),power['Qpower']
#p = np.asarray(power['index'])
#Qp = power['Qpower']
gauss = cPickle.load(open('Files/gaussmieff.pkl'))
wg = np.asarray(gauss['wave'][0:121])#len 121
mu = np.asarray(gauss['mu'])[np.arange(0,121*lenm,121)]
#Qg = gauss['Q_gauss']
#print len(Qg), len(gauss['mu'])
#quit()

data = asciitable.read('Files/fieldred/fieldred.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']
color = data['J-Ks']
name = data['name']
Qp = []
Qchip = []
for i in range(0,lenp):
#    loc = np.arange(i,121*lenp,121) #rather than loop, just pick out the right locations                                                             
    Qp.append(power['Qpower'][i*121:i*121+121]) 
    Qap = []
    for k in range(0,len(Qp[i])):
        Qap.append(Qp[i][k]*np.pi*10**2)
    Qchip.append(Qap)
Qg = []
Qchig = []
for i in range(lenm):
    loc = np.arange(i,121*lenm,121) #rather than loop, just pick out the right locations                                                             
    Qg.append(gauss['Q_gauss'][i*121:i*121+121]) 
        
    Qag = []
    for k in range(len(Qg[i])):
        Qag.append(Qg[i][k]*np.pi*mu[i]**2)
    Qchig.append(Qag)

# Define log probablity, prior, and likelihood functios
def logprior(theta):
    tol, n, c = theta
#    r = 10
#    if 0.05 < r < 0.4 and 0.1 < v < 1.0 and 
    if n > 0.0:
        return 0.0
    return -np.inf

# set initial values for parameters
t_init = 0
n_init = 0 #1e8 cm^-2 in um^-2
#t = -2.0
#n = 2 #1e8 cm^-2 in um^-2
#c_init = 2

def chisq(ydata,ymod,sd):
#    if sd==None:
#        chisq=np.sum(((ydata-ymod)**2))
#    else:
    chisq=np.sum( ((ydata-ymod)/sd)**2 )

    return chisq

ap_emc = []
bp_emc = []
Np_emc = []
Cp_emc = []
ag_emc = []
bg_emc = []
Ng_emc = []
Cg_emc = []
#print type(a), type(b)
#r = 1.45 #weighted mean effective radius for r^-3
#r = 0.316 #weighted mean effective radius for r^-3.5
r = 0.037 # RMS radius for r^-3
#r=0.022# RMS radius for r^-3.5
#r = 0.02 mean radius for r^-3
#r = 0.017 mean radius for r^-3.5
def loglike(theta,data):
    tol, n, c = theta
    s = np.exp(tol)
    #interpolate Q_ext*area and calculate likelihood
#    Qint=griddata((a,b), Q_grid, (r,v), method='linear')
#    Qint = Qint[0]
#    print 'Qint is', type(Qint), shape(Qint)
   # Qint = np.interp(pp, p, Qp[0])
    Qint = np.interp(wave, wp, Qp[0])
    Qint = Qint * n * np.pi * r**2 + c # Q*scattering area + constant
    return -0.5 * np.sum(((frat-Qint)**2/((error_frat**2+s**2)))+np.log(2*np.pi*(error_frat**2+s**2)))

def logprob(theta, data):
#    r, v, tol, n, c = theta
 #   r = theta
    if not np.isfinite(logprior(theta)):
        return -np.inf
#    if not np.isfinite(loglike(tol,data)):
#        return -np.inf
    return logprior(theta) + loglike(theta, data)

ndim, nwalkers = 3, 100

# Define a function to streamline the central tendency and error bar plots
def ctend_plot(point, ci, y, color, label):
    plt.plot(ci, [y, y], "-", color=color, linewidth=4, label=label)
    plt.plot(point, y, "o", color=color, markersize=10)

#### Input field L dwarf name ####            
field_name = []
field_name.append('2M0345_L0')
field_name.append('2M2130_L1_kc')
field_name.append('Kelu-1_L2')
field_name.append('u11291_1506+1321_050323')
field_name.append('2M2158-15_L4')
field_name.append('2M1507')

#### Input red L dwarf name (low-g)####                                                          
newred = []
newred.append('spex_prism_0032-4405_U20026')#L0                                                                                                                                                                                                                 
newred.append('spex_prism_0141-4633_040905') #L0                                                                                        
newred.append('U20098')#L0/J0210-3015                                                                                                                                                                                                                              
newred.append('spex_prism_U10141_060821') #L0                                                                                           
newred.append('spex_prism_0323-4631_U20157') #L0                                                                                        
newred.append('j2213-2136prismtc1')#L0                                                                                                                                                                                                                          
newred.append('spex_prism_2315+0617_U20993')#L0                                                                                                                                                                                                             
newred.append('U10381_JHK_2011dec08')#L1/0518-2756                                                                                                                                                                                                           
newred.append('U20048')#L2 / 0055+0134                                                                                                                                                                                                                    
newred.append('U10397')#L2 /0536-1920                                                                                                                                                                                                                  
newred.append('spex_prism_1726+1538_080713')#3                                                                                          
newred.append('U20001')#L3 / 0001+1535                                                                                                                                                                                                                            
newred.append('spex_prism_2208+2921_U40004')#L3                                                                                                                                                                                                                   
newred.append('2massj0126+1428spex')#L4                                                                                                                                                                                                                            
newred.append('U20198_2M0501-0010_dagny')#L4                                                                                                                                                                                                                     
newred.append('U20622')#L4 / 1538-1953                                                                                                  
  
newred.append('spex_prism_U20636_1551+0941_080713')#L4                                                                                  
newred.append('u11538_050323')#L4 / 1615+4953                                                                                                                                                                                                                              
newred.append('spex_prism_2206-4217_080714')#L4                                                                                                                                                                         
newred.append('U40005_2249+0044_katelyn')#L4                                                                                                                                                  
newred.append('U20171_0355+1133')#L5                                                                                                                                                                                                                                          
newred.append('U10372_JHK_2011dec08')#L5 / 0512-2949                                                                                    
newred.append('spex_prism_U10074_060821')#L1 / 0117-3403                                                                                
newred.append('U20037_2M0045+1634_adam')#L2                                                                                             
newred.append('spex_prism_2002-0521_080908')#5                                           

J = [14.776,14.832,15.066,15.799,15.389,15.376,15.861,15.262,16.436,15.768,15.669,15.522,15.797,17.1,14.982,15.934,16.319,16.789,15.555,16.587,14.05,15.463,15.178,13.059,15.316]
K = [13.269,13.097,13.5,14.035,13.702,13.756,14.065,13.615,14.438,13.854,13.659,13.71,14.148,15.3,12.963,14.004,14.31,14.306,13.609,14.358,11.526,13.285,13.49,11.37,13.417]
JK = []
JKav = [1.33, 1.33, 1.67, 1.63, 1.86, 1.46]
stype = []
l = 0
for k in range(0, 1):
    l = l + 1
    red = newred[k]
#    stype.append(4)                                                                                                                
#    field = field_name[5]
#    JKavg = JKav[5]
                                                                                                                      
    if k == 0 or k == 1 or k == 2 or k == 3 or k == 4 or k == 5 or k == 6:
        stype.append(0)
        field = field_name[0]
        JKavg = JKav[0]
    if k == 7 or k == 22:
        stype.append(1)
        field = field_name[1]
        JKavg = JKav[1]
    if k == 8 or k == 9 or k == 23:
        stype.append(2)
        field = field_name[2]
        JKavg = JKav[2]
    if k == 10 or k == 11 or k == 12:
        stype.append(3)
        field = field_name[3]
        JKavg = JKav[3]
    if k == 13 or k == 14 or k == 15 or k == 16 or k == 17 or k == 18 or k == 19:
        stype.append(4)
        field = field_name[4]
        JKavg = JKav[4]
    if k == 20 or k == 21 or k == 24:
        stype.append(5)                                                                                               
        field = field_name[5]
        JKavg = JKav[5]

    if J[k]-K[k]-JKavg < 0.1:
        continue

    print stype[l-1], red
#    data = asciitable.read('../../Documents/ASTRO/Mie/Miesults/normalized/extinction_%s.unnorm'%(red))
    data = asciitable.read('inputs/extinction_%s.unnorm'%(red))
    wave = data['wave']
    frat = data['normalized flux ratio']
    error_frat = data['flux ratio error']
    c_init = np.mean(frat)
#    c_init = 0

    frat = np.delete(frat,np.where((wave>1.35) & (wave<1.45)))
    error_frat = np.delete(error_frat,np.where((wave>1.35) & (wave<1.45)))
    wave = np.delete(wave,np.where((wave>1.35) & (wave<1.45))) 
    frat = np.delete(frat,np.where((wave>1.8) & (wave<2.0)))
    error_frat = np.delete(error_frat,np.where((wave>1.8) & (wave<2.0)))
    wave = np.delete(wave,np.where((wave>1.8) & (wave<2.0))) 

    chi2p = []
    chi2g = []
    for j in range(0, 2):
        Qchi2p = np.interp(wave, wp, Qchip[j])
        Qchi2p = n_init*Qchi2p + c_init
        chi2p.append(chisq(frat, Qchi2p, sd = error_frat))
    
#    r_init = a[np.where(chi2 == np.amin(chi2))]
#    v_init = b[np.where(chi2 == np.amin(chi2))]
    p_init = p[np.where(chi2p == np.amin(chi2p))]
    init = t_init, n_init, c_init # for random values
    print p_init, np.amin(chi2p)
    
    p0 = [init + np.random.rand(ndim)*1e-3 for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logprob, args=[data])
    print "set sampler"
    pos,probs,state = sampler.run_mcmc(p0, 2000)
    best_p = sampler.flatchain[sampler.flatlnprobability.argmax()]
    new_p0 = np.random.normal(best_p, 1e-3*np.abs(best_p), size=(nwalkers,ndim))
    pos,probs,state = sampler.run_mcmc(new_p0, 200)
    sampler.reset()
    pos, probs, state = sampler.run_mcmc(pos, 2000)
    cropchain = sampler.chain[:,:,:].reshape((-1,ndim))
#    random_samp = cropchain[np.random.randint(len(cropchain),size=200)]
#    param = 'r'
#    quit()
    output = open('pickles/power35w_%s_L%s.pkl'%(red,stype[l-1]), 'wb')
    cPickle.dump(cropchain, output)
    output.close()
    chainfig = plt.figure() 
    ax = chainfig.add_subplot(111)
    for tol, n, c in cropchain[np.random.randint(len(cropchain), size=100)]:
#        Qchain= np.interp(pp, Qp[0], p, method='linear')
        Qchain = np.interp(wave, wp, Qp[0])
        Qchain = Qchain * n * np.pi * r**2 + c
        ax.plot(wave,Qchain,'g',alpha=0.1)
        
    ax.plot(wave,frat,'k')
    ax.set_xlabel('$\lambda [\mu m]$')
    ax.set_ylabel('Observed Reddening, Extinction')
#    plt.title('%s, L%s'%(red,stype[l-1]))
#    .savefig('../../../Mie/Miesults/Plots/emceechain_area_s_L%s_%s_%s.pdf'%(stype[l],red, timestr)) 
#    chainfig.savefig('/Plots/emceechain_area_s_L%s_%s_%s.pdf'%(stype[k],red, timestr)) 
#    chainfig.savefig('../../../Mie/Miesults/Plots/emceechain_area_s_L%s_%s_%s.pdf'%(st,des, timestr)) 
    chainfig.savefig('Plots/emceechain_power_L%s_%s_%s.pdf'%(stype[l-1],red, timestr)) 
#    pdff.savefig(chainfig) 
    chainfig.clf()

    figure = triangle.corner(cropchain, labels=['log(s)', 'N [$10^8 cm^{-2}$]', 'C'], bins=40, quantiles=[0.16,0.5,0.84], verbose=True, show_titles=True)
#    figure = triangle.corner(cropchain, labels=['log(s)', 'N [$10^8 cm^{-2}$]', 'C'], bins=40, quantiles=[0.16,0.5,0.84], verbose=True, show_titles=True)
#        plt.tight_layout() 
#        plt.suptitle('%s, L%s'%(red,stype[l]), y=0.9, x=0.8)
#    plt.savefig('../../../Mie/Miesults/Plots/emcee_s_area_L%s_%s_%s.pdf'%(stype[l],red, timestr))
#    plt.savefig('../../../../../Dropbox/Codes/Plots/emcee_s_area_L%s_%s_%s.pdf'%(stype[l],red, timestr))
#    plt.suptitle('%s, L%s'%(red,stype[l-1]),y=0.9, x=0.8)
#    plt.savefig('../../../Mie/Miesults/Plots/emcee_s_area_L%s_%s_%s.pdf'%(st,des, timestr))
#        plt.savefig('Paper/triangle_L%s_%s.eps'%(stype[l],red))
    plt.subplots_adjust(hspace=0.1,wspace=0.1, bottom=0.15)
#    plt.tick_params(axis='both', which='major', labelsize=6)
#    plt.setp(figure.get_xticklabels(), fontsize=6)
#    plt.setp(figure.get_yticklabels(), fontsize=6)
#    plt.show()
    plt.savefig('Plots/emcee_power_L%s_%s_%s.pdf'%(stype[l-1],red, timestr))

#    quit()

#    pdft.savefig()
    plt.clf()
quit()
#    quit()

for l in range (0, len(desi)):
    if desi[l] == 'None' or desi[l] == 'J2224438-015852':
        continue
    if l > 0 and desi[l] == desi[l-1]:
        continue
    des = desi[l]
    spt = spty[l]
    nam = name[l]
    col = color[l]

    if spt == 10 or spt == 10.5:
        stype.append(0)
        st = 0
        field = field_name[0]
    if spt == 11 or spt == 11.5:
        st = 1
        stype.append(1)
        field = field_name[1]
    if spt == 12 or spt == 12.5:
        st = 2
        stype.append(2)
        field = field_name[2]
    if spt == 13 or spt == 13.5:
        st = 3
        stype.append(3)
        field = field_name[3]
    if spt == 14 or spt == 14.5:
        st = 4
        stype.append(4)
        field = field_name[4]
    if spt == 15 or spt == 15.5:

        st = 5
        stype.append(5)
        field = field_name[5]

#    if des == 'None':
#        des = nam
    if col-JKav[st] < 0.1:
        continue

    print st, des
#        data = asciitable.read('../../Documents/ASTRO/Mie/Miesults/normalized/extinction_%s.unnorm'%(red))
    data = asciitable.read('inputs/extinction_%s.unnorm'%(des))
    wave = data['wave']
    frat = data['normalized flux ratio']
    error_frat = data['flux ratio error']
    frat = np.delete(frat,np.where((wave>1.35) & (wave<1.45)))
    error_frat = np.delete(error_frat,np.where((wave>1.35) & (wave<1.45)))
    wave = np.delete(wave,np.where((wave>1.35) & (wave<1.45))) 
    frat = np.delete(frat,np.where((wave>1.8) & (wave<2.0)))
    error_frat = np.delete(error_frat,np.where((wave>1.8) & (wave<2.0)))
    wave = np.delete(wave,np.where((wave>1.8) & (wave<2.0))) 

    c_init = np.mean(frat)
    chi2 = []
    for j in range(0, ab):
        Qchi2 = np.interp(wave, w, Qchi[j])
        Qchi2 = n_init*Qchi2 + c_init
        chi2.append(chisq(frat, Qchi2, sd = error_frat))
        
    r_init = a[np.where(chi2 == np.amin(chi2))]
    v_init = b[np.where(chi2 == np.amin(chi2))]
#    print r_init[0], v_init[0], c_init
#    quit()
    init = r_init[0], v_init[0], t_init, n_init, c_init
    
    p0 = [init + np.random.rand(ndim)*1e-3 for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, logprob, args=[data])
    print "set sampler"
    pos,probs,state = sampler.run_mcmc(p0, 200)
    best_p = sampler.flatchain[sampler.flatlnprobability.argmax()]
    new_p0 = np.random.normal(best_p, 1e-3*np.abs(best_p), size=(nwalkers,ndim))
    pos,probs,state = sampler.run_mcmc(new_p0, 200)
    sampler.reset()
    pos, probs, state = sampler.run_mcmc(pos, 2000)
    cropchain = sampler.chain[:,:,:].reshape((-1,ndim))
    random_samp = cropchain[np.random.randint(len(cropchain),size=200)]
#        cPickle.dump(cropchain, open('../../Documents/ASTRO/codes/kay-repo/Python/%s_%s_area_cropchain.pkl'%(red,stype[l]), 'wb'))
    outputr = open('pickles/%s_%s_area_cropchain_%s.pkl'%(des,st,timestr), 'wb')
    cPickle.dump(cropchain, outputr)
    outputr.close()

    chainfig2 = plt.figure() 
    ax = chainfig2.add_subplot(111)
    for r, v, tol, n, c in cropchain[np.random.randint(len(cropchain), size=100)]:
        Qchain=griddata((a,b), Q_grid, (r,v), method='linear')
        Qchain = np.interp(wave, w, Qchain)
        Qchain = Qchain * n * np.pi * r**2 + c
        ax.plot(wave,Qchain,'g',alpha=0.1)
        
    ax.plot(wave,frat,'k')
    ax.set_xlabel('$\lambda [\mu m]$')
    ax.set_ylabel('Observed Reddening + Extinction')
#        plt.title('%s, L%s'%(red,stype[l]))
    plt.title('%s, L%s'%(des,st))
#    chainfig.savefig('../../../Mie/Miesults/Plots/emceechain_area_s_L%s_%s_%s.pdf'%(stype[l],red, timestr)) 
#    chainfig.savefig('../../../../../Dropbox/Codes/Plots/emceechain_area_s_L%s_%s_%s.pdf'%(stype[l],red, timestr)) 
#    chainfig.savefig('../../../Mie/Miesults/Plots/emceechain_area_s_L%s_%s_%s.pdf'%(st,des, timestr)) 
    chainfig2.savefig('Plots/emceechain_area_s_L%s_%s_%s.pdf'%(st,des, timestr)) 
    pdff.savefig(chainfig2) 
              
    figure = triangle.corner(cropchain, labels=['a [$\mu m$]', 'b', 'log(s)', 'N [$10^8 cm^{-2}$]', 'C'], bins=40, quantiles=[0.16,0.5,0.84], verbose=True, show_titles=True)
#        plt.tight_layout() 
#        plt.suptitle('%s, L%s'%(red,stype[l]), y=0.9, x=0.8)
#    plt.savefig('../../../Mie/Miesults/Plots/emcee_s_area_L%s_%s_%s.pdf'%(stype[l],red, timestr))
#    plt.savefig('../../../../../Dropbox/Codes/Plots/emcee_s_area_L%s_%s_%s.pdf'%(stype[l],red, timestr))
    plt.suptitle('%s, L%s'%(des,st),y=0.9, x=0.8)
#    plt.savefig('../../../Mie/Miesults/Plots/emcee_s_area_L%s_%s_%s.pdf'%(st,des, timestr))
#        plt.savefig('Paper/triangle_L%s_%s.eps'%(stype[l],red))
    plt.subplots_adjust(hspace=0.1,wspace=0.1, bottom=0.15)
#    plt.tick_params(axis='both', which='major', labelsize=6)
#    plt.setp(figure.get_xticklabels(), fontsize=6)
#    plt.setp(figure.get_yticklabels(), fontsize=6)
#    plt.show(
    plt.savefig('Plots/emcee_s_area_L%s_%s_%s.pdf'%(st,des, timestr))
    plt.close()
#    quit()

    pdft.savefig()
#plt.clf()
pdff.close()
pdft.close()
