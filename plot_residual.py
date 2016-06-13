#################
# Compare MCMC model fits with different conditions (particle size distributions etc)
# Plot best fit models and residuals
# Report chi2 and degrees of freedom
#################

import numpy as np
import matplotlib.pyplot as plt
import cPickle
from scipy.interpolate import griddata
import asciitable
import time
from scipy.stats import nanmean
timestr = time.strftime("%Y%m%d-%H%M%S")

# Define chi squared function                                                                                                 
def chisq(ydata,ymod,sd,s=None):
    ymod = ymod[np.where(np.isfinite(ymod))]
    ydata = ydata[np.where(np.isfinite(ymod))]
    sd = sd[np.where(np.isfinite(ymod))]
    ymod = ymod[np.where(np.isfinite(ydata))]
    ydata = ydata[np.where(np.isfinite(ydata))]
    sd = sd[np.where(np.isfinite(ydata))]
    sd[np.where(np.isnan(sd))] = 1.0
    if s==None:
        chisq=np.sum( ((ydata-ymod)/sd)**2)
    else:
        chisq=np.sum((ydata-ymod)**2/(sd**2+s**2))
    return chisq

#### Input field standard L dwarf name ####                                                                                       
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

stype = []

# extinction coeffs with power law
lenp = 2
power = cPickle.load(open('Files/powermieff.pkl'))
wp = np.asarray(power['wavelength'][0:121]) # len 121                                                                  
p = np.asarray(power['index'])[np.arange(0,121*lenp,121)]
Qp = []
#for i in range(0,lenp):

#    Qp.append(power['Qpower'][i*121:i*121+121])

Qp3 = np.asarray(power['Qpower'][0:121])
Qp35 = np.asarray(power['Qpower'][121:])
#r = 1.45
#r = 0.316
r = 0.037

lenm = 2
gauss = cPickle.load(open('Files/gaussmieff.pkl'))

wg = np.asarray(gauss['wave'][0:121])#len 121                                                                                 
mu = np.asarray(gauss['mu'])[np.arange(0,121*lenm,121)]
Qg05 = np.asarray(gauss['Q_gauss'][0:121])
Qg03 = np.asarray(gauss['Q_gauss'][121:])        

time = '20151025-190534'

# extinction coeffs with hansen 
mie = cPickle.load(open('pickles/hansenmieff.pkl'))
ab = 152
w = np.asarray(mie['wave'][0:184])
a = np.asarray(mie['a'])[np.arange(0,184*ab,184)]
b = np.asarray(mie['b'])[np.arange(0,184*ab,184)]
Q_grid = []
#Qchi = []
for i in range(ab):
#    loc = np.arange(i,184*ab,184) #rather than loop, just pick out the right locations                                 

    Q_grid.append(mie['Qex_averaged'][i*184:i*184+184])
 #   Qa = []

timeg = '20151106-140944'
timeh = '20150612-163558'
for l in range(0, 1):
    stype.append(0)
    red = newred[0]
    field = field_name[0]
                  
    if red == 'spex_prism_0141-4633_040905' or red == 'spex_prism_U10074_060821' or red == 'U20037_2M0045+1634_adam' or red == 'spex_prism_1726+1538_080713' or red == 'spex_prism_U20636_1551+0941_080713' or red == 'spex_prism_2002-0521_080908':
                      
        datar = asciitable.read('../Mie/ApplytoData/%s.dat' %(red))
        war = datar['wave']
        flur = datar['flux']
        fer = datar['ferr']
    else:
#        continue                                                                                                        
#newred                                                                                               
        datar = asciitable.read('../../ASTRO/Mie/Archive/newdata/%s.fits.txt' %(red))    
# Newreds  
        war = datar['col1']  
        flur = datar['col2'] 
        fer = datar['col3']

# standard object spectrum                                                                         
    dataf = asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(field), data_start=95, data_end=540)

    waf = dataf['wave']
    fluf = dataf['flux']
    fef = dataf['ferr']
    fef[np.where(fef==1.0)]=np.nan
# Normalize if comparing to red object                                                                                 
    fluf_norm = fluf/nanmean(fluf)
    fef = fef/nanmean(fluf)
#    fef[np.where(np.isnan(fef))]=1.0  
# Interpolate red spectra onto field spectra                                                                            
    flur = np.interp(waf, war, flur)
#    fer = np.interp(waf, war, fer)                                                                                      
# Normalize if comparing to field object                                                                                
    flur_norm = flur/nanmean(flur)

    data = asciitable.read('inputs/extinction_%s.unnorm'%(red))
    wave = data['wave']
    error_frat = data['flux ratio error']
    frat = data['normalized flux ratio']
    frat = np.delete(frat,np.where((wave>1.35) & (wave<1.45)))
    error_frat = np.delete(error_frat,np.where((wave>1.35) & (wave<1.45)))
    wave = np.delete(wave,np.where((wave>1.35) & (wave<1.45)))
    frat = np.delete(frat,np.where((wave>1.8) & (wave<2.0)))
    error_frat = np.delete(error_frat,np.where((wave>1.8) & (wave<2.0)))
    wave = np.delete(wave,np.where((wave>1.8) & (wave<2.0)))
    frat = np.interp(waf, wave, frat)
    error_frat = np.interp(waf,wave,error_frat)

    # Power law fits
    samples = cPickle.load(open('pickles/power3_spex_prism_0032-4405_U20026_L0.pkl'))
#    samples[:, 2] = np.exp(samples[:, 2])
    s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
#    s_emc.append(s_mcmc[0])
#    n_emc.append(n_mcmc[0])
#    c_emc.append(c_mcmc[0])
    Qint = np.interp(waf, wp, Qp3)
    Qint = n_mcmc[0]*np.pi*r**2*Qint + c_mcmc[0]
    flur_der = flur * np.exp(Qint)
#    flur_der = flur_der/nanmean(flur_der)
    
#    res = np.abs((fluf - flur_der)/fluf)*100
    res = ((frat - Qint)/frat)*100
    chip = chisq(frat,Qint,error_frat)
#,s_mcmc[0])
    dof = len(waf)-3

# gaussian fits
    samples = cPickle.load(open('pickles/gauss_%s_L%s_%s.pkl'%(red,stype[l],timeg)))
#    samples[:, 2] = np.exp(samples[:, 2])
    s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    Qintg = np.interp(waf, wg, Qg05)
    Qintg = n_mcmc[0]*np.pi*0.5**2*Qintg + c_mcmc[0]
    flur_der = flur * np.exp(Qintg)
#    flur_der = flur_der/nanmean(flur_der)
    
#    res = np.abs((fluf - flur_der)/fluf)*100
    resg = ((frat - Qintg)/frat)*100
    chig = chisq(frat,Qintg,error_frat)
#,s_mcmc[0])

# hansen
    samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,stype[l],timeh), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmch, b_mcmch, s_mcmch, n_mcmch, c_mcmch = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
#    a_emc.append(a_mcmc[0])
#    b_emc.append(b_mcmc[0])
    Qinth=griddata((a, b), Q_grid, (a_mcmch[0],b_mcmch[0]), method='linear')
    Qinth = n_mcmch[0]*np.pi*Qinth*a_mcmch[0]**2 + c_mcmch[0]
    Qinth = np.interp(waf, w, Qinth)
    flur_deh = flur * np.exp(Qinth)
 #   flur_deh = flur_deh/nanmean(flur_deh)
#    resh = np.abs((fluf - flur_deh)/fluf)*100
    resh = ((frat - Qinth)/frat)*100
    chih = chisq(frat,Qinth,error_frat)
#,s_mcmc[0])
    dofh = len(waf)-5

    fig1 = plt.figure(1)
    frame1=fig1.add_axes((.1,.3,.8,.6))
#    plt.plot(waf,fluf,c='k',label='field')
#    plt.plot(waf,flur_deh,c='g',label='de-reddened by hansen')
#    plt.plot(waf,flur_der,c='c',label='de-reddened by power law')
    plt.plot(waf,frat,c='k',label='Field')
    plt.plot(waf,Qinth,c='g',label='Hansen, $\chi^2$=%s, dof=%i'%(format(float('%.3g'%(chih)),'n'),dofh))
    plt.plot(waf,Qint,c='c',label='Power Law, $\chi^2$=%s, dof=%i'%(format(float('%.3g'%(chip)),'n'),dof))
    plt.plot(waf,Qintg,c='m',label='Gaussian, $\chi^2$=%s, dof=%i'%(format(float('%.3g'%(chig)),'n'),dof))
    leg = plt.legend(fancybox=True, loc=1, prop={'size':10})
    leg.get_frame().set_alpha(0.8)
#    plt.savefig('/Users/paigegiorla/C#    plt.legend()
#    plt.ylabel('Normalized Flux')
    plt.ylabel('Extinction')
    frame1.set_xticklabels([]) #Remove x-tic labels for the first frame
   # grid()
    
    frame2=fig1.add_axes((.1,.1,.8,.2))        
    plt.plot(waf,res, c='c',label='power law')
    plt.plot(waf,resh, c = 'g',label='hansen')
    plt.plot(waf,resg, c='m',label='gaussian')
    plt.xlabel('Wavelength ($\mu m$)')
    plt.ylabel('Percent Difference')
    leg = plt.legend(fancybox=True, loc=4, prop={'size':8})
    leg.get_frame().set_alpha(0.8)
#    plt.savefig('/Users/paigegiorla/Code/Python/BDNYC/interpolation/residuals/trimmed_{}_{}'.format(params[i][0],params[i][1])+'.png')
    plt.savefig('Plots/residual_ext_pow3_gauss5_hansen.pdf')
#    plt.savefig('Plots/residual_spec_pow_hansen.pdf')
    plt.show()
    plt.clf()
