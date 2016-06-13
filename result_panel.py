############
# plot results 
# one obj from each spec type
#############

import numpy as np
import random
import asciitable
import matplotlib.pyplot as plt
import matplotlib as mpl
#from scipy.stats import mode
#import triangle
#import matplotlib.mlab as mlab
#import emcee
import cPickle
from scipy.interpolate import griddata
#import matplotlib.axis as axis
from scipy.stats import nanmean
import time
timestr = time.strftime("%Y%m%d-%H%M%S")

a_emc = []
b_emc = []
n_emc = []
c_emc = []
s_emc = []
ao_emc = []
bo_emc = []
no_emc = []
co_emc = []

# Define chi squared function
def chisq(ydata,ymod,sd,s=None):
#    ydata = np.isfinite(ydata)
    ymod = ymod[np.where(np.isfinite(ymod))]
    ydata = ydata[np.where(np.isfinite(ymod))]
    sd = sd[np.where(np.isfinite(ymod))]
    ymod = ymod[np.where(np.isfinite(ydata))]
    ydata = ydata[np.where(np.isfinite(ydata))]
    sd = sd[np.where(np.isfinite(ydata))]
    sd[np.where(np.isnan(sd))] = 1.0
#    if np.isnan(sd):                                                                                                                                                           
#        chisq=np.sum(((ydata-ymod)**2))                                                                                                                                  
#    else:                                                                                                                                                             
    if s==None:
        chisq=np.sum( ((ydata-ymod)/sd)**2)
    else:
        chisq=np.sum((ydata-ymod)**2/(sd**2+s**2))
#    print ydata, ymod, sd, chisq
    return chisq
# Load pickled calculated Q_ex file for all wavelengths and reff&veff                                                                                                     
mie = cPickle.load(open('pickles/hansenmieff.pkl'))
ab = 152
w = np.asarray(mie['wave'][0:184])
a = np.asarray(mie['a'])[np.arange(0,184*ab,184)]
b = np.asarray(mie['b'])[np.arange(0,184*ab,184)]
Q_grid = []
Qchi = []
for i in range(ab):
#    loc = np.arange(i,184*ab,184) #rather than loop, just pick out the right locations                                                                                   
    Q_grid.append(mie['Qex_averaged'][i*184:i*184+184])
    Qa = []
    for k in range(len(Q_grid[i])):
        Qa.append(Q_grid[i][k]*np.pi)
    Qchi.append(Qa)


#### Input field L dwarf name ####            
field_name = []
field_name.append('2M0345_L0')
field_name.append('2M2130_L1_kc')
field_name.append('Kelu-1')
field_name.append('u11291_1506+1321_050323')
field_name.append('U12101_2158-1550_davy.fits')
field_name.append('2M1507')
field_name.append('2MASSIJ1010148-040649')
# proper field name
fro_name = []
fro_name.append('2M0345')
fro_name.append('2M2130')
fro_name.append('Kelu-1')
fro_name.append('2M1506')
fro_name.append('2M2158')
fro_name.append('2M1507')
fro_name.append('2M1010')
#### Input red L dwarf name ####                                                                                 
red_name = []
red_name.append('spex_prism_0141-4633_040905')
red_name.append('spex_prism_U10074_060821')
red_name.append('U20048')
red_name.append('spex_prism_1726+1538_080713')
red_name.append('U40005_2249+0044_katelyn')
red_name.append('J032642.32-210207.2')
red_name.append('2MASSJ21481628+4003593')
red_name.append('U20171_0355+1133')
# proper names                                                                                                   
pro_name=[]
pro_name.append('2M0141-4633')
pro_name.append('2M0117-3403')
pro_name.append('2M0055+0134')
pro_name.append('2M1726+1538')
pro_name.append('2M2249+0044')
pro_name.append('2M0326-2102')
pro_name.append('2148+4003')
pro_name.append('0355+1133')

stype = []

# field red Ls
data = asciitable.read('Files/fieldred/fieldred_refined.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']
color =data['J-Ks']
name = data['name']

chi_lg=[]
chi_lg_der=[]

label_size = 7
mpl.rcParams['ytick.labelsize'] = label_size
label_size = 7
mpl.rcParams['xtick.labelsize'] = label_size
combofig = plt.figure(figsize=(8,11))
#combofig2 = plt.figure(figsize=(8,11))
j=1
for l in range (0, 4):
    red = red_name[l]
#    JKavg = JKav[l]

# Red L dwarf                                                                                                                                                                 
    if l == 7:
        stype.append(6)
    else:
        stype.append(l)
    if l == 4 or l == 7: 
        datar = asciitable.read('../../ASTRO/Mie/Archive/newdata/%s.fits.txt' %(red)) 
        war = datar['col1']
        flur = datar['col2']
        fer = datar['col3']
    elif l == 6:
        datar = asciitable.read('Files/fieldred/%s.txt'%(red))
        war = datar['col1']
        flur = datar['col2']
        fer = datar['col3']
    elif l == 5:
        datar = asciitable.read('Files/fieldred/%s.txt'%(red))
        war = datar['wavelength']
        flur = datar['flux']
        fer = datar['ferr']
    else:
        datar = asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(red))
        war = datar['wave']
        flur = datar['flux']
        fer = datar['ferr']
#    ll = l - 4
    ll = l
    print stype[ll], red
    field = field_name[stype[ll]]

    if l == 2:
        time = '20160111-233120'
    elif l == 4:
        time = '20160115-134851'
    elif l == 6:
         time = '20151229-125357'
    elif l == 7:
        time = '20151215-161107'
    else:
        time = '20150612-163558'
    if l == 5 or l == 6:
        samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(red,stype[ll],time), 'rb'))
    else:
        samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,stype[ll],time), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    a_emc.append(a_mcmc[0])
    b_emc.append(b_mcmc[0])
    s_emc.append(s_mcmc[0])
# field object spectrum
    if l not in (2,4,6,7):
        dataf = asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(field), data_start=95, data_end=540)
        waf = dataf['wave']
        fluf = dataf['flux']
        fef = dataf['ferr']
    else:
        dataf = asciitable.read('Files/%s.txt' %(field), data_start=95, data_end=540)
        waf = dataf['col1']
        fluf = dataf['col2']
        fef = dataf['col3']

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
#    fer = fer/nanmean(flur)
#    print fef, flur_norm
# Calculate chi squared 
    chi2_lg = chisq(fluf_norm,flur_norm,fef,s_mcmc[0])
#,len(waf)-1)
    if np.isnan(chi2_lg):
        chi2_lg = 0
#    chi2_lg[np.where(np.isnan(chi2_lg))]=0
    
    chi_lg.append(chi2_lg)

#,len(waf)-1))

# Deredden low-g object & calculate chi squared
#    for j in range(0, ab):
#        Qchi2 = np.interp(waf, w, Qchi[j])
    Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
    Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
    Qint = np.interp(waf, w, Qint)    
    flur_der = flur * np.exp(Qint)
    flur_der = flur_der/nanmean(flur_der)
    chi2_lg_der = chisq(fluf_norm,flur_der,fef,s_mcmc[0])
#,len(waf)-5)
    if np.isnan(chi2_lg_der):
        chi2_lg_der = 0

    chi_lg_der.append(chi2_lg_der)

    ext = asciitable.read('inputs/extinction_%s.unnorm'%(red))
    wavee = ext['wave']
    frate = ext['normalized flux ratio']

#    if l in (0,1,2,3):
#    ll = l - 4
    axe = combofig.add_subplot(4,2,2*(ll+1))
    axer = combofig.add_subplot(4,2,2*(ll+1)-1)
#    else:
#    axe = combofig.add_subplot(len(red_name),2,2*(l-3))
#    axer = combofig.add_subplot(len(red_name),2,2*(l-3)-1)
    axe.locator_params(axis='y',nbins=4)
    axe.plot(waf,Qint,'g',linewidth=0.8)
    dif = max(Qint)-min(Qint)
    axe.text(1.9, max(Qint)+0.1*dif,'a=%1.2f(+%1.2f,-%1.2f) $\mu$m'%(a_mcmc[0],a_mcmc[1],a_mcmc[2]),fontsize=5)
    axe.text(1.9, max(Qint)-0.1*dif,'b=%1.2f(+%1.2f,-%1.2f)'%(b_mcmc[0],b_mcmc[1],b_mcmc[2]),fontsize=5)
    axe.text(1.9, max(Qint)-0.3*dif,'N=%1.2f(+%1.2f,-%1.2f) $10^8/cm^2$'%(n_mcmc[0],n_mcmc[1],n_mcmc[2]),fontsize=5)
    axe.plot(wavee, frate, 'k', linewidth=0.7)
    axe.set_xlim(0.9, 2.5)
    if l in (3,7):
        axe.set_xlabel('wavelength [$\mu$m]')
    if l in (2,6):
        axe.set_ylabel('Observed Reddening')
    axer.plot(waf, fluf_norm, 'k', label = 'L%s standard %s'%(stype[ll],fro_name[stype[ll]]),linewidth=0.7)

    if l == 7:
        axer.plot(waf, flur_norm, 'r', label='0355-type low-g %s, $\chi^2$=%s'%(pro_name[l],float('%.3g'%(chi_lg[ll]))),linewidth=0.7)
    else:
        axer.plot(waf, flur_norm, 'r', label='L%s low-g %s, $\chi^2$=%s'%(stype[ll],pro_name[l],float('%.3g'%(chi_lg[ll]))),linewidth=0.7)
    axer.plot(waf, flur_der, 'g', label='dereddened %s, $\chi^2$=%s'%(pro_name[l],float('%.3g'%(chi_lg_der[ll]))),linewidth=0.7)
    axer.locator_params(axis='y',nbins=4)
#    axe.text(2.0,0.9*max(fluf_norm),'L%s:%s'%(l,pro_name[l]),fontsize=7)
    axer.set_xlim(0.9, 2.5)
    if l in(2,6):
        axer.set_ylabel('Normalized Flux')
    if l in (3,7):
        axer.set_xlabel('wavelength [$\mu$m]')
       
    axer.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
    leg = plt.legend(fancybox=True, loc=1, prop={'size':5})
    leg.get_frame().set_alpha(0.8)
#    plt.title('%s, L%s'%(red, stype[l]))
#combofig1.savefig('Plots/combo_s1.pdf')
#combofig2.savefig('Plots/combo_s2.pdf')
plt.savefig('Plots/combo_s1.pdf')
#plt.show()
plt.clf()
#    quit()
 #   dchi.append(chi2_lg - chi2_lg_der)
 #   chirat.append(chi2_lg/chi2_lg_der)

