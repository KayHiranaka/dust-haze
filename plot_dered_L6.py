############
# plot chi squared improvement vs delta(J-K)
# compare original and dereddended spectra
#############

import numpy as np
import random
import asciitable
import matplotlib.pyplot as plt
import cPickle
from scipy.interpolate import griddata
from scipy.stats import nanmean
import locale
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
locale.setlocale(locale.LC_ALL, '')

a_emc = []
b_emc = []
n_emc = []
c_emc = []
ao_emc = []
bo_emc = []
no_emc = []
co_emc = []

# Define chi squared function
def chisq(ydata,ymod,sd):
    ymod = ymod[np.where(np.isfinite(ymod))]
    ydata = ydata[np.where(np.isfinite(ymod))]
    sd = sd[np.where(np.isfinite(ymod))]
    ymod = ymod[np.where(np.isfinite(ydata))]
    ydata = ydata[np.where(np.isfinite(ydata))]
    sd = sd[np.where(np.isfinite(ydata))]
    sd[np.where(np.isnan(sd))] = 1.0
    chisq=np.sum( ((ydata-ymod)/sd)**2)
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
    Q_grid.append(mie['Qex_averaged'][i*184:i*184+184])
    Qa = []
    for k in range(len(Q_grid[i])):
        Qa.append(Q_grid[i][k]*np.pi)
    Qchi.append(Qa)


#### Input field L dwarf name ####            
field_name = []
field_name.append('2M0345_L0')
field_name.append('2M2130_L1_kc')
field_name.append('Kelu-1_L2')
field_name.append('u11291_1506+1321_050323')
field_name.append('2M2158-15_L4')
field_name.append('2M1507')
field_name.append('2MASSIJ1010148-040649')

# low-g Ls
lg = asciitable.read('Files/lowgred/lowgred.txt')
lgname = lg['filename']
lgst = lg['sptype']
lgcol = lg['J-Ks']

# field red Ls
data = asciitable.read('Files/fieldred/fieldred_refined.txt')
desi = data['designation']
spty = data['spt']
color =data['J-Ks']
name = data['name']

# L6s
data6 = asciitable.read('Files/redL6.txt')
desi6 = data6['filename']
spty6 = data6['spt']
color6 = data6['J-Ks']

dJK = []
JKav = [1.33, 1.33, 1.67, 1.63, 1.86, 1.46, 1.89]
sh = ['o', 's', 'v', 'd', 'p', '8']
marky = ['mo', 'm^', 'md', 'mp', 'mh', 'm8']
markf = ['go', 'g^', 'gd', 'gp', 'gh', 'g8']
shape = []
JK0 = []
JK1 = []
JK2 = []
JK3 = []
JK4 = []
JK5 = []
chi_lg = []
chi_lg_der = []
dchi = []
chirat = []

plota = plt.subplot(111)
#time = '20150612-163558'
time = '20151229-125357'
for ll in range(0, len(desi6)):
    red = desi6[ll]
    spt = str(spty6[ll])
    stype = int(spt[1])
    col = color6[ll]
    if col-JKav[stype] < 0.1:
        continue

    samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(red,stype,time), 'rb'))
    if red == '2MASSIJ0103320+193536':
        datar = asciitable.read('Files/lowgred/%s.txt' %(red))
    else:
        datar = asciitable.read('Files/fieldred/%s.txt' %(red))
    war = datar['col1']
    flur = datar['col2']
    fer = datar['col3']
    
    field = field_name[stype]
    colav = JKav[stype]
    dJK.append(col-colav)
    print  red, stype

    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    a_emc.append(a_mcmc[0])
    b_emc.append(b_mcmc[0])

# field object spectrum
    if stype < 6: 
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

# Calculate chi squared 
    chi2_lg = chisq(fluf_norm,flur_norm,fef)
    if np.isnan(chi2_lg):
        chi2_lg = 0

    locale.format("%d", chi2_lg, grouping=True)
    chi_lg.append(chi2_lg)

# Deredden low-g object & calculate chi squared
    Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
    Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
    Qint = np.interp(waf, w, Qint)    
    flur_der = flur * np.exp(Qint)
    flur_der_norm = flur_der/nanmean(flur_der)
    chi2_lg_der = chisq(fluf_norm,flur_der_norm,fef)
    dofh = len(waf)-5
    if np.isnan(chi2_lg_der):
        chi2_lg_der = 0

    chi_lg_der.append(chi2_lg_der)

    plt.plot(waf, fluf_norm, 'k', label = 'L%s standard'%(stype))
    plt.plot(waf, flur_norm, 'r', label='low-g, $\chi^2$=%3s, DOF=%s'%(format(float('%.3g'%(chi2_lg)),'n'),len(waf)))
    plt.plot(waf, flur_der_norm, 'g', label='dereddened low-g, $\chi^2$=%3s, DOF=%s'%(format(float('%.3g'%(chi2_lg_der)),'n'),dofh))
    plt.xlim(0.8, 2.5)
    plt.xlabel('wavelength [${\mu}m]$')
    plt.ylabel('Normalized Flux')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
    leg = plt.legend(fancybox=True, loc=1, prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    plt.savefig('Plots/dered_%s_L%s_%s.pdf'%(red,stype, timestr))
    plt.clf()
#    quit()
    dchi.append(chi2_lg - chi2_lg_der)
    chirat.append(chi2_lg/chi2_lg_der)
quit()
plota.scatter(dJK, dchi, color='0.5')
#plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
#plota.plot(dJK, chi_lg_der, 'g.')
#plota.set_xlim(-0.5, 3)
#plota.set_ylim(0, 0.5)
plota.axes.grid(True, linestyle = '-', color = '0.75')
#plota.set_title('$\Delta\chi^2$ vs $\Delta$(J - K)', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
plota.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
plota.set_ylabel('$\Delta\chi^2$',fontsize=14)
#plota.plot(dcol, dchi_f, 'mo')
plt.savefig('Plots/dchisq_djk_%s.pdf'%timestr)
plt.clf()
#plt.savefig('Plots/dchisq_lg_djk_%s.pdf'%timestr)
#plt.clf()
plotb = plt.subplot(111)
#plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
plotb.scatter(dJK, chirat, color='0.5')
#plota.plot(dJK, chi_lg_der, 'g.')
#plota.set_xlim(-0.5, 3)
#plota.set_ylim(0, 0.5)
plotb.axes.grid(True, linestyle = '-', color = '0.75')
#plotb.set_title('$\chi^2$ ratio vs $\Delta$(J - K)', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
plotb.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
plotb.set_ylabel('$\chi^2_{before}/\chi^2_{after}$',fontsize=14)
#plotb = plt.subplot(111)
#plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
#plotb.plot(dcol, chirat_f, 'mo')
#plota.plot(dJK, chi_lg_der, 'g.')
#plota.set_xlim(-0.5, 3)
#plota.set_ylim(0, 0.5)
#plotb.axes.grid(True, linestyle = '-', color = '0.75')
#plotb.set_title('$\chi^2$ ratio vs $\Delta$(J - K)', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
#plotb.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
#plotb.set_ylabel('$\chi^2$ ratio',fontsize=14)
plt.savefig('Plots/chisqratio_djk_%s.pdf'%timestr)
plt.clf()


