############
# plot chi squared improvement vs delta(J-K)
#############

import numpy as np
import random
import asciitable
import matplotlib.pyplot as plt
#from scipy.stats import mode
#import triangle
#import matplotlib.mlab as mlab
#import emcee
import cPickle
from scipy.interpolate import griddata
#import matplotlib.axis as axis
from scipy.stats import nanmean
import locale
import time
timestr = time.strftime("%Y%m%d-%H%M%S")
locale.setlocale(locale.LC_ALL, '')
#print format(100000,'n')

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
    chisq=np.sum( ((ydata-ymod)/sd)**2)
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
#### Input red L dwarf name ####                                                          
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
newred.append('spex_prism_U20636_1551+0941_080713')#L3 # redo with L3 instead of L4!!!
newred.append('u11538_050323')#L4 / 1615+4953 #L5                                                                                                                    
newred.append('spex_prism_2206-4217_080714')#L4                                                                                                       
newred.append('U40005_2249+0044_katelyn')#L4                                                                                                          
newred.append('U20171_0355+1133')#L5                                                                                                                  
newred.append('U10372_JHK_2011dec08')#L5 / 0512-2949      
newred.append('spex_prism_U10074_060821')#L1 / 0117-3403
newred.append('U20037_2M0045+1634_adam')#L2
newred.append('spex_prism_2002-0521_080908')#5
newred.append('J215434.68-105530.7')#L5
newred.append('J032642.32-210207.2')#L5
newred.append('J155258.94+294847.9')#L0
newred.append('J035727.02-441730.6')#L0
newred.append('J171113.48+232632.8')#L0
stype = []

# field red Ls
data = asciitable.read('Files/fieldred/fieldred_refined.txt')
desi = data['designation']
spty = data['spt']
color =data['J-Ks']
name = data['name']
#data = asciitable.read('Files/redL6.txt')
#desi = data['filename']
#spty = data['spt']
#color = data['J-Ks']

J = [14.776,14.832,15.066,15.799,15.389,15.376,15.861,15.262,16.436,15.768,15.669,15.522,15.797,17.1,14.982,15.934,16.319,16.789,15.555,16.587,14.05,15.463,15.178,13.059,15.316,16.44,16.134,13.48,14.367,14.499]
K = [13.269,13.097,13.5,14.035,13.702,13.756,14.065,13.615,14.438,13.854,13.659,13.71,14.148,15.3,12.963,14.004,14.31,14.306,13.609,14.358,11.526,13.285,13.49,11.37,13.417,14.2,13.92,12.02,12.907,13.056]
JK = []
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

st_final=[]
name_final=[]
J_final=[]
K_final=[]
JK_final=[]
#time = '20150612-163558'
#time = '20160113-220717'
time = '20160115-134851'
for l in (14,19):
    red = newred[l]
    
# For newreds                                                           

    if l == 0 or l == 1 or l == 2 or l == 3 or l == 4 or l == 5 or l == 6 or l == 27 or l == 28 or l == 29:
        stype.append(0)                                                
        field = field_name[0]                                          
        JKavg = JKav[0]
        shape.append(sh[0])
    if l == 7 or l == 22:
        stype.append(1)                                                
        field = field_name[1]
        JKavg = JKav[1]
        shape.append(sh[1])
    if l == 8 or l == 9 or l == 23:
        sty=2
        stype.append(2)                                                
        field = field_name[2]
        JKavg = JKav[2]
        shape.append(sh[2])
    if l == 10 or l == 11 or l == 12 or l == 16:                                             
        stype.append(3)                                                
        field = field_name[3]
        JKavg = JKav[3]
        shape.append(sh[3])
    if l == 13 or l == 14 or l == 15 or l == 18 or l == 19: 
        sty=4
        stype.append(4)      
        field = field_name[4]
        JKavg = JKav[4]
        shape.append(sh[4])
    if l == 20 or l == 21 or l == 24 or l == 25 or l == 26 or l == 17:   
        stype.append(5)                                                                                                                         
        field = field_name[5]
        JKavg = JKav[5]
        shape.append(sh[5])

#    stype.append(6)
#    field = field_name[6]
    

    if red == 'U20001':
        continue

#    if J[l]-K[l]-JKavg < 0.1:
#        continue
    
# Red L dwarf                                                                                                                                                                 
    if red == 'spex_prism_0141-4633_040905' or red == 'spex_prism_U10074_060821' or red == 'U20037_2M0045+1634_adam' or red == 'spex_prism_1726+1538_080713' or red == 'spex_prism_U20636_1551+0941_080713' or red == 'spex_prism_2002-0521_080908':
        datar = asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(red))
        war = datar['wave']
        flur = datar['flux']
        fer = datar['ferr']
    elif red == 'J215434.68-105530.7' or red == 'J032642.32-210207.2' or red ==  'J155258.94+294847.9' or red == 'J035727.02-441730.6' or red == 'J171113.48+232632.8':
        datar = asciitable.read('Files/fieldred/%s.txt'%(red))
        war = datar['wavelength']
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

    
    print  red,sty
    name_final.append(red)
#    st_final.append(stype[l])
#    J_final.append(J[l])
#    K_final.append(K[l])
#    JK_final.append(J[l]-K[l]-JKavg)
#    if newred[l] == 'u11538_050323' or newred[l] == 'U20171_0355+1133':
#        continue

    if red == 'spex_prism_U20636_1551+0941_080713':
        samples = cPickle.load(open('pickles/spex_prism_U20636_1551+0941_080713_L3_like_20151214-143921.pkl','rb'))
    elif red == 'u11538_050323':
#        samples = cPickle.load(open('pickles/u11538_050323_L5_like_20151215-202111.pkl','rb'))
        samples = cPickle.load(open('pickles/u11538_050323_L6_like_20151215-161107.pkl','rb'))
    elif red == 'J215434.68-105530.7':
        samples = cPickle.load(open('pickles/J215434.68-105530.7_L5_like_20151215-224505.pkl','rb'))
    elif red == 'J032642.32-210207.2' or red ==  'J155258.94+294847.9' or red == 'J035727.02-441730.6' or red == 'J171113.48+232632.8':
        samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(red,stype[l],time), 'rb'))
    else:

#        samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,stype[l],time), 'rb'))
        samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,sty,time),'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    a_emc.append(a_mcmc[0])
    b_emc.append(b_mcmc[0])
#    print a_mcmc[0], b_mcmc[0]
# Calculate J-K color
#    dJK.append(J[l] - K[l]-JKavg)#delta J-K
#    JK.append(J[l] - K[l])#J-K

# field object spectrum
#    dataf = asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(field), data_start=95, data_end=540)
    if sty==4:
        dataf = asciitable.read('Files/%s.txt' %(field),data_start=70, data_end=520)
    else:    
        dataf = asciitable.read('Files/%s.txt' %(field), data_start=95, data_end=540)
    waf = dataf['col1']                                                                                                           
    fluf = dataf['col2']                                                                                                          
    fef = dataf['col3']     
#    waf = dataf['wave']
#    fluf = dataf['flux']
#    fef = dataf['ferr']
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
    chi2_lg = chisq(fluf_norm,flur_norm,fef)
#    chi2_lg = 
    if np.isnan(chi2_lg):
        chi2_lg = 0
#    chi2_lg[np.where(np.isnan(chi2_lg))]=0

    locale.format("%d", chi2_lg, grouping=True)
    chi_lg.append(chi2_lg)

#,len(waf)-1))

# Deredden low-g object & calculate chi squared
#    for j in range(0, ab):
#        Qchi2 = np.interp(waf, w, Qchi[j])
    Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
    Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
    Qint = np.interp(waf, w, Qint)    
    flur_der = flur * np.exp(Qint)
    flur_der_norm = flur_der/nanmean(flur_der)
    chi2_lg_der = chisq(fluf_norm,flur_der_norm,fef)
    dofh = len(waf)-5
#,len(waf)-5)
    if np.isnan(chi2_lg_der):
        chi2_lg_der = 0

    chi_lg_der.append(chi2_lg_der)

#,len(waf)-5))
#print len(JK)
    plt.plot(waf, fluf_norm, 'k', label = 'standard')
    plt.plot(waf, flur_norm, 'r', label='low-g, $\chi^2$=%3s, DOF=%s'%(format(float('%.3g'%(chi2_lg)),'n'),len(waf)))
    plt.plot(waf, flur_der_norm, 'g', label='dereddened low-g, $\chi^2$=%3s, DOF=%s'%(format(float('%.3g'%(chi2_lg_der)),'n'),dofh))
    plt.xlim(0.8, 2.5)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
    leg = plt.legend(fancybox=True, loc=1, prop={'size':10})
    leg.get_frame().set_alpha(0.5)
#    plt.title('%s, L%s'%(red, stype[l]))
    plt.savefig('Plots/dered_%s_L%s_%s.pdf'%(red, sty, timestr))
#    plt.show()
    plt.clf()

    dchi.append(chi2_lg - chi2_lg_der)
    chirat.append(chi2_lg/chi2_lg_der)
quit()
#Lowgdata={'filename':name_final,'sptype':st_final,'J':J_final,'Ks':K_final,'J-Ks':JK_final}
#asciitable.write(lowgdata,'Files/lowgred/lowgred.txt')
#quit()
dchi_f = []
chirat_f = []
shap = []
stypeo = []
chi_fg = []
dcol = []
chi_fg_der = []
# red field obj
#timer = '20150505-145343'
timer = '20160114-200357'
for ll in range (0, len(desi)):

    if desi[ll] == 'None' or desi[ll] == 'J2224438-015852':
#or desi[ll] == 'J121233.92+020626.7' or desi[ll] == 'J005110.83-154417.1' or desi[ll] == 'J110009.49+495745.4' or desi[ll] == 'J213952.22+214839.4' or desi[ll] == 'J033703.73-175806.6' or desi[ll] == 'J020823.83+273738.8' or desi[ll] == 'J062445.89-452150.9':
        continue

    if ll > 0 and desi[ll] == desi[ll-1]:
        continue
    des = desi[ll]
    spt = str(spty[ll])
    st = int(spt[1])
    col = color[ll]
    nam = name[ll]
    if st != 4:
        continue
    
    if des == 'None':
        des = nam
    
    if col-JKav[st] < 0.1:
        continue

    print st, des
#    continue
#    quit()

    samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(des,st,timer), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0\
)))
    ao_emc.append(a_mcmc[0])
    bo_emc.append(b_mcmc[0])
    co_emc.append(c_mcmc[0])
    no_emc.append(n_mcmc[0])

# for delta J-K
    dcol.append(col-JKav[st])
#    colo.append(col)

# field object spectrum
    dataf = asciitable.read('Files/%s.txt' %(field),data_start=35)
    waf = dataf['col1']                                                                                                           
    fluf = dataf['col2']                                                                                                          
    fef = dataf['col3']     
#    dataf=asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(field), data_start=95, data_end=540)
#    waf = dataf['wave']
#    fluf = dataf['flux']
#    fef = dataf['ferr']
    fef[np.where(fef==1.0)]=np.nan
# Normalize if comparing to red object                                                                                                                                    
    fluf_norm = fluf/nanmean(fluf)
    fef = fef/nanmean(fluf)
#    fef[np.where(np.isnan(fef))]=1.0

# Red L dwarf                                                                                                                                                             
    datar = asciitable.read('Files/fieldred/%s.txt' %(des))
    war = datar['wavelength']
    flur = datar['flux']
    fer = datar['ferr']

# Interpolate red spectra onto field spectra                                                                                                                              
    flur = np.interp(waf, war, flur)
    fer = np.interp(waf, war, fer)
# Normalize if comparing to field object                                                                                                                                  
    flur_norm = flur/nanmean(flur)
    fer = fer/nanmean(flur)

# Calculate chi squared 
    chi2_fg = chisq(fluf_norm,flur_norm,fef)
#,len(waf)-1)
    if np.isnan(chi2_fg):
        chi2_fg = 0
#    chi2_lg[np.where(np.isnan(chi2_lg))]=0
    locale.format("%d", chi2_fg, grouping=True)    
    chi_fg.append(chi2_fg)

#,len(waf)-1))

# Deredden low-g object & calculate chi squared
#    for j in range(0, ab):
#        Qchi2 = np.interp(waf, w, Qchi[j])
    Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
    Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
    Qint = np.interp(waf, w, Qint)    
    flur_der = flur * np.exp(Qint)
    flur_der = flur_der/nanmean(flur_der)
    chi2_fg_der = chisq(fluf_norm,flur_der,fef)
#,len(waf)-5)
    if np.isnan(chi2_fg_der):
        chi2_fg_der = 0

    chi_fg_der.append(chi2_fg_der)

#,len(waf)-5))
#print len(JK)
    plt.plot(waf, fluf_norm, 'k', label = 'standard')
    plt.plot(waf, flur_norm, 'r', label='field-g, $\chi^2$=%3s, DOF=%s'%(format(float('%.3g'%(chi2_fg)),'n'),len(waf)))
    plt.plot(waf, flur_der, 'g', label='dereddened field-g, $\chi^2$=%3s, DOF=%s'%(format(float('%.3g'%(chi2_fg_der)),'n'),dofh))
    plt.xlim(0.8, 2.5)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
    leg = plt.legend(fancybox=True, loc=1, prop={'size':10})
    leg.get_frame().set_alpha(0.5)
#    plt.title('%s, L%s'%(des, st))
    plt.savefig('Plots/dered_%s_%s_%s.pdf'%(des, st, timestr))
#    plt.show()
    plt.clf()
#    quit()
    dchi_f.append(chi2_fg - chi2_fg_der)
    chirat_f.append(chi2_fg/chi2_fg_der)

quit()
plota = plt.subplot(111)
#plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
plota.plot(dJK, dchi, 'go')
#plota.plot(dJK, chi_lg_der, 'g.')
#plota.set_xlim(-0.5, 3)
#plota.set_ylim(0, 0.5)
plota.axes.grid(True, linestyle = '-', color = '0.75')
#plota.set_title('$\Delta\chi^2$ vs $\Delta$(J - K)', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
plota.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
plota.set_ylabel('$\Delta\chi^2$',fontsize=14)
plota.plot(dcol, dchi_f, 'mo')
plt.savefig('Plots/dchisq_djk_%s.pdf'%timestr)
plt.clf()
#plt.savefig('Plots/dchisq_lg_djk_%s.pdf'%timestr)
#plt.clf()
plotb = plt.subplot(111)
#plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
plotb.plot(dJK, chirat, 'go')
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
plotb.plot(dcol, chirat_f, 'mo')
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

