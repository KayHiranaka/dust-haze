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
import time
timestr = time.strftime("%Y%m%d-%H%M%S")

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
field_name.append('Kelu-1_L2')
field_name.append('u11291_1506+1321_050323')
field_name.append('2M2158-15_L4')
field_name.append('2M1507')

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
newred.append('spex_prism_U20636_1551+0941_080713')#L4
newred.append('u11538_050323')#L4 / 1615+4953                                                                                                                    
newred.append('spex_prism_2206-4217_080714')#L4                                                                                                       
newred.append('U40005_2249+0044_katelyn')#L4                                                                                                          
newred.append('U20171_0355+1133')#L5                                                                                                                  
newred.append('U10372_JHK_2011dec08')#L5 / 0512-2949      
newred.append('spex_prism_U10074_060821')#L1 / 0117-3403
newred.append('U20037_2M0045+1634_adam')#L2
newred.append('spex_prism_2002-0521_080908')#5
stype = []

# field red Ls
data = asciitable.read('Files/fieldred/fieldred.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']
color =data['J-Ks']
name = data['name']

J = [14.776,14.832,15.066,15.799,15.389,15.376,15.861,15.262,16.436,15.768,15.669,15.522,15.797,17.1,14.982,15.934,16.319,16.789,15.555,16.587,14.05,15.463,15.178,13.059,15.316]
K = [13.269,13.097,13.5,14.035,13.702,13.756,14.065,13.615,14.438,13.854,13.659,13.71,14.148,15.3,12.963,14.004,14.31,14.306,13.609,14.358,11.526,13.285,13.49,11.37,13.417]
JK = []
dJK = []
JKav = [1.33, 1.33, 1.67, 1.63, 1.86, 1.46]
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
chi_lg_ab = []
chi_lg_der_ab = []
chi_lg_b = []
chi_lg_der_b = []
chi_lg_bs = []
chi_lg_der_bs = []
chi_lg_w = []
chi_lg_der_w = []
dchi = []
chirat = []
dchi_ab = []
chirat_ab = []
dchi_b = []
chirat_b = []
dchi_bs = []
chirat_bs = []
dchi_w = []
chirat_w = []
a_ab = []
a_abp = []
a_abm = []
a_a = []
a_ap = []
a_am = []
a_b = []
a_bp = []
a_bm = []
a_bs = []
a_bsp = []
a_bsm = []
a_abw = []
a_abwp = []
a_abwm = []
dJKab = []
dJKb = []
dJKbs = []
dJKw = []
n_ab = []
n_abp = []
n_abm = []
n_a = []
n_ap = []
n_am = []
n_b = []
n_bp = []
n_bm = []
n_bs = []
n_bsp = []
n_bsm = []
n_abw = []
n_abwp = []
n_abwm = []

time = '20150612-163558'
for l in range (0, len(newred)):
    red = newred[l]

# For newreds                                                           

    if l == 0 or l == 1 or l == 2 or l == 3 or l == 4 or l == 5 or l == 6:
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
        stype.append(2)                                                
        field = field_name[2]
        JKavg = JKav[2]
        shape.append(sh[2])
    if l == 10 or l == 11 or l == 12:                                             
        stype.append(3)                                                
        field = field_name[3]
        JKavg = JKav[3]
        shape.append(sh[3])
    if l == 13 or l == 14 or l == 15 or l == 16 or l == 17 or l == 18 or l == 19: 
        stype.append(4)      
        field = field_name[4]
        JKavg = JKav[4]
        shape.append(sh[4])
    if l == 20 or l == 21 or l == 24:   
        stype.append(5)                                                                                                                         
        field = field_name[5]
        JKavg = JKav[5]
        shape.append(sh[5])

    if J[l]-K[l]-JKavg < 0.1:
        continue
# Red L dwarf                                                                                                                                                                 
    if red == 'spex_prism_0141-4633_040905' or red == 'spex_prism_U10074_060821' or red == 'U20037_2M0045+1634_adam' or red == 'spex_prism_1726+1538_080713' or red == 'spex_prism_U20636_1551+0941_080713' or red == 'spex_prism_2002-0521_080908':
        datar = asciitable.read('../Mie/ApplytoData/%s.dat' %(red))
        war = datar['wave']
        flur = datar['flux']
        fer = datar['ferr']
    else:
#        continue
#newred                                                                                                                                                                   
        datar = asciitable.read('../Mie/Archive/newdata/%s.fits.txt' %(red))                                                                                           

# Newreds                                                                                                                                                                 
        war = datar['col1']                                                                                                                                                  
        flur = datar['col2']                                                                                                                                                 
        fer = datar['col3']             

    print stype[l], red
#    if newred[l] == 'u11538_050323' or newred[l] == 'U20171_0355+1133':
#        continue
# field object spectrum
    dataf = asciitable.read('../Mie/ApplytoData/%s.dat' %(field), data_start=95, data_end=540)

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
#    fer = fer/nanmean(flur)
#    print fef, flur_norm
# Calculate chi squared 

# emcee chain
    samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,stype[l],time), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    a_emc.append(a_mcmc[0])
    b_emc.append(b_mcmc[0])
#    print a_mcmc[0], b_mcmc[0]
# Calculate J-K color
#    dJK.append(J[l] - K[l]-JKavg)#delta J-K
#    JK.append(J[l] - K[l])#J-K
# both a and b hitting limits                                                                                                                            
    if newred[l] == 'U40005_2249+0044_katelyn':
        symb = 'gv'
        dJKab.append(J[l]-K[l]-JKavg)
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_lg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_lg_der):
            chi2_lg_der = 0
        chi_lg_der_ab.append(chi2_lg_der)
        chi2_lg_ab = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_lg_ab):
            chi2_lg_ab = 0    
        chi_lg_ab.append(chi2_lg_ab)
        dchi_ab.append(chi2_lg - chi2_lg_der)
        chirat_ab.append(chi2_lg/chi2_lg_der)
                                                                                                          
# b hitting limit hard                                                                                                                                   
    elif newred[l] == 'U10372_JHK_2011dec08' or newred[l] == 'U10381_JHK_2011dec08' or newred[l] == 'U10397' or newred[l] == 'U20171_0355+1133' or newred[l] == 'U20198_2M0501-0010_dagny' or newred[l] == 'spex_prism_0323-4631_U20157' or newred[l] == 'u11538_050323' or newred[l] == 'spex_prism_U20636_1551+0941_080713' or newred[l] == 'U20001' or newred[l] == 'U20048' or newred[l] == 'U20098' or newred[l] == 'spex_prism_2315+0617_U20993' or newred[l] == 'spex_prism_U10141_060821' or newred[l] == 'spex_prism_0141-4633_040905' or newred[l] == 'spex_prism_0032-4405_U20026':
        symb = 'g>'
        dJKb.append(J[l]-K[l]-JKavg)
        chi2_lg_b = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_lg_b):
            chi2_lg_b = 0    
        chi_lg_b.append(chi2_lg_b)
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_lg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_lg_der):
            chi2_lg_der = 0
        chi_lg_der_b.append(chi2_lg_der)
        dchi_b.append(chi2_lg_b - chi2_lg_der)
        chirat_b.append(chi2_lg_b/chi2_lg_der)

# both a & b have double peaks                                                                                                                           
    elif newred[l] == 'spex_prism_1726+1538_080713':
        symb = 'gh'
        dJKw.append(J[l]-K[l]-JKavg)
        chi2_lg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_lg):
            chi2_lg = 0    
        chi_lg_w.append(chi2_lg)
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_lg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_lg_der):
            chi2_lg_der = 0
        chi_lg_der_w.append(chi2_lg_der)
        dchi_w.append(chi2_lg - chi2_lg_der)
        chirat_w.append(chi2_lg/chi2_lg_der)

    else:
        symb = 'go'
        dJK.append(J[l]-K[l]-JKavg)
        chi2_lg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_lg):
            chi2_lg = 0    
        chi_lg.append(chi2_lg)
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_lg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_lg_der):
            chi2_lg_der = 0
        chi_lg_der.append(chi2_lg_der)
        dchi.append(chi2_lg - chi2_lg_der)
        chirat.append(chi2_lg/chi2_lg_der)


#,len(waf)-1))

# Deredden low-g object & calculate chi squared
#    for j in range(0, ab):
#        Qchi2 = np.interp(waf, w, Qchi[j])

#,len(waf)-5))
#print len(JK)
#    plt.plot(waf, fluf_norm, 'k', label = 'standard')
#    plt.plot(waf, flur_norm, 'r', label='low-g, $\chi^2$=%s'%(chi2_lg))
#    plt.plot(waf, flur_der, 'g', label='dereddened low-g, $\chi^2$=%s'%(chi2_lg_der))
#    plt.xlim(0.8, 2.5)
#    plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
#    leg = plt.legend(fancybox=True, loc=1, prop={'size':10})
#    leg.get_frame().set_alpha(0.5)
#    plt.title('%s, L%s'%(red, stype[l]))
#    plt.savefig('Plots/dered_%s_%s_%s.pdf'%(red, stype[l], timestr))
#    plt.show()
#    plt.clf()
#    quit()

#quit()
dcolo0 = []
dcolo1 = []
dcolo2 = []
dcolo3 = []
dcolo4 = []
dcolo5 = []
dcolo6 = []
dcolo7 = []
dchi_f0 = []
chirat_f0 = []
dchi_f1 = []
chirat_f1 = []
dchi_f2 = []
chirat_f2 = []
dchi_f3 = []
chirat_f3 = []
dchi_f4 = []
chirat_f4 = []
dchi_f5 = []
chirat_f5 = []
dchi_f6 = []
chirat_f6 = []
dchi_f7 = []
chirat_f7 = []
shap = []
stypeo = []
chi_fg0 = []
chi_fg_der0 = []
chi_fg1 = []
chi_fg_der1 = []
chi_fg2 = []
chi_fg_der2 = []
chi_fg3 = []
chi_fg_der3 = []
chi_fg4 = []
chi_fg_der4 = []
chi_fg5 = []
chi_fg_der5 = []
chi_fg6 = []
chi_fg_der6 = []
chi_fg7 = []
chi_fg_der7 = []
dcol = []
# red field obj
#timer = '20150505-145343'
for ll in range (0, len(desi)):

    if desi[ll] == 'None' or desi[ll] == 'J2224438-015852' or desi[ll] == 'J213952.22+214839.4' or desi[ll] == 'J231531.39+061714.2':
#or desi[ll] == 'J121233.92+020626.7' or desi[ll] == 'J005110.83-154417.1' or desi[ll] == 'J110009.49+495745.4' or desi[ll] == 'J213952.22+214839.4' or desi[ll] == 'J033703.73-175806.6' or desi[ll] == 'J020823.83+273738.8' or desi[ll] == 'J062445.89-452150.9':
        continue

    if ll > 0 and desi[ll] == desi[ll-1]:
        continue
    des = desi[ll]
    spt = spty[ll]
    col = color[ll]
    nam = name[ll]
    if des == 'None':
        des = nam

    if spt == 10 or spt == 10.5:
        stypeo.append(0)
        st = 0
        field = field_name[0]
        shap.append(sh[0])
    if spt == 11 or spt == 11.5:
        st = 1
        stypeo.append(1)
        field = field_name[1]
        shap.append(sh[1])
    if spt == 12 or spt == 12.5:
        st = 2
        stypeo.append(2)
        field = field_name[2]
        shap.append(sh[2])
    if spt == 13 or spt == 13.5:
        st = 3
        stypeo.append(3)
        field = field_name[3]
        shap.append(sh[3])
    if spt == 14 or spt == 14.5:
        st = 4
        stypeo.append(4)
        field = field_name[4]
        shap.append(sh[4])
    if spt == 15 or spt == 15.5:
        st = 5
        stypeo.append(5)
        field = field_name[5]
        shap.append(sh[5])
    
    if col-JKav[st] < 0.1:
        continue

    print st, des
#    continue
#    quit()
# field object spectrum
    dataf = asciitable.read('../Mie/ApplytoData/%s.dat' %(field), data_start=95, data_end=540)
    waf = dataf['wave']
    fluf = dataf['flux']
    fef = dataf['ferr']
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

# emcee chain
    samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(des,st,time), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0\
)))
    ao_emc.append(a_mcmc[0])
    bo_emc.append(b_mcmc[0])
    co_emc.append(c_mcmc[0])
    no_emc.append(n_mcmc[0])

# a hitting limit hard                                                                                                                                      
    if desi[ll] == 'J132629.64-003832.6' or desi[ll] == 'J010752.42+004156.3':
        symb = 'm<'
        dcolo1.append(col-JKav[st])
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg1.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der1.append(chi2_fg_der)
        dchi_f1.append(chi2_fg - chi2_fg_der)
        chirat_f1.append(chi2_fg/chi2_fg_der)

# b hitting limit hard                                                                                                                                      
    elif desi[ll] == 'J032642.32-210207.2' or desi[ll] == 'J024111.57-032659.1' or desi[ll] == 'J221207.17+343033.3' or desi[ll] == 'J143832.74+572216.8' or desi[ll] == 'J151240.82+340350.0' or desi[ll] == 'J005110.83-154417.1' or desi[ll] == '153113.38+164128.7' or desi[ll] == 'J021922.07+050630.9' or desi[ll] == 'J171113.48+232632.8' or desi[ll] == 'J065230.58+471036.5' or desi[ll] == 'J121233.92+020626.7' or desi[ll] == 'J215434.68-105530.7' or desi[ll] == 'J122815.34-154736.4' or desi[ll] == 'J033703.73-175806.6' or desi[ll] == 'J035727.02-441730.6' or desi[ll] == 'J233925.43+350716.3' or desi[ll] == 'J155258.94+294847.9':
        dcolo2.append(col-JKav[st])
        symb = 'm>'
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg2.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der2.append(chi2_fg_der)
        dchi_f2.append(chi2_fg - chi2_fg_der)
        chirat_f2.append(chi2_fg/chi2_fg_der)

# both a & b hitting limits softly                                                                                                                          
    elif desi[ll] == 'J083506.10+195303.7' or desi[ll] == 'J083542.14-081920.1':
        symb = 'ms'
        dcolo3.append(col-JKav[st])
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg3.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der3.append(chi2_fg_der)
        dchi_f3.append(chi2_fg - chi2_fg_der)
        chirat_f3.append(chi2_fg/chi2_fg_der)

# a limits softly                                                                                                                                           
    elif desi[ll] == 'J035822.60-411604.9':
        symb = 'md'
        dcolo4.append(col-JKav[st])
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg4.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der4.append(chi2_fg_der)
        dchi_f4.append(chi2_fg - chi2_fg_der)
        chirat_f4.append(chi2_fg/chi2_fg_der)

# b hitting limits softly                                                                                                                                   
    elif desi[ll] == 'J090546.55+562312.9' or desi[ll] == 'J123927.44+551537.3' or desi[ll] == 'J031013.90-275645.8' or desi[ll] == 'J110009.49+495745.4' or desi[ll] == 'J044743.13-193603.7':
        symb = 'mD'
        dcolo5.append(col-JKav[st])
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg5.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der5.append(chi2_fg_der)
        dchi_f5.append(chi2_fg - chi2_fg_der)
        chirat_f5.append(chi2_fg/chi2_fg_der)

# b have double peaks                                                                                                                                       
    elif desi[ll] == 'J051616.03-333202.6' or  desi[ll] == 'J023559.93-233120.5':
        symb = 'mh'
        dcolo6.append(col-JKav[st])
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg6.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der6.append(chi2_fg_der)
        dchi_f6.append(chi2_fg - chi2_fg_der)
        chirat_f6.append(chi2_fg/chi2_fg_der)

# a & b has twin peaks                                                                                                                                      
    elif desi[ll] == 'J060222.12+633636.5':
        symb = 'mp'
        dcolo7.append(col-JKav[st])
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg7.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der7.append(chi2_fg_der)
        dchi_f7.append(chi2_fg - chi2_fg_der)
        chirat_f7.append(chi2_fg/chi2_fg_der)

    else:
        dcolo0.append(col-JKav[st])
# Calculate chi squared between standard and red L 
        chi2_fg = chisq(fluf_norm,flur_norm,fef)
        if np.isnan(chi2_fg):
            chi2_fg = 0
        chi_fg0.append(chi2_fg)
# Deredden low-g object & calculate chi squared
        Qint=griddata((a, b), Q_grid, (a_mcmc[0],b_mcmc[0]), method='linear')
        Qint = n_mcmc[0]*np.pi*Qint*a_mcmc[0]**2 + c_mcmc[0]
        Qint = np.interp(waf, w, Qint)    
        flur_der = flur * np.exp(Qint)
        flur_der = flur_der/nanmean(flur_der)
        chi2_fg_der = chisq(fluf_norm,flur_der,fef)
        if np.isnan(chi2_fg_der):
            chi2_fg_der = 0
        chi_fg_der0.append(chi2_fg_der)
        dchi_f0.append(chi2_fg - chi2_fg_der)
        chirat_f0.append(chi2_fg/chi2_fg_der)

# for delta J-K
#    dcol.append(col-JKav[st])
#    colo.append(col)


#,len(waf)-5))
#print len(JK)
#    plt.plot(waf, fluf_norm, 'k', label = 'standard')
#    plt.plot(waf, flur_norm, 'r', label='field-g, $\chi^2$=%s'%(chi2_fg))
#    plt.plot(waf, flur_der, 'g', label='dereddened field-g, $\chi^2$=%s'%(chi2_fg_der))
#    plt.xlim(0.8, 2.5)
#    plt.legend(bbox_to_anchor=(1.05, 1), loc=1, borderaxespad=0.)
#    leg = plt.legend(fancybox=True, loc=1, prop={'size':10})
#    leg.get_frame().set_alpha(0.5)
#    plt.title('%s, L%s'%(des, st))
#    plt.savefig('Plots/dered_%s_%s_%s.pdf'%(des, st, timestr))
#    plt.show()
#    plt.clf()
#    quit()
plota = plt.subplot(111)
#plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)#
plota.plot(dJK, dchi, 'go', alpha = 0.8, markersize = 15)
plota.plot(dJKab, dchi_ab, 'gv', alpha = 0.8, markersize = 5)
plota.plot(dJKb, dchi_b, 'g>', alpha = 0.8, markersize = 5)
plota.plot(dJKw, dchi_w, 'gh', alpha = 0.8, markersize = 10)
plota.plot(dcolo0, dchi_f0, 'mo', alpha = 0.8, markersize = 15)
plota.plot(dcolo1, dchi_f1,  'm<', alpha=0.5, markersize = 5)
plota.plot(dcolo2, dchi_f2, 'm>', alpha=0.5, markersize = 5)
plota.plot(dcolo3, dchi_f3, 'ms', alpha=0.5, markersize = 6)
plota.plot(dcolo4, dchi_f4, 'md', alpha=0.5, markersize = 7)
plota.plot(dcolo5, dchi_f5, 'mD', alpha=0.5, markersize = 6)
plota.plot(dcolo6, dchi_f6, 'mh', alpha=0.5, markersize = 10) 
plota.plot(dcolo7, dchi_f7, 'mp', alpha=0.5, markersize = 9) 
#plota.plot(dJK, chi_lg_der, 'g.')
#plota.set_xlim(-0.5, 3)
#plota.set_ylim(0, 0.5)
plota.axes.grid(True, linestyle = '-', color = '0.75')
plota.set_title('$\Delta\chi^2$ vs $\Delta$(J - K)', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
plota.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
plota.set_ylabel('$\Delta\chi^2$',fontsize=14)
plt.savefig('Plots/dchisq_djk_%s.pdf'%timestr)
plt.clf()
#plt.savefig('Plots/dchisq_lg_djk_%s.pdf'%timestr)
#plt.clf()
plotb = plt.subplot(111)
#plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
plotb.plot(dJK, chirat, 'go', alpha = 0.8, markersize = 15)
plotb.plot(dJKab, chirat_ab, 'gv', alpha = 0.8, markersize = 5)
plotb.plot(dJKb, chirat_b, 'g>', alpha = 0.8, markersize = 5)
plotb.plot(dJKw, chirat_w, 'gh', alpha = 0.8, markersize = 10)
plotb.plot(dcolo0, chirat_f0, 'mo', alpha = 0.8, markersize = 15)
plotb.plot(dcolo1, chirat_f1,  'm<', alpha=0.5, markersize = 5)
plotb.plot(dcolo2, chirat_f2, 'm>', alpha=0.5, markersize = 5)
plotb.plot(dcolo3, chirat_f3, 'ms', alpha=0.5, markersize = 6)
plotb.plot(dcolo4, chirat_f4, 'md', alpha=0.5, markersize = 7)
plotb.plot(dcolo5, chirat_f5, 'mD', alpha=0.5, markersize = 6)
plotb.plot(dcolo6, chirat_f6, 'mh', alpha=0.5, markersize = 10) 
plotb.plot(dcolo7, chirat_f7, 'mp', alpha=0.5, markersize = 9) 
#plotb.plot(dJK, chirat, 'go')
#plota.plot(dJK, chi_lg_der, 'g.')
#plota.set_xlim(-0.5, 3)
#plota.set_ylim(0, 0.5)
plotb.axes.grid(True, linestyle = '-', color = '0.75')
plotb.set_title('$\chi^2$ ratio vs $\Delta$(J - K)', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
plotb.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
plotb.set_ylabel('$\chi^2$ ratio',fontsize=14)
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

