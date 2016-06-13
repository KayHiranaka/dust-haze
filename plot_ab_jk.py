############
# plot a, b vs J-K color using emcee results
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
#from scipy.interpolate import griddata
#import matplotlib.axis as axis
import time
timestr = time.strftime("%Y%m%d-%H%M%S")

a_emc = []
b_emc = []
a_plus = []
b_plus = []
a_minus = []
b_minus = []
n_emc = []
c_emc = []
n_plus = []
c_plus = []
n_minus = []
c_minus = []
ao_emc = []
bo_emc = []
ao_plus = []
bo_plus = []
ao_minus = []
bo_minus = []
no_emc = []
co_emc = []
no_plus = []
co_plus = []
no_minus = []
co_minus = []
a_emc0 = []
a_emc1 = []
a_emc2 = []
a_emc3 = []
a_emc4 = []
a_emc5 = []
n_emc0 = []
n_emc1 = []
n_emc2 = []
n_emc3 = []
n_emc4 = []
n_emc5 = []
a_plus0 = []
a_plus1 = []
a_plus2 = []
a_plus3 = []
a_plus4 = []
a_plus5 = []
n_plus0 = []
n_plus1 = []
n_plus2 = []
n_plus3 = []
n_plus4 = []
n_plus5 = []
a_minus0 = []
a_minus1 = []
a_minus2 = []
a_minus3 = []
a_minus4 = []
a_minus5 = []
n_minus0 = []
n_minus1 = []
n_minus2 = []
n_minus3 = []
n_minus4 = []
n_minus5 = []

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
#data = cPickle.load(open('Files/fieldred/fieldred.pkl'))
data = asciitable.read('Files/fieldred/fieldred.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']
color =data['J-Ks']
name = data['name']
#for k in range(len(desi)):
#    if desi[k] == None:
#        desi[k] = 'None'
#        print desi[k]
#    if name[k] == None:
#        name[k]= 'None'
#        print name[k]
#tx = {'designation':desi, 'name':name, 'spt':spty, 'J-Ks':color}
#asciitable.write(tx, 'Files/fieldred/fieldred.txt')
#quit()
#print newred[13] # negative delta J-K
# Read J & K colors
#file = asciitable.read('../../../Mie/Archive/newdata/KaysData.txt')
#J = data['J']
#K = data['K']
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

time = '20150612-163558'
for l in range (0, len(newred)):
    red = newred[l]
#    if l > 0 and red[l] == red[m] for m in range(l):
#        continue
    
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
    print stype[l], red, J[l]-K[l]-JKavg
    if newred[l] == 'u11538_050323' or newred[l] == 'U20171_0355+1133':
        continue
#    samples = cPickle.load(open('../../Documents/ASTRO/codes/kay-repo/Python/%s_%s_area_cropchain.pkl'%(des,stype[l]), 'rb'))
    samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(red,stype[l],time), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    a_emc.append(a_mcmc[0])
    b_emc.append(b_mcmc[0])
    a_plus.append(a_mcmc[1])
    a_minus.append(a_mcmc[2])
    b_plus.append(b_mcmc[1])
    b_minus.append(b_mcmc[2])
    c_emc.append(c_mcmc[0])
    n_emc.append(n_mcmc[0])
    n_plus.append(n_mcmc[1])
    n_minus.append(n_mcmc[2])
    c_plus.append(c_mcmc[1])
    c_minus.append(c_mcmc[2])
    print a_mcmc[0], a_mcmc[1], a_mcmc[2]
# Calculate J-K color
    dJK.append(J[l] - K[l]-JKavg)#delta J-K
    JK.append(J[l] - K[l])#J-K
print len(JK)

stypeo = []
colo = []
dcol = []
shap = []
colo0 = []
colo1 = []
colo2 = []
colo3 = []
colo4 = []
colo5 = []
ao_emc0 = []
ao_emc1 = []
ao_emc2 = []
ao_emc3 = []
ao_emc4 = []
ao_emc5 = []
ao_plus0 = []
ao_plus1 = []
ao_plus2 = []
ao_plus3 = []
ao_plus4 = []
ao_plus5 = []
ao_minus0 = []
ao_minus1 = []
ao_minus2 = []
ao_minus3 = []
ao_minus4 = []
ao_minus5 = []
no_emc0 = []
no_emc1 = []
no_emc2 = []
no_emc3 = []
no_emc4 = []
no_emc5 = []
no_plus0 = []
no_plus1 = []
no_plus2 = []
no_plus3 = []
no_plus4 = []
no_plus5 = []
no_minus0 = []
no_minus1 = []
no_minus2 = []
no_minus3 = []
no_minus4 = []
no_minus5 = []


# red field obj
timer = '20150505-145343'
for ll in range (0, len(desi)):

    if desi[ll] == 'None' or desi[ll] == 'J2224438-015852':
        continue
    if desi[ll] == 'None' or desi[ll] == 'J2224438-015852' or desi[ll] == 'J121233.92+020626.7' or desi[ll] == 'J005110.83-154417.1' or desi[ll] == 'J110009.49+495745.4' or desi[ll] == 'J213952.22+214839.4' or desi[ll] == 'J033703.73-175806.6' or desi[ll] == 'J020823.83+273738.8' or desi[ll] == 'J062445.89-452150.9':
        continue

    if ll > 0 and desi[ll] == desi[ll-1]:
        continue
    des = desi[ll]
    spt = spty[ll]
    col = color[ll]
    nam = name[ll]

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

    print st, des, col-JKav[st]

    samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(des,st,timer), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    ao_emc.append(a_mcmc[0])
    bo_emc.append(b_mcmc[0])
    ao_plus.append(a_mcmc[1])
    ao_minus.append(a_mcmc[2])
    bo_plus.append(b_mcmc[1])
    bo_minus.append(b_mcmc[2])
    co_emc.append(c_mcmc[0])
    no_emc.append(n_mcmc[0])
    no_plus.append(n_mcmc[1])
    no_minus.append(n_mcmc[2])
    co_plus.append(c_mcmc[1])
    co_minus.append(c_mcmc[2])


    if des == 'None':
        des = nam
    print st, des

# for delta J-K
    dcol.append(col-JKav[st])
    colo.append(col)
#    print col-JKav[st], des, st
print len(dcol)
#quit()
plota = plt.subplot(111)
#plota.errorbar(colo, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
#plota.errorbar(JK, a_emc, yerr = (a_plus, a_minus), fmt = 'gs', alpha = 0.8)
#plota.errorbar(colo, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
#plota.errorbar(JK, a_emc, yerr = (a_plus, a_minus), fmt = 'gs', alpha = 0.8)
plota.errorbar(dcol, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
plota.errorbar(dJK, a_emc, yerr = (a_plus, a_minus), fmt = 'gs', alpha = 0.8)
#plota.set_xlim(-0.5, 3)
plota.set_ylim(0, 0.5)
plota.axes.grid(True, linestyle = '-', color = '0.75')
plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
plota.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
plota.set_ylabel('Mean Effective Radius [${\mu}m$]',fontsize=14)
#plt.savefig('../../../Mie/Miesults/Plots/emc_a_jk_%s.pdf'%timestr)
plt.savefig('Plots/emc_a_djk_%s.pdf'%timestr)
#plt.savefig('Paper/emc_a_jk.eps')
#plt.show()
plt.clf()
plotb = plt.subplot(111)
#plotb.errorbar(colo, bo_emc, yerr = (bo_plus, bo_minus), fmt = 'md', alpha=0.5)
#plotb.errorbar(JK, b_emc, yerr = (b_plus, b_minus), fmt = 'gd', alpha = 0.8)
plotb.errorbar(dcol, bo_emc, yerr = (bo_plus, bo_minus), fmt = 'md', alpha=0.5)
plotb.errorbar(dJK, b_emc, yerr = (b_plus, b_minus), fmt = 'gd', alpha = 0.8)
#plotb.set_xlim(-0.5, 3)
plotb.set_ylim(0.0, 1.1)
plotb.axes.grid(True, linestyle = '-', color = '0.75')
plotb.set_title('b vs $\Delta$J - K',fontsize=15)
plotb.set_xlabel('$\Delta$(J - K)',fontsize=14)
plotb.set_ylabel('Effective Variance b',fontsize=14)
#plt.savefig('../../../Mie/Miesults/Plots/emc_b_jk_%s.pdf'%timestr)
#plt.savefig('Plots/emc_b_djk_%s.pdf'%timestr)
#plt.show()
plt.clf()
plotn = plt.subplot(111)
#plotn.errorbar(colo, no_emc, yerr = (no_plus, no_minus), fmt = 'mp', alpha=0.5)
#plotn.errorbar(JK, n_emc, yerr = (n_plus, n_minus), fmt = 'gp', alpha = 0.8)
plotn.errorbar(dcol, no_emc, yerr = (no_plus, no_minus), fmt = 'mp', alpha=0.5)
plotn.errorbar(dJK, n_emc, yerr = (n_plus, n_minus), fmt = 'gp', alpha = 0.8)
#plota.set_xlim(-0.5, 3)
plotn.set_ylim(0, 10)
plotn.axes.grid(True, linestyle = '-', color = '0.75')
#plotn.set_title('Column Density vs J - K', fontsize=15)
plotn.set_title('Column Density vs J - K', fontsize=15)
plotn.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plotn.set_xlabel('J - K',fontsize=14)
plotn.set_ylabel('Column Density [10^8 $cm^{-2}$]',fontsize=14)
#plt.savefig('../../../Mie/Miesults/Plots/emc_n_jk_%s.pdf'%timestr)
plt.savefig('Plots/emc_n_djk_%s.pdf'%timestr)
#plt.savefig('Paper/emc_n_jk.eps')
#plt.show()
plt.clf()
plotc = plt.subplot(111)
#plotc.errorbar(colo, co_emc, yerr = (co_plus, co_minus), fmt = 'mo', alpha=0.5)
#plotc.errorbar(JK, c_emc, yerr = (c_plus, c_minus), fmt = 'go', alpha = 0.8)
#plotb.set_xlim(-0.5, 3)
#plotb.set_ylim(0.0, 1.1)
plotc.axes.grid(True, linestyle = '-', color = '0.75')
plotc.set_title('offset vs $\Delta$J - K',fontsize=15)
plotc.set_xlabel('$\Delta$J - K',fontsize=14)
plotc.set_ylabel('Vertical Offset',fontsize=14)
#plt.savefig('../../../Mie/Miesults/Plots/emc_c_jk_%s.pdf'%timestr)
#plt.savefig('Plots/emc_c_jk_%s.pdf'%timestr)
#plt.show()
plt.clf()
