############
# plot MCMC parameter vs delta(J-K) color using emcee results
# Use different symbol for parameters exceeding limits
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

# input lowg Ls
lowg = asciitable.read('Files/lowgred/lowgred.txt')
lgname = lowg['filename']
lgst = lowg['sptype']
lgJK = lowg['J-Ks']
lgJ = lowg['J']
lgK = lowg['Ks']
# field red Ls
#data = cPickle.load(open('Files/fieldred/fieldred.pkl'))
data = asciitable.read('Files/fieldred/fieldred_refined.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']
color =data['J-Ks']
name = data['name']
J = [14.776,14.832,15.066,15.799,15.389,15.376,15.861,15.262,16.436,15.768,15.669,15.522,15.797,17.1,14.982,15.934,16.319,16.789,15.555,16.587,14.05,15.463,15.178,13.059,15.316]
K = [13.269,13.097,13.5,14.035,13.702,13.756,14.065,13.615,14.438,13.854,13.659,13.71,14.148,15.3,12.963,14.004,14.31,14.306,13.609,14.358,11.526,13.285,13.49,11.37,13.417]
JK = []
dJK = []
JKav = [1.33, 1.33, 1.67, 1.63, 1.86, 1.46, 1.89]
sh = ['o', 's', 'v', 'd', 'p', '8']
marky = ['mo', 'm^', 'md', 'mp', 'mh', 'm8']
markf = ['go', 'g^', 'gd', 'gp', 'gh', 'g8']
shape = []

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
#lgname = []
lgname0 = []
lgname_conf = []
a_emc0 = []
n_emc0 = []
dJK0 = []
#time = '20150612-163558'  
#time = '20160113-220717'
for l in range(0, len(lgname)):
    red = lgname[l]
    stype.append(lgst[l])
    J.append(lgJ[l])
    K.append(lgK[l])
    JKavg = JKav[stype[l]]
#    if J[l]-K[l]-JKavg < 0.1:
#        continue
    
    print stype[l], red
    if red == 'spex_prism_U20636_1551+0941_080713':
        samples = cPickle.load(open('pickles/spex_prism_U20636_1551+0941_080713_L3_like_20151214-143921.pkl','rb'))
    elif red == 'u11538_050323':
#        samples = cPickle.load(open('pickles/u11538_050323_L5_like_20151215-202111.pkl','rb'))
        samples = cPickle.load(open('pickles/u11538_050323_L6_like_20151215-161107.pkl','rb'))
    elif red == 'U20171_0355+1133':
        samples = cPickle.load(open('pickles/U20171_0355+1133_L6_like_20151215-161107.pkl','rb'))
    elif red == 'J215434.68-105530.7':
        samples = cPickle.load(open('pickles/J215434.68-105530.7_L5_like_20151215-224505.pkl','rb'))
    elif red == 'J032642.32-210207.2' or red ==  'J155258.94+294847.9' or red == 'J035727.02-441730.6' or red == 'J171113.48+232632.8':
        samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(red,stype[l],time), 'rb'))
    elif stype[l] == 4:
        time = '20160115-134851'
        samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,stype[l],time), 'rb'))
    elif stype[l] == 2:
        time = '20160111-233120'
        samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,stype[l],time), 'rb'))
    else:
        time = '20150612-163558'
        samples = cPickle.load(open('pickles/%s_L%s_like_%s.pkl'%(red,stype[l],time), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    a_emc0.append(a_mcmc[0])
#    b_emc.append(b_mcmc[0])
#    a_plus.append(a_mcmc[1])
#    a_minus.append(a_mcmc[2])
#    b_plus.append(b_mcmc[1])
#    b_minus.append(b_mcmc[2])
#    c_emc.append(c_mcmc[0])
    n_emc0.append(n_mcmc[0])
#    n_plus.append(n_mcmc[1])
#    n_minus.append(n_mcmc[2])
#    c_plus.append(c_mcmc[1])
#    c_minus.append(c_mcmc[2])
#    print a_mcmc[0], a_mcmc[1], a_mcmc[2]
# Calculate J-K color
#    dJK0.append(J[l] - K[l]-JKavg)#delta J-K
    dJK0.append(lgJK[l])
#    JK.append(J[l] - K[l])#J-K
#    print (n_mcmc[0]),(n_mcmc[1]),(n_mcmc[2])

#    continue

# both a and b hitting limits
    if red == 'U40005_2249+0044_katelyn':
        symb = 'gv'
        dJKab.append(J[l]-K[l]-JKavg)
        a_ab.append(a_mcmc[0])
        a_abp.append(a_mcmc[1])
        a_abm.append(a_mcmc[2])
        n_ab.append(n_mcmc[0])
        n_abp.append(n_mcmc[1])
        n_abm.append(n_mcmc[2])
#        plota.errorbar(dJK[l], a_mcmc[0], yerr = (a_mcmc[1], a_mcmc[2]), fmt = 'gv', alpha = 0.8)
#    print dJK[l], a_emc[l], a_plus[l], a_minus[l]
#        plt.show()
#        quit()
# b hitting limit hard
    elif red == 'U10372_JHK_2011dec08' or red == 'U10381_JHK_2011dec08' or red == 'U10397' or red == 'U20171_0355+1133' or red == 'U20198_2M0501-0010_dagny' or red == 'spex_prism_0323-4631_U20157' or red == 'u11538_050323' or red == 'spex_prism_U20636_1551+0941_080713' or red == 'U20001' or red == 'U20048' or red == 'U20098' or red == 'spex_prism_2315+0617_U20993' or red == 'spex_prism_U10141_060821' or red == 'spex_prism_0141-4633_040905' or red == 'spex_prism_0032-4405_U20026' or red == 'J215434.68-105530.7' or red == 'J032642.32-210207.2' or red == 'J155258.94+294847.9' or red == 'J035727.02-441730.6' or red == 'J171113.48+232632.8' or red == 'j2213-2136prismtc1' or red == 'u11538_050323':
        symb = 'g>'
        dJKb.append(J[l]-K[l]-JKavg)
        a_b.append(a_mcmc[0])
        a_bp.append(a_mcmc[1])
        a_bm.append(a_mcmc[2])
        n_b.append(n_mcmc[0])
        n_bp.append(n_mcmc[1])
        n_bm.append(n_mcmc[2])
#        plota.errorbar(dJK[l], a_emc[l], yerr = (a_plus[l], a_minus[l]), fmt = 'gD', alpha = 0.8)

# both a & b have double peaks
    elif red == 'spex_prism_1726+1538_080713':
        symb = 'gp'
        dJKw.append(J[l]-K[l]-JKavg)
        a_abw.append(a_mcmc[0])
        a_abwp.append(a_mcmc[1])
        a_abwm.append(a_mcmc[2])
        n_abw.append(n_mcmc[0])
        n_abwp.append(n_mcmc[1])
        n_abwm.append(n_mcmc[2])
#        plota.errorbar(dJK[l], a_emc[l], yerr = (a_plus[l], a_minus[l]), fmt = 'gh', alpha = 0.8)
    
    else:
        symb = 'go'
        dJK.append(J[l]-K[l]-JKavg)
        a_emc.append(a_mcmc[0])
        a_plus.append(a_mcmc[1])
        a_minus.append(a_mcmc[2])
        n_emc.append(n_mcmc[0])
        n_plus.append(n_mcmc[1])
        n_minus.append(n_mcmc[2])
        lgname_conf.append(red)
        print red
    lgname0.append(red)
    shape.append(symb)

# Write delta(J-K),a,N in a file for future use
data = {'filename': lgname0 , 'd(J-K)': dJK0, 'a': a_emc0, 'N': n_emc0}
asciitable.write(data, 'Files/emceeresults_lg.txt')

#        plota.errorbar(dJK[l], a_emc[l], yerr = (a_plus[l], a_minus[l]), fmt = 'go', alpha = 0.8)
# Write confident delta(J-K),a,N in a file for future use
dataa = {'filename': lgname_conf , 'd(J-K)': dJK, 'a': a_emc, 'N': n_emc}
asciitable.write(dataa, 'Files/emceeresults_lg_conf.txt')

#quit()
   
#quit()
print len(lgname0)
#plt.show()
#quit()
stypeo = []
colo = []
dcol = []
shap = []
dcolo0 = []
dcolo1 = []
dcolo2 = []
dcolo3 = []
dcolo4 = []
dcolo5 = []
dcolo6 = []
dcolo7 = []
ao_emc0 = []
ao_emc1 = []
ao_emc2 = []
ao_emc3 = []
ao_emc4 = []
ao_emc5 = []
ao_emc6 = []
ao_emc7 = []
ao_plus0 = []
ao_plus1 = []
ao_plus2 = []
ao_plus3 = []
ao_plus4 = []
ao_plus5 = []
ao_plus6 = []
ao_plus7 = []
ao_minus0 = []
ao_minus1 = []
ao_minus2 = []
ao_minus3 = []
ao_minus4 = []
ao_minus5 = []
ao_minus6 = []
ao_minus7 = []
no_emc0 = []
no_emc1 = []
no_emc2 = []
no_emc3 = []
no_emc4 = []
no_emc5 = []
no_emc6 = []
no_emc7 = []
no_plus0 = []
no_plus1 = []
no_plus2 = []
no_plus3 = []
no_plus4 = []
no_plus5 = []
no_plus6 = []
no_plus7 = []
no_minus0 = []
no_minus1 = []
no_minus2 = []
no_minus3 = []
no_minus4 = []
no_minus5 = []
no_minus6 = []
no_minus7 = []
fgname = []
fgname0 = []
# red field obj
timer = '20150505-145343'
for ll in range (0, len(desi)):

    if desi[ll] == 'None' or desi[ll] == 'J2224438-015852' or desi[ll] == 'J213952.22+214839.4' or desi[ll] == 'J231531.39+061714.2':
        continue


    if ll > 0 and desi[ll] == desi[ll-1]:
        continue
    des = desi[ll]
    spt = str(spty[ll])
    st = int(spt[1])
    col = color[ll]
    nam = name[ll]

    if col-JKav[st] < 0.1:
        continue

    print st, des, col-JKav[st]

    if st == 4:
        time = '20160114-200357'
    elif st == 6:
        time = '20151229-125357'
    else:
        time = timer
    samples = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(des,st,time), 'rb'))
    samples[:, 2] = np.exp(samples[:, 2])
    a_mcmc, b_mcmc, s_mcmc, n_mcmc, c_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    ao_emc.append(a_mcmc[0])
#    bo_emc.append(b_mcmc[0])
#    ao_plus.append(a_mcmc[1])
#    ao_minus.append(a_mcmc[2])
#    bo_plus.append(b_mcmc[1])
 #   bo_minus.append(b_mcmc[2])
 #   co_emc.append(c_mcmc[0])
    no_emc.append(n_mcmc[0])
 #   no_plus.append(n_mcmc[1])
 #   no_minus.append(n_mcmc[2])
 #   co_plus.append(c_mcmc[1])
 #   co_minus.append(c_mcmc[2])

# a hitting limit hard
    if desi[ll] == 'J132629.64-003832.6' or desi[ll] == 'J010752.42+004156.3' or desi[ll] == '2MASSWJ2244316+204343':
        symb = 'm<'
        dcolo1.append(col-JKav[st])
        ao_emc1.append(a_mcmc[0])
        ao_plus1.append(a_mcmc[1])
        ao_minus1.append(a_mcmc[2])
        no_emc1.append(n_mcmc[0])
        no_plus1.append(n_mcmc[1])
        no_minus1.append(n_mcmc[2])
    
# b hitting limit hard
    elif desi[ll] == 'J032642.32-210207.2' or desi[ll] == 'J024111.57-032659.1' or desi[ll] == 'J221207.17+343033.3' or desi[ll] == 'J143832.74+572216.8' or desi[ll] == 'J151240.82+340350.0' or desi[ll] == 'J005110.83-154417.1' or desi[ll] == '153113.38+164128.7' or desi[ll] == 'J021922.07+050630.9' or desi[ll] == 'J171113.48+232632.8' or desi[ll] == 'J065230.58+471036.5' or desi[ll] == 'J121233.92+020626.7' or desi[ll] == 'J215434.68-105530.7' or desi[ll] == 'J122815.34-154736.4' or desi[ll] == 'J033703.73-175806.6' or desi[ll] == 'J035727.02-441730.6' or desi[ll] == 'J233925.43+350716.3' or desi[ll] == 'J155258.94+294847.9' or desi[ll] == '2MASSJ21481628+400359':
        dcolo2.append(col-JKav[st])
        ao_emc2.append(a_mcmc[0])
        ao_plus2.append(a_mcmc[1])
        ao_minus2.append(a_mcmc[2])
        no_emc2.append(n_mcmc[0])
        no_plus2.append(n_mcmc[1])
        no_minus2.append(n_mcmc[2])
        symb = 'm>'

# both a & b hitting limits softly
    elif desi[ll] == 'J083506.10+195303.7' or desi[ll] == 'J083542.14-081920.1':
        symb = 'ms'
        dcolo3.append(col-JKav[st])
        ao_emc3.append(a_mcmc[0])
        ao_plus3.append(a_mcmc[1])
        ao_minus3.append(a_mcmc[2])
        no_emc3.append(n_mcmc[0])
        no_plus3.append(n_mcmc[1])
        no_minus3.append(n_mcmc[2])

# a limits softly
    elif desi[ll] == 'J035822.60-411604.9':
#        symb = 'md'
        symb = 'm<'
        dcolo4.append(col-JKav[st])
        ao_emc4.append(a_mcmc[0])
        ao_plus4.append(a_mcmc[1])
        ao_minus4.append(a_mcmc[2])
        no_emc4.append(n_mcmc[0])
        no_plus4.append(n_mcmc[1])
        no_minus4.append(n_mcmc[2])

# b hitting limits softly
    elif desi[ll] == 'J090546.55+562312.9' or desi[ll] == 'J123927.44+551537.3' or desi[ll] == 'J031013.90-275645.8' or desi[ll] == 'J110009.49+495745.4' or desi[ll] == 'J044743.13-193603.7' or desi[ll] == 'J032740.93-314815.5':
# or desi[ll] == 'J231747.36-483849.4':
#        symb = 'mD'
        symb = 'm>'
        dcolo5.append(col-JKav[st])
        ao_emc5.append(a_mcmc[0])
        ao_plus5.append(a_mcmc[1])
        ao_minus5.append(a_mcmc[2])
        no_emc5.append(n_mcmc[0])
        no_plus5.append(n_mcmc[1])
        no_minus5.append(n_mcmc[2])

# b have double peaks
    elif desi[ll] == 'J051616.03-333202.6' or  desi[ll] == 'J023559.93-233120.5':
        symb = 'mh'
        dcolo6.append(col-JKav[st])
        ao_emc6.append(a_mcmc[0])
        ao_plus6.append(a_mcmc[1])
        ao_minus6.append(a_mcmc[2])
        no_emc6.append(n_mcmc[0])
        no_plus6.append(n_mcmc[1])
        no_minus6.append(n_mcmc[2])

# a & b has twin peaks
    elif desi[ll] == 'J060222.12+633636.5':
        symb = 'mp'
        dcolo7.append(col-JKav[st])
        ao_emc7.append(a_mcmc[0])
        ao_plus7.append(a_mcmc[1])
        ao_minus7.append(a_mcmc[2])
        no_emc7.append(n_mcmc[0])
        no_plus7.append(n_mcmc[1])
        no_minus7.append(n_mcmc[2])

    else:
        dcolo0.append(col-JKav[st])
        ao_emc0.append(a_mcmc[0])
        ao_plus0.append(a_mcmc[1])
        ao_minus0.append(a_mcmc[2])
        no_emc0.append(n_mcmc[0])
        no_plus0.append(n_mcmc[1])
        no_minus0.append(n_mcmc[2])
        fgname0.append(desi[ll])
        
    
    if des == 'None':
        des = nam
#    print st, des
    fgname.append(desi[ll])
# for delta J-K
    dcol.append(col-JKav[st])
    colo.append(col)
#    print col-JKav[st], des, st
#    fgname.append(desi[ll])
# Write confident delta(J-K),a,N in a file for future use
dataoo = {'filename': fgname0 , 'd(J-K)': dcolo0, 'a': ao_emc0, 'N': no_emc0}
asciitable.write(dataoo, 'Files/emceeresults_fg_conf.txt')

# Write delta(J-K),a,N in a file for future use
datao = {'filename': fgname , 'd(J-K)': dcol, 'a': ao_emc, 'N': no_emc}
asciitable.write(datao, 'Files/emceeresults_fg.txt')

#quit()

plota = plt.subplot(111)
plota.axes.grid(True, linestyle = '-', color = '0.75')
plota.errorbar(dJKab, a_ab, yerr = (a_abp, a_abm), fmt = 'gs', alpha = 0.7, markersize = 5)
plota.errorbar(dJKb, a_b, yerr = (a_bp, a_bm), fmt = 'gD', alpha = 0.7, markersize = 5)
#plota.errorbar(dJKbs, a_bs, yerr = (a_bsp, a_bsm), fmt = 'gd', alpha = 0.7, markersize = 6)
#plota.errorbar(dJKbs, a_bs, yerr = (a_bsp, a_bsm), fmt = 'gD', alpha = 0.8, markersize = 6)
#plota.errorbar(dJKw, a_abw, yerr = (a_abwp, a_abwm), fmt = 'gh', alpha = 0.8, markersize = 10)
plota.errorbar(dJK, a_emc, yerr = (a_plus, a_minus), fmt = 'go', alpha = 0.7, markersize = 12)
#plota.errorbar(dcolo6, ao_emc6, yerr = (ao_plus6, ao_minus6), fmt = 'mh', alpha=0.5, markersize = 10)
#plota.errorbar(dcolo7, ao_emc7, yerr = (ao_plus7, ao_minus7), fmt = 'mp', alpha=0.5, markersize = 9)
plota.errorbar(dcolo5, ao_emc5, yerr = (ao_plus5, ao_minus5), fmt = 'mD', alpha=0.5, markersize = 5)
plota.errorbar(dcolo4, ao_emc4, yerr = (ao_plus4, ao_minus4), fmt = 'md', alpha=0.5, markersize = 6)
plota.errorbar(dcolo3, ao_emc3, yerr = (ao_plus3, ao_minus3), fmt = 'ms', alpha=0.5, markersize = 5)
#plota.errorbar(dcolo5, ao_emc5, yerr = (ao_plus5, ao_minus5), fmt = 'mD', alpha=0.5, markersize = 6)
#plota.errorbar(dcolo4, ao_emc4, yerr = (ao_plus4, ao_minus4), fmt = 'md', alpha=0.5, markersize = 7)
#plota.errorbar(dcolo3, ao_emc3, yerr = (ao_plus3, ao_minus3), fmt = 'ms', alpha=0.5, markersize = 6)
plota.errorbar(dcolo2, ao_emc2, yerr = (ao_plus2, ao_minus2), fmt = 'mD', alpha=0.5, markersize = 5)
plota.errorbar(dcolo1, ao_emc1, yerr = (ao_plus1, ao_minus1), fmt = 'md', alpha=0.5, markersize = 6)
plota.errorbar(dcolo0, ao_emc0, yerr = (ao_plus0, ao_minus0), fmt = 'mo', alpha=0.5, markersize = 12)
plota.set_ylim(0.1, 0.45)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
#plota.set_title('Mean Effective Radius vs J - K', fontsize=15)
plota.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plota.set_xlabel('J - K',fontsize=14)
plota.set_ylabel('Mean Effective Radius [${\mu}m$]',fontsize=14)
plt.savefig('Plots/a_djk_%s.pdf'%timestr)

plt.clf()

plotn = plt.subplot(111)
plotn.axes.grid(True, linestyle = '-', color = '0.75')
plotn.errorbar(dJKab, n_ab, yerr = (n_abp, n_abm), fmt = 'gs', alpha = 0.7, markersize = 5)
plotn.errorbar(dJKb, n_b, yerr = (n_bp, n_bm), fmt = 'gD', alpha = 0.7, markersize = 5)
#plotn.errorbar(dJKbs, n_bs, yerr = (n_bsp, n_bsm), fmt = 'gd', alpha = 0.7, markersize = 6)
plotn.errorbar(dJK, n_emc, yerr = (n_plus, n_minus), fmt = 'go', alpha = 0.7, markersize = 12)
#plotn.errorbar(dJKbs, n_bs, yerr = (n_bsp, n_bsm), fmt = 'gD', alpha = 0.8, markersize = 6)
#plotn.errorbar(dJKw, n_abw, yerr = (n_abwp, n_abwm), fmt = 'gh', alpha = 0.8, markersize = 10)
plotn.errorbar(dcolo0, no_emc0, yerr = (no_plus0, no_minus0), fmt = 'mo', alpha=0.5, markersize = 12)
#plotn.errorbar(dcolo6, no_emc6, yerr = (no_plus6, no_minus6), fmt = 'mh', alpha=0.5, markersize = 10)
#plotn.errorbar(dcolo7, no_emc7, yerr = (no_plus7, no_minus7), fmt = 'mp', alpha=0.5, markersize = 9)
plotn.errorbar(dcolo5, no_emc5, yerr = (no_plus5, no_minus5), fmt = 'mD', alpha=0.5, markersize = 5)
plotn.errorbar(dcolo4, no_emc4, yerr = (no_plus4, no_minus4), fmt = 'md', alpha=0.5, markersize = 6)
plotn.errorbar(dcolo3, no_emc3, yerr = (no_plus3, no_minus3), fmt = 'ms', alpha=0.5, markersize = 5)
#plotn.errorbar(dcolo5, no_emc5, yerr = (no_plus5, no_minus5), fmt = 'mD', alpha=0.5, markersize = 6)
#plotn.errorbar(dcolo4, no_emc4, yerr = (no_plus4, no_minus4), fmt = 'md', alpha=0.5, markersize = 7)
#plotn.errorbar(dcolo3, no_emc3, yerr = (no_plus3, no_minus3), fmt = 'ms', alpha=0.5, markersize = 6)
plotn.errorbar(dcolo2, no_emc2, yerr = (no_plus2, no_minus2), fmt = 'mD', alpha=0.5, markersize = 5)
plotn.errorbar(dcolo1, no_emc1, yerr = (no_plus1, no_minus1), fmt = 'md', alpha=0.5, markersize = 6)
plotn.set_ylim(0, 8)
#plotn.set_title('Column Density vs J - K', fontsize=15)
#plotn.set_title('Column Density vs J - K', fontsize=15)
plotn.set_xlabel('$\Delta$(J - K)',fontsize=14)
#plotn.set_xlabel('J - K',fontsize=14)
plotn.set_ylabel('Column Density [$10^8 cm^{-2}$]',fontsize=14)

plt.savefig('Plots/n_djk_%s.pdf'%timestr)
plt.clf()


print len(colo)
quit()
#plota = plt.subplot(111)
#plota.errorbar(colo, ao_emc, yerr = (ao_plus, ao_minus), fmt = 'ms', alpha=0.5)
v#plota.errorbar(JK, a_emc, yerr = (a_plus, a_minus), fmt = 'gs', alpha = 0.8)
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
