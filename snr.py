############
# Calculate and plot S/N for observed reddening
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

data = asciitable.read('Files/fieldred/fieldred.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']
color = data['J-Ks']
name = data['name']
Q_grid = []
Qchi = []

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
#newred.append('spex_prism_0141-4633_040905') #L0
#newred.append('spex_prism_0323-4631_U20157') #L0                                                                                 
#newred.append('spex_prism_U10141_060821') #L0          
#newred.append('spex_prism_0032-4405_U20026')#L0                                                                                                       
#newred.append('U20098')#L0                                                                                                                            
#newred.append('j2213-2136prismtc1')#L0                                                                                                                
#newred.append('spex_prism_2315+0617_U20993')#L0                                                                                                       
#newred.append('spex_prism_U10074_060821') #L1
#newred.append('U10381_JHK_2011dec08')#L1                                                                                                              
#newred.append('U10397')#L2                                                                                                                            
#newred.append('U20048')#L2                                                                                                                            
#newred.append('U20037_2M0045+1634_adam') #L2
#newred.append('spex_prism_1726+1538_080713')#3
#newred.append('U20001')#L3                                                                                                                            
#newred.append('spex_prism_2208+2921_U40004')#L3                                                                                                       
#newred.append('spex_prism_U20636_1551+0941_080713')#L4
#newred.append('2massj0126+1428spex')#L4                                                                                                               
#newred.append('U20198_2M0501-0010_dagny')#L4                                                                                                          
#newred.append('U20622')#L4                                                                                                                            
#newred.append('u11538_050323')#L4                                                                                                                     
#newred.append('spex_prism_2206-4217_080714')#L4                                                                      
#newred.append('U40005_2249+0044_katelyn')#L4                                                                                                          
#newred.append('spex_prism_2002-0521_080908')#L5
#newred.append('U20171_0355+1133')#L5                                                                                                                  
#newred.append('U10372_JHK_2011dec08')#L5       
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

J = [14.776,14.832,15.066,15.799,15.389,15.376,15.861,15.262,16.436,15.768,15.669,15.522,15.797,17.1,14.982,15.934,16.319,16.789,15.555
,16.587,14.05,15.463,15.178,13.059,15.316]
K = [13.269,13.097,13.5,14.035,13.702,13.756,14.065,13.615,14.438,13.854,13.659,13.71,14.148,15.3,12.963,14.004,14.31,14.306,13.609,14.358,11.526,13.285,13.49,11.37,13.417]
JK = []
JKav = [1.33, 1.33, 1.67, 1.63, 1.86, 1.46]

#### Red Field L dwarfs ####                                                                                                      
#oldred = []
#oldred.append('spex_prism_U20791_060821')#L4   
#oldred.append('2MASS1239+5515')#L5
#oldred.append('2MASSJ02355993-2331205')#L1
#oldred.append('2MASSJ00250365+4759191')#L4
#oldred.append('2MASSJ06244595-4521548')#L5
#oldred.append('2MASSWJ1841086+311727')#L4
#oldred.append('2MASSJ08354256-0819237')#5

#stype = []
#pdff = PdfPages('fitmodels.pdf')
#with PdfPages('triangle.pdf') as pdft:
#pdft = PdfPages('triangle.pdf')
for k in range(0, len(newred)):
    red = newred[k]
# For newreds                                                            
#    if k == 0 or k == 1 or k == 2 or k == 3 or k == 4 or k == 5 or k == 6:                          
#        stype.append(0)                                                
#        field = field_name[0]                                          
#    if k == 7 or k == 8:                                                         
#        stype.append(1)                                                
#        field = field_name[1]                                          
#    if k ==9 or k == 10 or k == 11: 
#        stype.append(2)                                                
#        field = field_name[2]                                          
#    if k == 12 or k == 13 or k == 14:               
 #       sp = 3
#        stype.append(3)                                                
#        field = field_name[3]                                          
#    if k == 15 or k == 16 or k == 17 or k == 18 or k == 19 or k == 20 or k == 21: 
#        sp = 4
#        stype.append(4)      
#        field = field_name[4]
#    if k == 22 or k == 23 or k == 24:   
#        sp = 5
#        stype.append(5)                                                                                                                         
#        field = field_name[5]                                                                                                                                                          

    if k == 0 or k == 1 or k == 2 or k == 3 or k == 4 or k == 5 or k == 6:
        stype=0
        field = field_name[0]
        JKavg = JKav[0]
    if k == 7 or k == 22:
        stype=1
        field = field_name[1]
        JKavg = JKav[1]
    if k == 8 or k == 9 or k == 23:
        stype=2
        field = field_name[2]
        JKavg = JKav[2]
    if k == 10 or k == 11 or k == 12:
        stype=3
        field = field_name[3]
        JKavg = JKav[3]
    if k == 13 or k == 14 or k == 15 or k == 16 or k == 17 or k == 18 or k == 19:
        stype=4
        field = field_name[4]
        JKavg = JKav[4]
    if k == 20 or k == 21 or k == 24:
        stype=5                                                                                                                
        field = field_name[5]
        JKavg = JKav[5]
    
#    if k != 10 or k != 17 or k != 20:
#        continue
#    if red != 'spex_prism_1726+1538_080713' or red != 'u11538_050323' or red != 'U20171_0355+1133':
 #       continue

    if J[k]-K[k]-JKavg < 0.1:
        continue

    print stype, red
#    data = asciitable.read('../../Documents/ASTRO/Mie/Miesults/normalized/extinction_%s.unnorm'%(red))
    data = asciitable.read('inputs/extinction_%s.unnorm'%(red))
    wave = data['wave']
    frat = data['normalized flux ratio']
    error_frat = data['flux ratio error']
    frat_norm = frat - np.mean(frat)

    figy = plt.figure() 
    ax = figy.add_subplot(211)        
#    ax.plot(wave,frat_norm,'k')
    ax.plot(wave,frat,'k')
    ax.plot(wave, error_frat, '0.5')
#    ax.set_xlabel('$\lambda~[\mu m]$')
    ax.set_ylabel('Observed Reddening (unnormalized)')
#    ax.set_ylabel('Observed Reddening (normalized)')
    plt.title('%s, L%s'%(red,stype))

#    figy.savefig('Plots/rederror_s_%s_%s_%s.pdf'%(red,stype, timestr)) 
#    pdff.savefig(chainfig) 
#    plt.close()
    
    snr = np.abs(frat / error_frat)
#    figyy = plt.figure()
    ax1 = figy.add_subplot(212)
    ax1.plot(wave, snr, 'b')
#    ax1.set_xlabel('$\lambda~[\mu m]$')
    ax1.set_ylabel('S/N')
    figy.savefig('Plots/snr_s_%s_%s_%s.pdf'%(red,stype, timestr))
#    plt.show()
    plt.close()

stype = []
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
    frat_norm = frat - np.mean(frat)

    figr = plt.figure() 
    ax2 = figr.add_subplot(211)
        
    ax2.plot(wave,frat,'k')
#    ax2.plot(wave,frat_norm,'k')
    ax2.plot(wave, error_frat, '0.5')
#    ax2.set_xlabel('$\lambda~[\mu m]$')
#    ax2.set_ylabel('Observed Reddening (normalized)')
    ax2.set_ylabel('Observed Reddening (unnormalized)')
    plt.title('%s, L%s'%(des,st))
#    figr.savefig('Plots/rederror_s_%s_%s_%s.pdf'%(des,st, timestr)) 
#    pdff.savefig(chainfig) 
#    plt.close()
    
    snr = np.abs(frat / error_frat)
#    figrr = plt.figure()
    ax3 = figr.add_subplot(212)
    ax3.plot(wave, snr, 'b')
#    ax3.set_xlabel('$\lambda~[\mu m]$')
    ax3.set_ylabel('S/N')
    figr.savefig('Plots/snr_s_%s_%s_%s.pdf'%(des,st, timestr))
    plt.close()
#    pdff.savefig(chainfig2) 

#    pdft.savefig()
#plt.clf()
#pdff.close()
#pdft.close()
