##################
# Check if MCMC parameters are properly burned-in
##################

from pylab import *
import numpy as np
import asciitable
import matplotlib.pyplot as plt
import emcee
import cPickle
import time
timestr = time.strftime("%Y%m%d-%H%M%S")

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

##### red L dwarfs (field g)

data = asciitable.read('Files/fieldred/fieldred.txt')
desi = data['designation']
spty = data['spt']
#sid = data['source_id']                                                                                                 
color = data['J-Ks']
name = data['name']

# input time
time = '20150422-200357'
J = [14.776,14.832,15.066,15.799,15.389,15.376,15.861,15.262,16.436,15.768,15.669,15.522,15.797,17.1,14.982,15.934,16.319,16.789,15.555
,16.587,14.05,15.463,15.178,13.059,15.316]
K = [13.269,13.097,13.5,14.035,13.702,13.756,14.065,13.615,14.438,13.854,13.659,13.71,14.148,15.3,12.963,14.004,14.31,14.306,13.609,14.358,11.526,13.285,13.49,11.37,13.417]

JKav = [1.33, 1.33, 1.67, 1.63, 1.86, 1.46]


stype = []
l = 0
for k in range(0, len(newred)):
    l = l + 1
    red = newred[k]
    if k == 0 or k == 1 or k == 2 or k == 3 or k == 4 or k == 5 or k == 6:
        stype.append(0)
    if k == 7 or k == 22:
        stype.append(1)
    if k == 8 or k == 9 or k == 23:
        stype.append(2)
    if k == 10 or k == 11 or k == 12:
        stype.append(3)
    if k == 13 or k == 14 or k == 15 or k == 16 or k == 17 or k == 18 or k == 19:
        stype.append(4)
    if k == 20 or k == 21 or k == 24:
        stype.append(5)                                                                                                  
        
    JKavg = JKav[stype[k]]
    if J[k]-K[k]-JKavg < 0.1:
        continue

    cropchain = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(red,stype[l-1],time), 'rb'))
    print cropchain[:,0][0]
    quit()
    plt.plot(cropchain[:,0])
    plt.title('%s,%s'%(red,stype[l-1]))
    plt.ylabel('a')
    plt.savefig('Plots/burnin_a_%s_%s.pdf'%(stype[l-1],red))
    plt.ylim(0, 0.5)
    plt.close()
    plt.plot(cropchain[:,1])
    plt.title('%s,%s'%(red,stype[l-1]))
    plt.ylabel('b')
    plt.ylim(0, 1)
    plt.savefig('Plots/burnin_b_%s_%s.pdf'%(stype[l-1],red))
    plt.close()
    plt.plot(cropchain[:,2])
    plt.title('%s,%s'%(red,stype[l-1]))
    plt.ylabel('log(s)')
    plt.savefig('Plots/burnin_logs_%s_%s.pdf'%(stype[l-1],red))
    plt.close()
    plt.plot(cropchain[:,3])
    plt.title('%s,%s'%(red,stype[l-1]))
    plt.ylabel('N')
    plt.savefig('Plots/burnin_N_%s_%s.pdf'%(stype[l-1],red))
    plt.close()
    plt.plot(cropchain[:,4])
    plt.title('%s,%s'%(red,stype[l-1]))
    plt.ylabel('C')
    plt.savefig('Plots/burnin_C_%s_%s.pdf'%(stype[l-1],red))

    plt.close()

quit()
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
    if spt == 11 or spt == 11.5:
        st = 1
        stype.append(1)
    if spt == 12 or spt == 12.5:
        st = 2
        stype.append(2)
    if spt == 13 or spt == 13.5:
        st = 3
        stype.append(3)
    if spt == 14 or spt == 14.5:
        st = 4
        stype.append(4)
    if spt == 15 or spt == 15.5:
        st = 5
        stype.append(5)
    if col-JKav[st] < 0.1:
        continue
    cropchain = cPickle.load(open('pickles/%s_%s_area_cropchain_%s.pkl'%(des,st,time), 'rb'))
#    print len(cropchain[:,0])
    plt.plot(cropchain[:,0])
    plt.title('%s, %s'%(des,st))
    plt.ylabel('a')
    plt.ylim(0, 0.5)
    plt.savefig('Plots/burnin_a_%s_%s.pdf'%(st,des))
    plt.close()
    plt.plot(cropchain[:,1])
    plt.title('%s, %s'%(des,st))
    plt.ylabel('b')
    plt.ylim(0, 10)
    plt.savefig('Plots/burnin_b_%s_%s.pdf'%(st,des))
    plt.close()
    plt.plot(cropchain[:,2])
    plt.title('%s, %s'%(des,st))
    plt.ylabel('log(s)')
    plt.savefig('Plots/burnin_logs_%s_%s.pdf'%(st,des))
    plt.close()
    plt.plot(cropchain[:,3])
    plt.title('%s, %s'%(des,st))
    plt.ylabel('N')
    plt.savefig('Plots/burnin_N_%s_%s.pdf'%(st,des))
    plt.close()
    plt.plot(cropchain[:,4])
    plt.title('%s, %s'%(des,st))
    plt.ylabel('C')
    plt.savefig('Plots/burnin_C_%s_%s.pdf'%(st,des))
    plt.close()


