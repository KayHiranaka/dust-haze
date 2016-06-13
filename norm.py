# Divide field&red L dwarfs and normalize observed extinction and Mie extinction curves
##############

from pylab import *
import asciitable
import numpy
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import nanmean
import time
timestr = time.strftime("%Y%m%d-%H%M%S")

#### Input data ####                                                                                                               \
#### Input field L dwarf name ####                                                                                                 \
field_name = []
field_name.append('2M0345_L0')
field_name.append('2M2130_L1_kc')
field_name.append('Kelu-1')
field_name.append('u11291_1506+1321_050323')
#field_name.append('2MASSJ21580457-1550098')
field_name.append('U12101_2158-1550_davy.fits')
field_name.append('2M1507')
field_name.append('2MASSIJ1010148-040649')

#### Input red L dwarf name ####                                                                                                   \
red_name = []
red_name.append('spex_prism_0141-4633_040905')
red_name.append('spex_prism_U10074_060821')
red_name.append('U20037_2M0045+1634_adam')
red_name.append('spex_prism_1726+1538_080713')
red_name.append('spex_prism_U20636_1551+0941_080713')
red_name.append('spex_prism_2002-0521_080908')

#### New L dwarf gammas ####                                                                                            
newred = []
newred.append('spex_prism_0323-4631_U20157') #L0                                                         
newred.append('spex_prism_U10141_060821') #L0                                                                 
newred.append('U10397')#L2                                                                                             
newred.append('U20048')#L2                                                                                             
newred.append('spex_prism_0032-4405_U20026')#L0
newred.append('U20098')#L0
newred.append('j2213-2136prismtc1')#L0
newred.append('spex_prism_2315+0617_U20993')#L0
newred.append('U10381_JHK_2011dec08')#L1
newred.append('U10397')#L2
newred.append('U20001')#L3
newred.append('spex_prism_2208+2921_U40004')#L3
newred.append('2massj0126+1428spex')#L4
newred.append('U20198_2M0501-0010_dagny')#L4
newred.append('U20622')#L4
newred.append('u11538_050323')#L4#0355
newred.append('spex_prism_2206-4217_080714')#L4
newred.append('U40005_2249+0044_katelyn')#L4
newred.append('U20171_0355+1133')#L5#0355
newred.append('U10372_JHK_2011dec08')#L5

#### Red Field L dwarfs ####
oldred = []
oldred.append('J215434.68-105530.7')#L5
oldred.append('spex_prism_U20791_060821')#L4

data = asciitable.read('Files/fieldred/fieldred_refined.txt')
#desi = data['filename']
#spty = data['spt']
#color = data['J-Ks']
desi = data['designation']
spty = data['spt']
color = data['J-Ks']
name = data['name']

reff = np.mgrid[0.05:0.4:8j]
veff = np.mgrid[0.1:1.0:19j]

rgrid = np.mgrid[0.05:0.4:101j] # reff and veff on a fine grid                                                                                                                                     
vgrid = np.mgrid[0.1:0.6:101j]
stype = []
for l in range(0, len(desi)):
#    field = field_name[6]
#    red = red_name[l]
    #red = newred[l]
#    red = oldred[l]
#    stype.append(5)
#    field = field_name[5]
# For newreds
#    if l == 0 or l == 1:
#        stype.append(0)
#        field = field_name[0]
#    if l == 2 or l == 3:
#        stype.append(2)
#        field = field_name[2]
#    if l == 4 or l == 5 or l == 6 or l == 7:
#        stype.append(0)
#        field = field_name[0]
#    if l == 8:
#        stype.append(1)
#        field = field_name[1]
#    elif l == 9:
#        stype.append(2)
#        field = field_name[2]
#    if l == 10 or l == 11:
#        stype.append(3)
#        field = field_name[3]
#    if l == 12 or l == 13 or l == 14 or l == 15 or l == 16 or l == 17:
#        stype.append(4)
#        field = field_name[4]
#    if l == 18 or l == 19:
#        stype.append(5)
#        field = field_name[5]
#    else:
#        continue
#field red
    spt = str(spty[l])
    spt = int(spt[1])
    des = desi[l]
    col = color[l]
    red = des
    if spt != 4:
        continue
    field = field_name[spt]
# Read field L dwarf spectra 0.9-2.5um
#    dataf = asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(field), data_start=95, data_end=540)
    dataf = asciitable.read('Files/%s.txt' %(field), data_start=70, data_end=520)
    
#    waf = dataf['wave']
#    fluf = dataf['flux']
#    fef = dataf['ferr']
    waf = dataf['col1']
    fluf = dataf['col2']
    fef = dataf['col3']
    fef = fef[np.where(np.isfinite(fef))]
    fluf = fluf[np.where(np.isfinite(fef))]
    waf = waf[np.where(np.isfinite(fef))]
# Normalize if comparing to red object      
    fluf_norm = fluf/nanmean(fluf)
    fef_norm = fef/nanmean(fluf)
    print red
    
# Red L dwarf
#    datar = asciitable.read('../../ASTRO/Mie/ApplytoData/%s.dat' %(red))
#newred
#    datar = asciitable.read('../../ASTRO/Mie/Archive/newdata/%s.fits.txt' %(red))
#fieldred
    datar = asciitable.read('Files/fieldred/%s.txt' %(red), data_start=95, data_end=540)
# red L6
#    if red == '2MASSIJ0103320+193536':
#        datar = asciitable.read('Files/lowgred/%s.txt' %(red))
#    else:
#        datar = asciitable.read('Files/fieldred/%s.txt' %(red))

# Newreds
#    war = datar['col1']
#    flur = datar['col2']
#    fer = datar['col3']
    war = datar['wavelength']
    flur = datar['flux']
    fer = datar['ferr']
    fer= fer[np.where(np.isfinite(fer))]
    flur = flur[np.where(np.isfinite(fer))]
    war= war[np.where(np.isfinite(fer))]

# Interpolate red spectra onto field spectra
    flur = numpy.interp(waf, war, flur)
    fer = numpy.interp(waf, war, fer)
# Normalize if comparing to field object
#    flur_norm = flur/nanmean(flur)
#    fer_norm = fer/nanmean(flur)
# Flux ratio of field/red                                                                                                                                                                                                                     
    frat = numpy.log(fluf/flur)

    wav = waf[numpy.where(numpy.isfinite(frat))]
    frat = frat[numpy.where(numpy.isfinite(frat))]
    fef = fef[numpy.where(numpy.isfinite(frat))]
    fluf = fluf[numpy.where(numpy.isfinite(frat))]
    fer = fer[numpy.where(numpy.isfinite(frat))]
    flur = flur[numpy.where(numpy.isfinite(frat))]

# Error propagation 
    fef[numpy.where(numpy.isnan(fef))] = 0
    fer[numpy.where(numpy.isnan(fer))] = 0
    fef[numpy.where(fef == 1.0)] = 0
    fer[numpy.where(fer == 1.0)] = 0

    error_frat = numpy.sqrt((fef/fluf)**2 + (fer/flur)**2)
#    print numpy.mean(error_frat)

#    fratio = frat - numpy.mean(frat) # Normalize observed extinction by subtracting mean value (subtraction instead of division because it's ln(flux ratio))
    data = {'wave':wav, 'normalized flux ratio': frat, 'flux ratio error':error_frat}
    asciitable.write(data, 'inputs/extinction_%s_unc.unnorm'%(red))
    plot(wav, frat)
#    plot(wav, error_frat/nanmean(error_frat),'0.5')
    savefig('Plots/ext_%s_unc.pdf'%(red))
    clf()
#    combofig = plt.figure()
#    ax = combofig.add_subplot(1,1,1)
#    ax.plot(waf, flur_norm, 'r')
#    ax.plot(waf, fluf_norm, 'k')
#    ax.set_xlim(0.9, 2.5)
#    ax.set_ylabel('Normalized Flux',fontsize = 14)
#    ax.set_xlabel('Wavelength [$\mu$m]',fontsize = 14)
#    savefig('Plots/spec_%s_%s.pdf'%(red,timestr))
#    show()
#    quit()
#    ext = asciitable.read('inputs/extinction_%s.unnorm'%(red))
#    wavee = ext['wave']
#    frate = ext['normalized flux ratio']
#    axe = combofig.add_subplot(2,1,2)
#    axe.plot(wavee, frate, 'g')
#    axe.set_xlim(0.9, 2.5)
 #   axe.set_xlabel('wavelength [$\mu$m]',fontsize = 14)
#    axe.set_ylabel('Observed Reddening')
#    savefig('../../../../../Dropbox/Codes/Paper/spec.eps')
#    savefig('Plots/spec_ext_%s_%s.pdf'%(red,timestr))
#    show()
quit()
a = []
b = []
#    print type(fratio)
# Mie extinction coefficients
for i in range(0, len(rgrid)):
    for j in range(0, len(vgrid)):
        a.append(rgrid[i])
        b.append(vgrid[j])
            
# Mie Coefficients averaged over Hansen distribution and interpolated on grids                                                                                                                    
        data = asciitable.read('../../../Mie/Miesults/Hansen/bigrid/Mg2SiO4_a%sum_b%s.grid2'%(rgrid[i], vgrid[j]))
        wave = data['wave']
        Qex = data['Qext(interpolated on a grid)']

#    for h in range(0, len(reff)):
#        colors = ['b', 'c', 'g', 'y', 'm', 'r', 'k']
#        for k in range(0, len(veff)):

#            Qex = numpy.interp(wav, wave, Qex)
        Qex_norm = Qex - numpy.mean(Qex) # also normalizing Mie coefficients by subtractiong mean
#            print type(Qex_norm)
        data = {'wave':wave, 'normalized Qext':Qex_norm}
        asciitable.write(data, '../../../Mie/Miesults/Hansen/Norm/Mg2SiO4_a%sum_b%s.norm'%(rgrid[i], vgrid[j]))
            
#            xlim(0.8, 2.6)
#            xlabel('Wavelength (${\mu}m$)')
#            ylabel('$Q_{ext}$')
#            title('r=%sum v=%sum'%(rgrid[i], vgrid[j]))
#            plot(wav, Qex_norm, 'g', label='normalized')
#            plot(wav, Qex, 'g--', label='unnormalized')
#            leg = legend(fancybox=True, loc='center', prop={'size':12})
#            leg.get_frame().set_alpha(0.5)
#            savefig('../../../Mie/Miesults/Plots/Mie_normalized')
#            show()
#            print i, j
#quit()

#    xlim(0.8, 2.6)
#    xlabel('Wavelength (${\mu}m$)')
#    ylabel('$ln({I_{0}}/{I})$')
#    title('L%s'%(l))
#    plot(wav, fratio, 'k', label='normalized')
#    plot(wav, frat, 'k--', label='unnormalized')
#    savefig('../../../Mie/Miesults/Plots/L%s_normalized.pdf'%(l))
#    leg = legend(fancybox=True, loc='center', prop={'size':12})
#    leg.get_frame().set_alpha(0.5)
#    show()
