# Modified IDL code quad_fracs.pro
import asciitable
import numpy as np

#function quad_fracs, xc,yc,xx,yy
#;compute the fraction of points in each quadrant given the center
#;        point xc,yc.  There's an ambiguity which quadrant to stick
#; things if xc,yc lands on a data point.  If that happens, quad_fracs returns an array where it computes the fractions putting the xc,yc point in each quadrant. 

#;;the equations are not accurate if N < 20 or the returned probability
#;;    is > 0.2.  if the returned is greater than 0.2, then it is a
#;; sign that you can treat them as drawn from the same distribution.

def qf(xc, yc, xx, yy):

    r1=xx[np.where((xx > xc) & (yy > yc))]
    r2=xx[np.where((xx < xc) & (yy > yc))]
    r3=xx[np.where((xx < xc) & (yy < yc))]
    r4=xx[np.where((xx > xc) & (yy < yc))]

    nt=float(xx.size)

#    if r1.max() == -1:
#        n1=0 
#    else:
#        n1=r1.size
#    if r2.max() == -1:
#        n2=0 
#    else:
 #       n2=r2.size
 #   if r3.max() == -1:
 #       n3=0 
 #   else:
 #       n3=r3.size
 #   if r4.max() == -1:
 #       n4=0 
 #   else:
 #       n4=r4.size

    n1,n2,n3,n4 = len(r1),len(r2),len(r3),len(r4)

#    fracs=[n1, n2, n3,n4 ]/nt
    frac1,frac2,frac3,frac4 = n1/nt, n2/nt, n3/nt, n4/nt
#    print frac1,frac2,frac3,frac4
#    print np.sum([frac1,frac2,frac3,frac4])

    if np.sum([frac1,frac2,frac3,frac4]) != 1.: #this means xc,yc hit a point
#        fracs=[n1, n2, n3,n4]
#        diag=[ [1.,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
#        fracs=(fracs+diag)/nt
        n1,n2,n3,n4 = n1+1, n2+1, n3+1, n4+1
        frac1,frac2,frac3,frac4 = n1/nt, n2/nt, n3/nt, n4/nt
#        print 'Ambiguity due to (xc,yc) landing on a datapoint'
#        return
    return [frac1,frac2,frac3,frac4]
# end

