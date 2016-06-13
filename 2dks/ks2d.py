# Modified IDL code ks2d.pro
import asciitable
import numpy as np
import quad_fracs
import raw_prob_ks as rp
from scipy.stats import pearsonr

#;;compute the pobability that two 2-d arrays are drawn from the same
#;;        distribution.  

def ks2d(x1,y1,x2,y2):
#;x1,y1 points for model 1
#;x2,y2 points for model 2

    da1=0
#;loop through the first array as center points

    n1=x1.size
    n2=x2.size
    for i in range(0,n1):# do begin
        d1=quad_fracs.qf(x1[i],y1[i], x1,y1)
        d2=quad_fracs.qf(x1[i],y1[i],x2,y2)
#        da1=np.amax([da1,np.amax(abs(d1-d2))])
        da1 = np.amax([abs(d1[0]-d2[0]),abs(d1[1]-d2[1]),abs(d1[2]-d2[2]),abs(d1[3]-d2[3])])
#        print da1
#endfor

    da2=0
#;loop throught the second array
    for i in range(0,n2):# do begin
        d1=quad_fracs.qf(x2[i],y2[i], x1,y1)
        d2=quad_fracs.qf(x2[i],y2[i],x2,y2)
#        da2=np.amax([da2,np.amax(abs(d1-d2))])
        da2 = np.amax([abs(d1[0]-d2[0]),abs(d1[1]-d2[1]),abs(d1[2]-d2[2]),abs(d1[3]-d2[3])])
#endfor

#;average the d's
    d=np.mean(da1,da2)

    n=(n1*n2)/(n1+n2)

#;get linear correlation coef for each sample
    r1=pearsonr(x1,y1)
    r2=pearsonr(x2,y2)

    rr=np.sqrt(1.0-0.5*(r1[0]**2+r2[0]**2))

#;ESTIMATE probability using raw K-S prob

#;this line was ain error, I coppied teh wrong equation from NR.
#;     Caught by Stephane Blondin
#;lambda= ( sqrt(n)*d)/(1.+sqrt(1-rr^2)*(0.25-.75/sqrt(n)))

#;new version of Numerical Recipies says:
    lamb = ( np.sqrt(n)*d)/(1.+rr*(0.25-.75/np.sqrt(n)))
#    print rr, lamb
    

#;stop

    prob = rp.rp(lamb)

    return prob
#end



