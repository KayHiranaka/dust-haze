# Run 2-dimensional KS test on field and low-g L dwarfs
import asciitable
import numpy as np
import ks2d

# two distributions for field and low-g
lg = asciitable.read('../Files/emceeresults_lg.txt')
fg = asciitable.read('../Files/emceeresults_fg.txt')

x1 = lg['d(J-K)']       # Delta(J-K) for lg
y1 =  lg['a']      # a , N for lg
y11 = lg['N']
x2 =  fg['d(J-K)']      # Delta(J-K) for fg
y2 =  fg['a']
y22 = fg['N']      # a , N for fg 

# call ks2d function to compute probability
prob = ks2d.ks2d(x1,y1,x2,y2)
print 'probability drawn from same parent distribution (a)=', prob
#print '(should be close to 1, or at least not near zero)'

print ''
prob = ks2d.ks2d(x1,y11,x2,y22)
print 'probability drawn from same parent distribution (N)=', prob
quit()
# two distributions for *confident* field and low-g objects
lg = asciitable.read('../Files/emceeresults_lg_conf.txt')
fg = asciitable.read('../Files/emceeresults_fg_conf.txt')

x1 = lg['d(J-K)']       # Delta(J-K) for lg
y1 =  lg['a']      # a , N for lg
y11 = lg['N']
x2 =  fg['d(J-K)']      # Delta(J-K) for fg
y2 =  fg['a']
y22 = fg['N']      # a , N for fg 

# call ks2d function to compute probability
prob = ks2d.ks2d(x1,y1,x2,y2)
print 'probability drawn from same parent distribution (conf a)=', prob
#print '(should be close to 1, or at least not near zero)'

print ''
prob = ks2d.ks2d(x1,y11,x2,y22)
print 'probability drawn from same parent distribution (conf N)=', prob

quit()
#;make a subtle shift in one of the distributions
proba=ks2d.ks2d(x1+0.1,y1,x2,y2)
print 'probability drawn from same parent distribution=', proba
print '(should be a little less than last value)'

print ''

#;now see how they compare when distributions have different spreads
x2=np.random.standard_normal((2000,))*5
y2=np.random.standard_normal((2000,))*5

probab = ks2d.ks2d(x1,y1,x2,y2)
print 'probability drawn from same parent distribution=', probab
print '(should be much less than 1)'


#end

