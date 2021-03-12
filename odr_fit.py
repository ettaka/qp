from scipy.odr import *


def linear(p, x):
   k, b = p
   return k*x + b

def fit(x,y):
    linear_model = Model(linear)
    data = RealData(x, y)
    odr = ODR(data, linear_model, beta0=[0., 1.])
    out = odr.run()

    #print (dir(out))
    print ('coeff',out.beta)
    print ('cov coeff0',out.cov_beta[0,:])
    print ('cov coeff1',out.cov_beta[1,:])
    print ('sd coeff', out.sd_beta)
    print ('res var', out.res_var)

    return out.beta,out.sd_beta,out.res_var

