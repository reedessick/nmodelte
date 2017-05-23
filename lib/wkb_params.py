description = "storage of Nevin's WKB parameters for different models"
author      = "reed.essick@ligo.org, nevin@mit.edu"

#-------------------------------------------------

import nmode_utils as nmu

#-------------------------------------------------

### hard code in look-up tables for his estimates
wo = {
    0.5 : {
        1 : 1.46e-3,
        3 : 1.44e-3,
        5 : 1.43e-3,
        7 : 1.42e-3,
        9 : 1.41e-3,
    },
    0.6 : {
        1 : 1.18e-3,
        3 : 1.14e-3,
        5 : 1.12e-3,
        7 : 1.11e-3,
        9 : 1.10e-3,
    },
    0.8 : {
        1 : 0.90e-3,
        3 : 0.88e-3,
        5 : 0.86e-3,
        7 : 0.84e-3,
        9 : 0.82e-3,
    },
    1.0 : {
        1 : 0.71e-3,
        3 : 0.66e-3,
        5 : 0.60e-3,
        7 : 0.53e-3,
        9 : 0.43e-3,
    },
    1.1 : {
        1 : 0.60e-3,
        3 : 0.53e-3,
        5 : 0.44e-3,
    },
}

raise NotImplementedError, 'need to convert wo->(Rprim/Rsun)'
Rprim = {}

#---

### nevin's convention: w = a*l*wo/n
### my convention: w = alpha*l/n
nevin_a = {
    0.5 : {
        1 : 0.6,
        3 : 0.7,
        5 : 0.8,
        7 : 0.9,
        9 : 0.9,
    },
    0.6 : {
        1 : 1.1,
        3 : 1.3,
        5 : 1.4,
        7 : 1.5,
        9 : 1.7,
    },
    0.8 : {
        1 : 1.8,
        3 : 2.2,
        5 : 2.7,
        7 : 3.1,
        9 : 3.8
    },
    1.0 : {
        1 :  2.9,
        3 :  4.4,
        5 :  6.6,
        7 : 10.7,
        9 : 18.8,
    },
    1.1 : {
        1 :  4.0,
        3 :  7.4,
        5 : 14.2,
    },
}

raise NotImplementedError, 'need to convert a->alpha=a*wo'
alpha = {}

#---

c = {
    0.5 : {
        1 : 0.08e-12,
        3 : 0.08e-12,
        5 : 0.09e-12,
        7 : 0.10e-12,
        9 : 0.11e-12,
    },
    0.6 : {
        1 : 0.27e-12,
        3 : 0.30e-12,
        5 : 0.33e-12,
        7 : 0.37e-12,
        9 : 0.41e-12,
    },
    0.8 : {
        1 : 1.6e-12,
        3 : 1.9e-12,
        5 : 2.4e-12,
        7 : 3.3e-12,
        9 : 4.1e-12,
    },
    1.0 : {
        1 : 7.4e-12,
        3 : 12e-12,
        5 : 29e-12,
        7 : 120e-12,
        9 : 1200e-12,
    },
    1.1 : {
        1 : 18e-12,
        3 : 32e-12,
        5 : 790e-12,
    },
}

#---

Ialm_hat = {
    0.5 : {
        1 : 1.75e-2,
        5 : 1.50e-2,
        9 : 1.35e-2,
    },
    0.6 : {
        1 : 1.1e-2,
        5 : 1.0e-2,
        9 : 0.9e-2,
    },
    0.8 : {
        1 : 0.69e-2,
        5 : 0.55e-2,
        9 : 0.45e-2,
    },
    1.0 : {
        1 : 0.36e-2,
        5 : 0.21e-2,
        9 : 0.13e-2,
    },
    1.1 : {
        1 : 0.20e-2,
        5 : 0.10e-2,
    },
}

#--- 

### nevin's convention: kappa = kappa_hat*T*(P/day)**2
### my convention : kappa = k_hat*(T/0.2)*(P/10day)**2

nevin_k_hat = {
    0.5 : {
        1 : 0.03e3,
        3 : 1.0e3,
        5 : 1.7e3,
        7 : 2.2e3,
        9 : 4.1e3,
    },
    0.6 : {
        1 : 0.4e3,
        3 : 1.3e3,
        5 : 2.1e3,
        7 : 2.9e3,
        9 : 4.1e3,
    },
    0.8 : {
        1 : 0.8e3,
        3 : 2.9e3,
        5 : 4.8e3,
        7 : 6.9e3,
        9 : 11e3,
    },
    1.0 : {
        1 : 1.5e3,
        3 : 6.4e3,
        5 : 16e3,
        7 : 37e3,
        9 : 34e3,
    },
    1.1 : {
        1 : 1.4e3,
        3 : 16e3,
        5 : 21e3,
    },
}

raise NotImplementedError, 'need to convert from nevin\'s convention for k_hat to mine'
k_hat = {
}

#-------------------------------------------------

'''
masses = set()
ages   = set()
for d in [alpha_lookup, gamma_lookup, khat_lookup, Ialm_lookup, radius_lookup]:
    for mass in d.keys():
        masses.add( mass )
        for age in d[mass]:
            if not isinstance(age, str):
                ages.add( age )
masses = sorted(masses)
ages   = sorted(ages)
'''
