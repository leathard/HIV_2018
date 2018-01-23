"""Provides parameter sets (as dicts) for simulations.

shared_params holds the common values across models.
These are picked to match EHG paper,
with three very slight differences in parameterization
(chosen for robustness to possible experiments).

- produce nMales males and nFemales females
- set pct of M and F seeded with HIV (matching EHG paper discussion),
  instead of an integer value
- set transmission rates for M and for F separately
"""
import math

def f2(params): #chk chk
    print params['sScenarioName']


## Parameters values used in all scenarios (picked to match EHG paper)
#    epsilon def is in sim.py

#WARNING: shared_params is not a complete parameter set!
#         Some vals are None, indicating they are to be filled in.
shared_params = dict(
    sScenarioName = '', #name of simulation model, and usu. folder for output
    nSim = 100,  #the number of replications for each parameter set (each value of epsilon)
    #equal number of males & females
    nMF = (10000,10000),
    pSeedHIV = (0.01,0.01),
    nBurnDays = 365*5,
    nSimDays=365*250,
    nOutputInterval = 365,
    nRndSeed = 7930881,  #fr EHG, but will be augmented based on sim specification
    #daily infection rate during each stage: a sequence
    beta_M2F = None,
    beta_F2M = None,
    #durations during each stage: a sequence
    durationsHIV = None,
    #aggregate partnership formation paramter
    rho = None,
    #chk this partnership dissolution param used by Partnership at initialization
    sigma01 = None,
    )

def phi_ehg(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that ``male`` and ``female`` are selected as potential partners.
    Based on individuals' current number of partnerships, and mixing parameter
    :note: epsilon measures resistence to concurrency;
      it is 1 if no concurrency and 0 if no resistance to concurrency
    :note: we will use closure to set the `epsilon` local to the `sim`
    """
    if ( female.n_partners==0 and male.n_partners==0 ):
        return 1.0
    else:
        return (1.0-epsilon)  #epsilon controls concurrency

"""
I'd like to combine this with targets.
EHG have epsilon as resistence to concurrency,
starting at 1.
That is like everyone having a target of 1.
Consider eps of 0.8.
That is a 20% chance of forming if at least one of them already has a partner.
Suppose (MK 2000) there is a 90% chance at at least one has a partner.
P(only man has partner) is .68*.27
P(only woman has partner) is .32*.73
P(both have) =.68*.73
So conditional on at least one having a partner,
20% is roughly the probability that the woman has no partner.

Here is the obvious way.
The target becomes a second constraint.
Probability is unchnaged if target is above current.
Probability falls toward zero if target is below zero.
So we have
- an overall concurrency penalty, e.
- M & F penalties for exceeding target, eF and eM.
- return (1-e)(1-eF)(1-eM)

This seems to have two parts:
1. selection of target (shifted geometric?)
2. determination of penalty
"""

def phi2(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that `male` and `female` are selected as potential partners.
    Based on each individual's current number of partnerships,
    and mixing parameter
    :note:
      epsilon controls concurrency (1-> none, 0 -> no resistance),
      eF is female resistance to concurrency,
      eM is male resistance to concurrency,
      setting eF>eM (since epsilon in (0,1))
    """
    eF = 0 if female.n_partners==0 else math.sqrt(epsilon)
    eM = 0 if male.n_partners==0 else epsilon*epsilon
    return (1.0-eF)*(1-eM)  #epsilon controls concurrency

def phi3(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that `male` and `female` are selected as potential partners.
    Based on each individual's current number of partnerships,
    and mixing parameter ``epsilon``, which controls resistance
    to concurrency
    :note:
      Could take a random draw inside phi, but gain efficiency
      by drawing them all at once in ``random_pairings``
    :note:
      epsilon controls concurrency (1-> none, 0 -> no resistance),
    Here we implement the following gender asymmetry:
    it is the number of partners of the female that is controlling.
    Males will always accept a new partner who has no partners,
    and females without a partner will always accept a new partner.
    Otherwise partnerships form with prob 1-eps.
    """
    if(female.n_partners==0):
        return 1.0
    else:
        return (1.0-epsilon)  #epsilon controls concurrency


## partnership formation parameters used by Morris & Kretzschmar 2000
class mk2000():
    rho = 0.00010564 #partnership formation (MUCH lower than earlier KM!)
    #duration parameters deduced from Table 1
    #mean duration in days is 1/sigma so mean duration in months is 1/(30 sigma)
    #so, sigma = 1/(30*duration in months)
    sigma01 = 1 / (30 * 239.1)   #their sigma1
    sigma02 = 1 / (30 * 28.4)  #their sigma2

## partnership formation parameters used by Morris & Kretzschmar
class mk():
    pass
mk.rho = 0.01 
mk.sigma = 0.005

## factor by which to increase staged transmission parameters
## to acheive same R0 for serial monogamy
k = 1.435

## partnership parameters from Manicaland
class manic():
    pass
manic.rho = 0.007677224796228377
manic.sigma = 0.00195735077681542


# beta params are from EHG's code: concprimaryinf/R/conc.sim.R
# (because reported there with more precision than in EHG Table 1)
beta_p = 0.007315068
beta_a =0.0002904110
beta_s = 0.002082192
beta_0 = 0
# dur params are from EGH's Table 1
# (I made a note that dur_0 differs from their posted R code,
# but now I'm not seeing where.  It matches *some* of their R code.)
dur_p = 88
dur_a = 3054
dur_s = 274
dur_0 = 307

#for simulations without staged transmission:
dur_c = dur_p+dur_a+dur_s+dur_0
beta_c = (beta_p*dur_p+beta_a*dur_a+beta_s*dur_s)/dur_c


#get the shared params
ehg_staged01 = shared_params.copy()
#params for the core EHG staged simulations
beta_ehg = (beta_p, beta_a, beta_s, beta_0)
#update the sim specific params
ehg_staged01.update(
    sim_name = 'ehg-staged',  #CHK do I use this? overridden in get_param_set
    rho = mk.rho,       #aggregate partnership formation paramter
    sigma01 = mk.sigma,   #chk where is this partnership dissolution param used?
    #daily infection rate during each stage: a sequence
    beta_M2F = beta_ehg,
    beta_F2M = beta_ehg,
    #durations during each stage: a sequence
    durationsHIV = (dur_p, dur_a, dur_s, dur_0),
    )

ehg_constant = ehg_staged01.copy()
ehg_constant['sim_name'] = 'ehg-constant',
ehg_constant['beta_M2F'] = beta_c, 0, 0, 0
ehg_constant['beta_F2M'] = beta_c, 0, 0, 0
ehg_constant['durationsHIV'] = dur_c, 0, 0, 0

mk_constant = ehg_constant.copy()
mk_constant['beta_M2F'] = 0.05, 0, 0, 0
mk_constant['beta_F2M'] = 0.05, 0, 0, 0

ehg2sigma_constant = ehg_constant.copy()
ehg2sigma_constant['sigma01'] = mk2000.sigma01
ehg2sigma_constant['sigma02'] = mk2000.sigma02

ehg2sigma_staged = ehg_staged01.copy()
ehg2sigma_staged['sigma01'] = mk2000.sigma01
ehg2sigma_staged['sigma02'] = mk2000.sigma02

def avtrans2asymtrans(betas, reltransfm):
    """Return dict, f2m->list and m2f->list.
    For average (stage-dependent) transmission rates betas,
    computes the f2m and m2f rates that give that average
    but have relative f2m transmission of reltransfm.

    Each transmission rate beta_i needs to be split into
    a M -> F and F -> M version so that two otherwise
    identical partnerships, one male infected and one
    female infected, will transmit at this rate on average.

    E.g., solving for rfm=0.7:
    b = (m + f)/2  and .7m = f
    yields
    b = 1.7m/2 or m= 2b/1.7

    betas : [float]
      the stage-dependent average transmission rates
    reltransfm : [float]
      female to male transmission rate relative to male to female
    """
    result = dict()
    mscales = tuple(2./(1.+r) for r in reltransfm)
    m2f = tuple(mscale*bi for (mscale,bi) in zip(mscales,betas))
    result['m2f'] = m2f
    result['f2m'] = tuple(r*bi for (r,bi) in zip(reltransfm,m2f))
    return result

# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 autoindent
