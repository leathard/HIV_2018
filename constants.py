"""
Provides a place to set constatns for all various flavors of the simulation. Together with utility functions that are unlikely to change by their nature.
"""
import math


## partnership formation parameters used by Morris & Kretzschmar
## partnership formation parameters used by Morris & Kretzschmar
class mk:
    rho = 0.01
    sigma = 0.005


## factor by which to increase staged transmission parameters
## to acheive same R0 for serial monogamy
k = 1.435


## partnership parameters from Manicaland
class manic:
    rho = 0.007677224796228377
    sigma = 0.00195735077681542


# beta params are from EHG's code: concprimaryinf/R/conc.sim.R
# (because reported there with more precision than in EHG Table 1)
beta_p = 0.007315068
beta_a = 0.0002904110
beta_s = 0.002082192
beta_0 = 0
# dur params are from EGH's Table 1 (dur_0 differs from their posted code!)
dur_p = 88
dur_a = 3054
dur_s = 274
dur_0 = 307

beta_c = (beta_p * dur_p + beta_a * dur_a + beta_s * dur_s) / (dur_p + dur_a + dur_s + dur_0)
dur_c = dur_p + dur_a + dur_s + dur_0

# data from vandepitte and carael for subsaharan africa. 
# low value 1st quartile
# med value mean
# high value 3rd quartile 
# can add one extreme high?

# jk: betas for Immune agents Attia et al. 2008 
# jk: no data on different stages! 
# Attia et al. report 92% decrease in transmission rate 
beta_p_ART = .08 * beta_p
beta_a_ART = .08 * beta_a
beta_s_ART = .08 * beta_s

beta_c_ART = (beta_p_ART * dur_p + beta_a_ART * dur_a + beta_s_ART * dur_s) / (dur_p + dur_a + dur_s + dur_0)

beta_p_PREP = 0.4 * beta_p
beta_a_PREP = 0.4 * beta_a
beta_s_PREP = 0.4 * beta_s

p_nclients_low = 0.037
p_nclients_med = 0.087
p_nclients_high = 0.125
p_nsexworkers_low = 0.00065
p_nsexworkers_med = 0.017
p_nsexworkers_high = 0.0245

# jk: look for data 
p_nM_ART = 0.20
p_nF_ART = 0.20

#  11 different values of epsilon (the concurrency parameter)
#  EHG's R version: epsilon = 1-1.02*(exp(3.9*(seq(0,1,.1)-1))-.02)  
#  :note: lowest index is biggest (approx 1) epsilon
nOutGroups = 12  # based on EHG's Rsim.cpp; otherwise arbitrary
ehg_epsilons = tuple(1 - 1.02 * (math.exp(3.9 * (xi / 10. - 1)) - .02) for xi in range(nOutGroups - 1))


def phi_ehg(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that ``male`` and ``female`` are selected as potential partners.
    Based on individuals' current number of partnerships, and mixing parameter
    :note: epsilon is 1 if no concurrency and 0 if no resistance to concurrency
    :note: we will use closure to set the `epsilon` local to the `sim`
    """
    if (female.n_partners == 0 and male.n_partners == 0):
        return 1.0
    else:
        return (1.0 - epsilon)  # epsilon controls concurrency


def phi2(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that `male` and `female` are selected as potential partners.
    Based on each individual's current number of partnerships,
    and mixing parameter
    Sawers, Isaac, Stillwaggon - coital dilution example.
    :note:
      epsilon controls concurrency (1-> none, 0 -> no resistance),
      eF is female resistance to concurrency,
      eM is male resistance to concurrency,
      setting eF>eM (since epsilon in (0,1))
    """
    eF = 0 if female.n_partners == 0 else math.sqrt(epsilon)
    eM = 0 if male.n_partners == 0 else epsilon * epsilon
    return (1.0 - eF) * (1 - eM)  # epsilon controls concurrency


def phi3(male, female, epsilon):
    """Return float, the probability of partnership formation,
    **given** that `male` and `female` are selected as potential partners.
    Based on each individual's current number of partnerships,
    and mixing parameter ``epsilon``, which controls resistance
    to concurrency`
    note:
      epsilon controls concurrency (1-> none, 0 -> no resistance),
    jk: Here we add commercial sex workers with no
    it is the number of partners of the female that is controlling.
    Males will always accept a new partner who has no partners,
    and females without a partner will always accept a new partner.
    Otherwise partnerships form with prob 1-eps.
    """
    if (female.n_partners == 0):
        return 1.0
    else:
        return (1.0 - epsilon)  # epsilon controls concurrency


# jk: I do not understand why below averages are necessary?
def avtrans2asymtrans(beta, reltransfm):
    """Return dict, f2m->list and m2f->list.
    For average (stage-dependent) transmission rates beta,
    computes the f2m and m2f rates that give that average
    but have relative f2m transmission of reltransfm.

    Each transmission rate beta_i needs to be split into
    a M -> F and F -> M version so that two otherwise
    identical partnerships, one male infected and one
    female infected, will transmit at this rate on average.

    E.g., solving for reltransfm=0.7:
    b = (m + f)/2  and .7m = f
    yields
    b = 1.7m/2 or m= 2b/1.7

    beta : list of float
      the stage-dependent average transmission rates
    reltransfm : float
      female to male transmission rate relative to male to female

      #jk: fcsw to male transmission rate over different stages
      is the same as female to male transmission rate.
    """
    from copy import copy
    assert (reltransfm <= 1.0)

    mscale = 2. / (1 + reltransfm)
    m2f = tuple(mscale * bi for bi in beta)
    f2m = tuple(reltransfm * bi for bi in m2f)
    return {'m2f': m2f, 'f2m': tuple(reltransfm * bi for bi in m2f), 'fcsw2m': copy(f2m)}
