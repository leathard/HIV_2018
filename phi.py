"""Define probability of partnership `phi`. 

phi was defined by MK1996 and EHG and SIS implement. 
"""

def phi_ehg(male, female, epsilon):
    """Define the probability of partnership formation. 

    Parameter phi was defined by Morris&Kretzchmar 1996 paper and used by Eaton,
    Hallett, Garnett 2011 and Sawers, Isaac, Stillwaggon 2011. 
    This parameter is the probability of partnership formation given that male
    and female are selected as potential partners. It is based on individuals'
    current number of partnerships, and mixing parameter. 
    epsilon = 1 if no concurrency and 0 if no resistance to concurrency. 

    Returns: 
        the probability of partnership formation (float). 
    """
    if (female.n_partners == 0 and male.n_partners == 0):
        return 1.0
    else:
        return (1.0 - epsilon)  # epsilon controls concurrency


def phi2(male, female, epsilon):
    """Define the probability of partnership formation. 

    Parameter phi was defined by Morris&Kretzchmar 1996 paper and used by Eaton,
    Hallett, Garnett 2011 and Sawers, Isaac, Stillwaggon 2011. 
    This parameter is the probability of partnership formation given that male
    and female are selected as potential partners. It is based on individuals'
    current number of partnerships, and mixing parameter. 
    epsilon = 1 if no concurrency and 0 if no resistance to concurrency. 

    Sawers, Isaac, Stillwaggon - coital dilution example.
    :note:
      epsilon controls concurrency (1-> none, 0 -> no resistance),
      eF is female resistance to concurrency,
      eM is male resistance to concurrency,
      setting eF>eM (since epsilon in (0,1))

    Returns: 
        the probability of partnership formation (float). 
    """
    eF = 0 if female.n_partners == 0 else math.sqrt(epsilon)
    eM = 0 if male.n_partners == 0 else epsilon * epsilon
    return (1.0 - eF) * (1 - eM)  # epsilon controls concurrency


