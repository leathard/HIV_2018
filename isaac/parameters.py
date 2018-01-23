"""
Create scenarios (parameter dicts) for concurrency simulations.
IN PROGRESS
"""
import ConfigParser
import glob, math, os, random
import numpy as np

def phiEHG(male, female, epsilon):
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

class Scenarios(object):
    """Provides a class for representing the baseline
    (as defined in the file `fbaseline`)
    and associated scenarios
    (as defined in the file `fscenarios`).
    :note: output filenames are not part of a scenario
    """
    _baseline = None
    def __init__(self, fbaseline, fscenarios):
        """Return None.
        :side-effects: set `self.baseconfig` and `self.parsedScenarios`
          each to a `ConfigParser`.

        fbaseline : str
          filename for baseline parameters
        fscenarios : str
          filename for scenarios (deviations from baseline)
        """
        # get the baseline parameters (parameter names are the *sections*,
        #  parameter values are the "value" *option*)
        self.baseconfig = baseconfig = ConfigParser.ConfigParser()
        with open("baseline.ini","r") as fin:
            baseconfig.readfp(fin)
        # get the scenarios from `fscenarios` (scenario names are the *sections*;
        #  parameter names are the *options*)
        self.parsedScenarios = parsedScenarios = ConfigParser.ConfigParser()
        parsedScenarios.optionxform = str #otherwise options will be lower-cased
        with open(fscenarios,"r") as fin:
            parsedScenarios.readfp(fin)

    def __str__(self):
        """Return [str], the scenario names."""
        return "\n".join(s for s in self.parsedScenarios.sections()) 

    @property
    def baseline(self):
        """Return dict, the baseline scenario (from `fbaseline`).
        Note that baseconfig sections are parameter names.
        """
        if self._baseline is None:
            baselineParameters = dict(sScenarioName="baseline")
            baseconfig = self.baseconfig #get parsed fbaseline
            #each baseconfig section is a parameter
            for section in baseconfig.sections():
                ptype = baseconfig.get(section,"type")
                #might list as parameters items whose values are computed from parameters;
                #save these computations for later
                if baseconfig.get(section,"value") == "computed":
                    continue
                if ptype == "float":
                    val = baseconfig.getfloat(section,"value") 
                elif ptype == "int":
                    val = baseconfig.getint(section,"value")
                    assert isinstance(val, int) #chkchk superfluous
                elif ptype == "bool":
                    val = baseconfig.getboolean(section,"value")
                #may need to evaluate some parameters (DANGEROUS: evaluates arbitrary code!!!)
                else:
                    val = baseconfig.get(section,"value")
                    val = self.cast_config_val(section, val) #cast val to specified type
                baselineParameters[section] = val
            computeParams(baselineParameters) #chkchk
            self._baseline = baselineParameters
        return self._baseline


    def get_scenario(self, scenarioName):
        """Return dict, the parameter configuration for one scenario,
        as specified first in `fbaseline` and then in `fscenarios`.
        The scenario specifications in `fscenarios` are only
        parameter values that deviate from the baseline.
        Thus we being with (a copy of) the baseline, and modify it.

        scenarioName : str
          the name of the scenario (i.e., the section in `fscenarios`)
        :raises: AssertionError
        """
        parsedScenarios = self.parsedScenarios
        assert scenarioName in parsedScenarios.sections(), \
            "{} not available in scenarios.ini".format(scenarioName)
        #start with the baseline parameters:
        scenario = self.baseline.copy()
        if (scenarioName == "baseline"):
            return scenario
        scenario["sScenarioName"] = scenarioName
        #scenario["name"] = scenarioName  #chkchkchk
        #get the settings for this scenario, and update `scenario`
        for param, val in parsedScenarios.items(scenarioName):
            val = self.cast_config_val(param, val) #cast val to specified type
            scenario[param] = val
        computeParams(scenario) #chkchkchk
        return scenario

    def cast_config_val(self, name, val):
        """Return obj, of type determined by `fbaseline`.
        name : str
          the parameter name
        val : str
          a string representation of the value
        """
        basecfg = self.baseconfig
        ptype = basecfg.get(name,"type")
        #the first cases just check for correct type
        if ptype == "float":
            if not isinstance(val, float): val = float(val)
        elif ptype == "int":
            if not isinstance(val, int): val = int(val)
        elif ptype == "bool":
            if not isinstance(val, bool):
                val = False if (val.lower()=="false") else True
        #the remaining cases handle casting of string to value
        elif ptype == "intEval":
            val = eval(val)
            assert isinstance(val, int), "{} value {} is not int".format(section,val)
        elif ptype == "floatEval":
            val = eval(val)
            assert isinstance(val, float), "{} value {} is not float".format(section,val)
        elif ptype == "[float]":
            val =  tuple(float(x) for x in val.split())
            for v in val:
                assert isinstance(v, float), "{} value {} is not float".format(section,val)
        elif ptype == "[int]":
            val =  tuple(int(x) for x in val.split())
            for v in val:
                assert isinstance(v, int), "{} value {} is not float".format(section,val)
        #otherwise, keep string value
        elif ptype == "computed":
            val =  None
        else:
            assert isinstance(val, str), "{} value {} is not str".format(section,val)
        return val

def computeParams(pdict):
    """Return None.
    Adds to pdict the values of computed parameters.
    :side-effects: mutates pdict
    """
    #move eps computation into here? #chkchkchk

    #possibly lower duration (higher sigma) secondary partnerships
    pSigma02 = pdict['pSigma02']
    if pSigma02 < 1.0:
        raise ValueError('pSigma02: {}'.format(pSigma02))
    if (pdict['sScenarioName'] == "baseline"):
        assert pdict['sigma02'] is None, "{}".format(pdict['sigma02'])
    pdict['sigma02'] = pSigma02 * pdict['sigma01']

    #possibly asymmetric HIV transmission rates:
    avtrans = pdict['transmissionRatesHIV']
    reltrans = pdict['relativeTransmissionRatesHIV']
    asymtrans = avtrans2asymtrans(avtrans, reltrans)
    if (pdict['sScenarioName'] == "baseline"):
        assert pdict['beta_M2F'] is None
    pdict['beta_M2F'] = asymtrans['m2f']
    if (pdict['sScenarioName'] == "baseline"):
        assert pdict['beta_F2M'] is None
    pdict['beta_F2M'] = asymtrans['f2m']

    sPhi = pdict['sPhi']
    if (sPhi == 'phiEHG'):
        pdict['model_phi'] = phiEHG  #in this file
    else:
        raise ValueError('unknown phi name')

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



scenarios = Scenarios("baseline.ini", "scenarios.ini")
baseline = scenarios.baseline
snames = """duration0100 duration0200""".split()
for name in snames:
    exec("{name} = scenarios.get_scenario('{name}')".format(name=name))
#just testing that naming is working as expected
for name in snames:
    assert name == eval(name)['sScenarioName']



def get_param_set(params, n_sim, overwrite=True):
    """Yield dict, each a replicate specific set of parameters.
    There are `n_sim` (e.g., 100) replicates for the parameters `params`.
    Each replicate is given its own random seed, sim_num, and outfile names.
    This generator function should yield dicts with **pickleable** elements only,
    so that it can be used with multiprocessing.
    
    params : dict
      the scenario to be simulated
    n_sim : int
      the number of replicates
    overwrite : bool
      whether or not to overwrite previous results
    """
    name = params['sScenarioName']
    #all output goes below `out` folder, which should already exist
    outdir = os.path.join("out",name)
    if not os.path.exists(outdir):
        print "creating output directory"
        os.mkdir(outdir) 
    #"old" and "new" means without and with modern sector
    outfiles = ('oldmacro', 'newmacro', 'oldmicro', 'newmicro')
    previously_done = list()
    if not overwrite:
        for outname in outfiles:
            outfile_pattern = outname + '???.csv'
            search_pattern = os.path.normpath(os.path.join(outdir, outfile_pattern))
            previously_done.extend(glob.glob(search_pattern))
    previously_done = set(int(x[-7:-4]) for x in previously_done) #just the replicate numbers
    for n in range(n_sim):
        if n in previously_done:
            print 'Skip: simulation {0} (already completed).'.format(n)
            continue
        else:
            simparams = params.copy()  #fresh copy for each replicate (IMPORTANT!)
            simparams['sim_num'] = n
            seed = (params['nRndSeed'], n)   #replicate(n)-specific seed
            simparams['prng'] = random.Random(seed)
            loggername =  "{0}replicate{1:03d}".format(name,n) #replicate specific logger name
            simparams['loggername'] = loggername
            simparams['logfile'] = os.path.normpath(os.path.join(outdir, loggername+".txt"))
            for outfile in outfiles: #create replication specific output files
                fullname = "{0}{1:03d}.csv".format(outfile, n)
                simparams[outfile] = os.path.normpath(os.path.join(outdir, fullname))
            yield simparams

def get_lhs_paramset(baseparams, limitdict, nsamples, overwrite=True, prng=None):
    """Yield dict, a simulation specific set of parameters.
    There will be len(limitdict) dicts of parameters `params`.
    Each replicate has its own random seed, sim_num, and outfile names.
    This generator function should yield dicts with **pickleable** elements only,
    so that it can be used with multiprocessing.
    """
    if prng is None:
        prng = np.random.RandomState(314) #used by `lhs`
    exp_name = baseparams['name']
    outdir = os.path.join('out',exp_name)
    outfiles = ('oldmacro', 'newmacro', 'oldmicro', 'newmicro')
    """
    previously_done = list()
    for outname in outfiles:
        outfile_pattern = outname + '???.csv'
        search_pattern = os.path.normpath(os.path.join(outdir, outfile_pattern))
        previously_done.extend(glob.glob(search_pattern))
    if not overwrite:
        previously_done = set(int(x[-7:-4]) for x in previously_done)
    elif previously_done:
        doremove = raw_input("Remove files from last experiment of {}?(y,n)".format(exp_name))
        if doremove[0].lower()=='y':
            for fname in previously_done:
                os.remove(fname)
            previously_done = []
    """
    pnames, plims = zip(*limitdict.items())
    design = lhs(nsamples, plims, midpt=False, prng=prng)
    for n, row in enumerate(design):
        simparams = baseparams.copy()  #fresh copy for each sim (IMPORTANT!)
        simparams.update(zip(pnames,row))
        simparams['nAgents'] = int(round(simparams['nAgents'])) #restore int value
        simparams['sim_num'] = n
        seed = (baseparams['nRndSeed'], n)   #replicate(n)-specific seed
        simparams['prng'] = random.Random(seed)
        simparams['logfile'] = os.path.normpath(os.path.join(outdir, 'logfile.txt'))
        simparams['loggername'] = '{0}{1:03d}'.format(exp_name,n)
        for outfile in outfiles:
            fullname = "{0}{1:03d}.csv".format(outfile, n)
            simparams[outfile] = os.path.normpath(os.path.join(outdir, fullname))
        yield simparams



#for each param, divide range into n intervals, and use each interval *once*
def lhsdesign(n, p, midpt=False, prng=None):
    """Return array. Each of the `n` rows
    has `p` values.

    Parameters
    ----------
    n : int
      number of samples
    p : int
      the number of "parameters"
    midpt : bool
      Use midpoints if true. Else sample randomly
      from each interval.
    prng : RandomState
      A NumPy random number generator.
    """
    if prng is None:
        prng = np.random.RandomState()
    tpdesign = np.empty((p,n), dtype=float)
    for idx in range(p):
        #left bound
        smpl = np.arange(n, dtype=float)
        #add portion of interval
        if midpt:
            smpl += 0.5
        else:
            smpl += prng.random_sample(n)
        #scale
        smpl /= float(n)
        #shffle
        prng.shuffle(smpl)
        tpdesign[idx] = smpl
    return tpdesign.transpose()

def lhs(n, lohiseq, midpt=False, prng=None):
    """Return array. Each of the `n` rows
    has a value for each of the (len(lohiseq) parameters,
    whose ranges are provided by lohiseq.  

    Parameters
    ----------
    n : int
      number of samples
    lohiseq : sequence of 2 tuples
      the (low,high) range for each parameter
    midpt : bool
      Use midpoints if True. Else sample randomly
      from each interval.
    prng : RandomState
      A NumPy random number generator.
    """
    if prng is None:
        prng = np.random.RandomState()
    p = len(lohiseq)
    design = lhsdesign(n, p, midpt, prng)
    for idx, row in enumerate(lohiseq):
        low, high = row
        delta = (high-low)
        #scale
        design[:,idx] *= delta
        #translate
        design[:,idx] += low
    return design

if __name__ == '__main__':
    scenarios = Scenarios("baseline.ini", "scenarios.ini")
    for nm in str(scenarios).split():
        print(nm)
        s = scenarios.get_scenario(nm)
        for attr in ("sScenarioName", "nSim", "nMF", "pSeedHIV", "nBurnDays",
            "nSimDays", "nOutputInterval", "nRndSeed", "beta_M2F", "beta_F2M",
            "durationsHIV",  "rho",  "sigma01", "pSigma02"):
            assert attr in s, "{attr} not in {s}".format(attr=attr, s=s)



# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 autoindent
