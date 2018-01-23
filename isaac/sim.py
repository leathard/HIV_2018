"""Provides function ``onesim``, which can run a single concurrency simulation.
Compare to EHG's sim.cpp.
"""
import numbers
import multiprocessing as mp
import functools, glob, math, shutil, tempfile, time, pprint
import os.path as path
from itertools import izip as zip   #Python 2.7
from collections import defaultdict
import datetime, logging
#logging.basicConfig(level=logging.INFO)
#setting the level to CRITICAL effectively turns off logging
logging.basicConfig(level=logging.CRITICAL,filename='out/temp.log', filemode='w')
logging.info(datetime.datetime.today().isoformat())
import numpy as np
from numpy.random import RandomState

## imports from concurrency.py
from concurrency import Partnership, Person, stagedHIVfactory

##########  globals (globals should **only** be constants! for MP!)
#11 different values of epsilon (the concurrency parameter)
#  EHG's R version: epsilon = 1-1.02*(exp(3.9*(seq(0,1,.1)-1))-.02)  
#  :note: lowest index is biggest (approx 1) epsilon
ehg_epsilons =  tuple( 1 - 1.02*(math.exp(3.9*(xi/10. - 1))-.02) for xi in range(11) )
nOutGroups = 12  #based on EHG's Rsim.cpp; otherwise arbitrary

def random_pairings(n_pairs, males, females, params):
    """Return n_pairs random male-female pairs (**with** replacement).

    n_pairs : int
        the number of pairs to form
    males : list
        the potential male partners
    females : list
        the potential female partners
    params : dict
        params['prng'] is the random number generator;
        params['phi'] is a function that returns partnership formation probability 
    """
    _phi = params['sim_phi']  #e.g., phiEGH in parameters.py
    _prng = params['prng']
    _randint = _prng.randint
    _random_sample = _prng.random_sample
    if n_pairs < 0:
        raise ValueError('cannot form a negative number of pairs')
    nMales = len(males)
    nFemales = len(females)
    assert n_pairs < nMales and n_pairs < nFemales
    while (n_pairs > 0):
        # `n_pairs` random males and females (*with* replacement)
        males2pair = (males[idxM] for idxM in _randint(0, nMales, n_pairs))
        females2pair = (females[idxF] for idxF in _randint(0, nFemales, n_pairs))
        draws = _random_sample(n_pairs)
        #must form pairs sequentially due to concurrency resistence
        for (male, female, draw) in zip(males2pair, females2pair, draws):
            # form partnership with probability phi, which varies with concurrency
            if(draw < _phi(male, female)): #ai: EHG use binomial; why? chk
                #ai: check that partnership doesn't already exist
                if female not in male.partners:
                    assert male not in female.partners
                    n_pairs -= 1
                    yield male, female


class Scheduler(object):
    """Provides a scheduler object that
    maintains a schedule of future events
    and also keeps a registry of other objects.
    """
    def __init__(self, params):
        #TODO: change to weakrefs
        self.params = params
        self.deaths = defaultdict(list)  # multimap for scheduling deaths: map date -> list of Person
        self.cures = defaultdict(list)  # multimap for scheduling cures: map date -> list of Person
        self.transmissions = defaultdict(list)   # map **upcoming** transmission date to list of Partnership
        self.dissolutions = defaultdict(list)  # map each **end** date list of Partnership
        #list of objects (switch to weakref set?)
        self.n_people = 0 # initialize partnership counter
        self.males = list()
        self.females = list()
        self.n_partnerships = 0 # initialize partnership counter
        self.partnerships = list()
    def coresim(self, males, females, day, params):
        """Return None; modify `params`.
        Run the core stages of the simulation (one iteration):
        form partnerships, do HIV transmissions, 
        dissolve partnerships, do HIV deaths.
        These are run every iteration (including burn days).  
        """
        self.form_partnerships(males, females, day)
        self.hiv_transmissions(day)
        self.dissolve_partnerships(day)
        self.hiv_deaths(day)
    def register_person(self, person):
        if person.sex == 'M':
            self.males.append(person)
            self.n_people += 1
        elif person.sex == 'F':
            self.females.append(person)
            self.n_people += 1
        else:
            raise ValueError('unknown sex {}'.format(person.sex))
        assert self.n_people==len(self.males)+len(self.females)
    def register_infection(self, disease):
        end = disease.end_date
        host = disease.host
        if disease.is_fatal:
            if host in self.deaths[end]:
                raise ValueError('disease already registered')
            self.deaths[end].append(host)
        else:
            raise ValueError('all currently fatal')
            #self.cures[end].append(host)
    def register_transmission(self, partnership):
        self.transmissions[partnership.transmission_date].append(partnership)
    def register_partnership(self, pship):
        # add partnership to partnership multimap with end date
        self.partnerships.append(pship)
        self.n_partnerships += 1
        assert self.n_partnerships==len(self.partnerships),\
            '{} != {}'.format(self.n_partnerships,len(self.partnerships))
        self.dissolutions[pship.end_date].append(pship)
    def delete_partnership(self, pship):
        """Return None; delete a partnership from ``partnerships`` multimap
        *before* scheduled end date (e.g., when partner dies from HIV).
        Also remove any scheduled transmission.
        """
        self.dissolutions[pship.end_date].remove(pship)
        self.partnerships.remove(pship)
        self.n_partnerships -= 1
        #destroy the partnership (match EHG's destructor for the partnership class)
        male, female = pship._mf
        male.remove_partnership(pship) #removes from set *and* decrements this person's count
        female.remove_partnership(pship)
        # check if a future transmission is scheduled, and if so remove it
        if(pship.transmission_scheduled):
            self.transmissions[pship.transmission_date].remove(pship)

    def clear_partnerships(self):
        """Return None. Delete all partnerships (preparing for next simulation).
        It's just cleanup; compare EHG's DumpPartnerships #chk 
        This may work very differently than the EHG code, but that shd't matter.
        """
        for lst in self.dissolutions.values():
            for pship in lst[:]:
                self.delete_partnership(pship)
        assert self.n_partnerships == 0
        self.dissolutions.clear()
        self.transmissions.clear()

    ######################  ai: convenience methods #######################################
    def count_partnerships(self):
        return sum(len(lst) for lst in self.dissolutions.values())
    def count_transmissions(self):
        return sum(len(lst) for lst in self.transmissions.values())
    def count_scheduled_deaths(self):
        return sum(len(lst) for lst in self.deaths.values())
    def show_one_day(self, day):
        """Return str, the schedule (so far) for `day`.
        (Used only for logging.)
        """
        males = self.males
        females = self.females
        report = """
        Schedule for day {day}
        ======================
        Number infected: {ninf} of {pop}
        Transmissions: {ntrans_total} total and {ntrans_today} today
        Deaths: {d_total} total and {d_today} today
        Dissolutions: {p_total} total and {p_today} today
        """.format(
        day=day,
        ninf=sum(p.is_infected for p in (males+females)),
        pop=len(males)+len(females),
        ntrans_total=sum(len(lst) for lst in self.transmissions.values()),
        d_total=sum(len(lst) for lst in self.deaths.values()),
        p_total=sum(len(lst) for lst in self.dissolutions.values()),
        #scheduled for today (BUT: avoid adding keys!)
        ntrans_today= len(self.transmissions.get(day,[])),
        d_today= len(self.deaths.get(day,[])),
        p_today= len(self.dissolutions.get(day,[])),
        )
        return report

    #################   MAIN SCHEDULE   ##################################
    def form_partnerships(self, males, females, day): #chkchkchk schedule partnerships
        params = self.params
        ##### form partnerships
        logging.debug('coresim: BEGIN form partnerships')
        old_partnership_count = self.n_partnerships   #chkchkchk
        assert self.count_partnerships() == old_partnership_count
        # population sizes
        nMales = len(males); nFemales = len(females);
        #We follow EHG and impose a hard ceiling on the number of partnerships!
        # (There a lots of things not to like about this.)
        # There is a maximum of about one partnership per person. (vs MK2000)
        max_new_partnerships = (nMales+nFemales)//2 - self.n_partnerships
        assert type(max_new_partnerships) is int
        # determine how many new partnerships to form (binomial draw with probability rho)
        prng = params['prng']
        rho = params['rho']     #aggregate partnership formation paramter
        assert type(rho) is float
        nFormPartnerships = 0  # needed to match EHG behavior (R's rbinom(0, rho)=0)
        if max_new_partnerships > 0:
            nFormPartnerships = prng.binomial(n=max_new_partnerships, p=rho)
        # form exactly nFormPartnerships new partnerships
        npairings = 0  #counter for testing that `random_pairings` returns the right number
        #:note: random pairings uses phi, which is in `params`
        for male, female in random_pairings(nFormPartnerships, males, females, params):
            newpartnership = Partnership(male, female, day, registry=self, params=params)  # form partnership
            npairings += 1
        assert npairings == nFormPartnerships
        logging.debug('coresim: END form partnerships')
        logging.info('\t{0} partnerships formed'.format(nFormPartnerships))
        assert self.count_partnerships() == old_partnership_count + nFormPartnerships

    def hiv_transmissions(self, day):
        params = self.params
        ##### do HIV transmissions
        logging.debug('coresim: begin HIV transmission')
        #recall: schedule.transmissions maps days to lists of partnerships
        Disease = params['Disease']
        transmissions4today = self.transmissions[day]  #list of partnerships
        ntrans_today = len(transmissions4today)
        ntrans_waiting = self.count_transmissions()
        logging.info('do {0} of {1} transmissions'.format(ntrans_today, ntrans_waiting))
        for pship in transmissions4today:
            transmitter, disease = pship.transmit(Disease, day) #chk is there a better way?
            if transmitter is not None:
                assert disease is not None
                tally_transmission(day=day, transmitter=transmitter, counter=params['counters'])
        del self.transmissions[day]
        del transmissions4today
        # lose ntrans_today future transmissions but may have gotten some new via infections at transmission
        assert self.count_transmissions() >= ntrans_waiting - ntrans_today
        logging.debug('coresim: end HIV transmission')

    def dissolve_partnerships(self, day):
        ##### dissolve partnerships
        logging.debug('begin: do scheduled partnership dissolutions')
        partnerships4today = self.dissolutions[day]
        logging.info('Do {0} dissolutions'.format(len(partnerships4today)))
        for pship in list(partnerships4today):  #*termination* dates map to partnership lists
            self.delete_partnership(pship)
        assert not self.dissolutions[day]  #shd be empty
        del self.dissolutions[day]
        del partnerships4today
        logging.debug('end: do scheduled partnership dissolutions')

    def hiv_deaths(self, day):
        ##### do HIV deaths
        logging.info('begin: HIV deaths')
        deaths4today = self.deaths[day]
        logging.info('do {0} deaths'.format(len(deaths4today)))
        assert len(deaths4today)==len(set(deaths4today))
        for indiv in tuple(deaths4today): #don't really need a copy (chk Person.die)
            self.kill(indiv)
        del self.deaths[day]
        del deaths4today
        logging.info('end: HIV deaths')

    def kill(self, indiv):  #compare EHG's `Kill`
        """Return None; "kills" the individual (i.e., reset, so population is constant).
        Terminate partnerships and reset individual to unpartnered susceptible.
        """
        logging.info('ENTER: kill indiv')
        #terminate all partnerships
        for pship in tuple(indiv.partnerships):
            self.delete_partnership(pship)
        assert len(indiv.partnerships)==0
        indiv.reset()  #reborn anew, uninfected, same sex
        logging.info('EXIT: kill indiv')



def run_simulations_mp(params, workerpool=None):
    if workerpool is None:
        n_processors = mp.cpu_count()
        print '{0} cores available'.format(n_processors)
        workerpool = mp.Pool(n_processors)
    workerpool.map(onesim, get_param_set(params))
    #for worker in workerpool: worker.join()
    return workerpool

def run_simulations(params):
    """Return None. Run all pending simulation for this model.
    Input should be a dict of model-specific parameters.
    """
    for simparams in get_param_set(params):
        onesim(simparams)

def get_param_set(params):
    """Yield dict, a replicate specific set of parameters.
    There are nSim (e.g., 100) replicates for each epsilon,
    and each one has its own random seed, sim_name, and outfilename.
    This generator function should yield dicts with **pickleable** elements only,
    so that it can be used with multiprocessing.
    """
    outfolder = params['outfolder']
    outfile_pattern = 'eps??sim??.out'
    search_pattern = path.normpath(path.join(outfolder, outfile_pattern))
    previously_done = glob.glob(search_pattern)
    for n in range(params['nSim']):
        for i,eps in enumerate(ehg_epsilons):  #ehg_epsilons is a global constant
            assert type(eps) is float
            sim_name = 'eps{0:02d}sim{1:02d}'.format(i,n)
            outfilename = path.normpath(path.join(outfolder, sim_name + '.out'))
            if path.normpath(outfilename) in previously_done:
                print 'Skip simulation ' + sim_name + ' (already completed).'
                continue
            else:
                seed = params['nRndSeed'], n, int(eps * 10**5)   #replicate(n)-specific seed
                simparams = params.copy()  #fresh copy for each replicate (IMPORTANT!)
                simparams.update(
                    prng = RandomState(seed=seed),
                    #get phi for this simulation (to be passed to `random_pairings` via `coresim`)
                    epsilon = eps,  # `phi` needs this!
                    sim_name = sim_name,
                    outfilename = outfilename,
                    )
                yield simparams

    

def old_run_simulations(params):
    """Return None.
    Run all the simulations for a model (determined by `params`).
    """
    for n in range(params['nSim']):
        for i,eps in enumerate(ehg_epsilons):
            assert type(eps) is float
            sim_name = 'eps{0:02d}sim{1:02d}'.format(i,n)
            outfolder = params['outfolder']
            outfile_pattern = params['outfile_pattern']
            #params['outfile_pattern'] = 'eps??sim??.out'
            outfilename = outfolder + sim_name + '.out'
            previously_done = glob.glob(path.normpath(outfolder+outfile_pattern))
            if path.normpath(outfilename) in previously_done:
                print 'Skip simulation ' + sim_name + ' (already completed).'
                continue
            else:
                print '\nDo simulation ' + sim_name
                #give every run a unique seed
                seed = [params['nRndSeed'], n, int(eps * 10**5)]
                params['prng'] = RandomState(seed=seed)
                #get phi for this simulation (to be passed to `random_pairings` via `coresim`)
                params['epsilon'] = eps  # `get_phi` needs this!
                params['phi'] = params['get_phi'](params)
                #create a temp file to hold results during one simulation
                tempfh = tempfile.NamedTemporaryFile(mode='w', suffix='.out', dir=outfolder, delete=False)
                tempfname = tempfh.name
                params['fout'] = tempfh
                #run one simulation
                t0 = time.clock()
                onesim(params=params)
                dt = time.clock() - t0
                print '\nSimulation ' + sim_name + ' completed successfully in {0:d} minutes'.format(int(dt)//60)
                print 
                tempfh.close()
                shutil.move(tempfname, outfilename)


def onesim(params):
    """Return None; setup and run one replicate
    (e.g., all daily iterations for 1 out of 100 replicates, for 1 value of epsilon)
    """
    sim_name = '{0}: {1}'.format(params['sScenarioName'], params['sim_name'])
    print '\n\nBegin ' + sim_name + '\n' + pprint.pformat(params)
    #create a temp file to hold the results of one simulation
    outfolder = params['outfolder']
    tempfh = tempfile.NamedTemporaryFile(mode='w', suffix='.out', dir=outfolder, delete=False)
    assert params.get('fout', None) is None
    params['fout'] = tempfh   #used by `record_output`
    tempfname = tempfh.name   #we'll rename this if the simulation runs to completion
    #time the simulation
    t0 = time.clock()
    schedule = Scheduler(params=params)  #event scheduler (transmissions, deaths, dissolutions)
    nMales, nFemales = params['nMF']
    nSimDays = params['nSimDays']
    nBurnDays = params['nBurnDays']
    nOutputInterval = params['nOutputInterval']  # usu. 365 (i.e., once a year)
    # type check parameters
    assert type(nMales) is int
    assert type(nFemales) is int
    assert type(nSimDays) is int
    assert type(nBurnDays) is int
    assert type(nOutputInterval) is int

    #create the disease for this model
    # note: don't add factory to model_params as it won't pickle,
    #       which create MP problems
    params['Disease'] = stagedHIVfactory(durations=params['durationsHIV'],
                                    transM2F=params['beta_M2F'],
                                    transF2M=params['beta_F2M'])

    #create the phi for this model
    # note: don't add to model_params as it won't pickle -> MP problems
    #params['sim_phi'] = functools.partial(params['model_phi'], epsilon=params['epsilon'])
    params['sim_phi'] = lambda male, female: params['model_phi'](male, female, params['epsilon'])

    #counter used to tally incidence transmission by stage of infection
    # note: don't make this global (MP!)
    params['counters'] = counter = defaultdict(int)

    #prepare outfile by writing header
    header = 'nMinfect,nFinfect,'
    header += 'MPrimaryTrans,FPrimaryTrans,MAsymptomaticTrans,FAsymptomaticTrans,MSymptomaticTrans,FSymptomaticTrans,'
    header += 'MPships,,,,,,,,,,,,FPships,,,,,,,,,,,,iMPships,,,,,,,,,,,,iFPships,,,,,,,,,,,,'
    header += 'MPrimary,FPrimary'
    tempfh.write(header)

    # Prepare variables (note that `males` and `females` not reused across simulations)
    males = list(Person(sex='M', registry=schedule, params=params) for i in range(nMales))
    females = list(Person(sex='F', registry=schedule, params=params) for i in range(nFemales))

    #begin simulation loop /* Do the simulations */
    for day in range(nSimDays+nBurnDays):
        logging.info('\nBEGIN ITERATION for day {0}.'.format(day))
        logging.debug(schedule.show_one_day(day))

        # Seed infections after nBurnDays have passed (ONCE)
        if(day == nBurnDays):
            assert schedule.count_scheduled_deaths() == 0
            diseases = seedInfections(males, females, day, schedule=schedule, params=params)
            assert schedule.count_scheduled_deaths() == len(diseases)  #for now, disease is fatal
            assert all(deathday >= day for deathday in schedule.deaths) #equal only possible on day==nBurnDays

        #run the core of the simulation (runs even during burn days)
        #params holds counter and fout
        schedule.coresim(
            males=males,
            females=females,
            day=day,
            params=params
            )

        # Record the output once a "period" (i.e., every nOutputInterval days)
        # :note: this won't record after last year is run (it needs one more day to pass the test).
        #        We keep it this way just to match EHG.
        if( day >= nBurnDays and (day-nBurnDays) % nOutputInterval == 0 ):
            print '.',
            outIndex = (day - nBurnDays) / nOutputInterval
            #ai: recording strategy
            #    params holds counter and fout
            record_output(males, females, params)
            #reset counter
            counter.clear()  #reuse counter each day (since wrote contents to disk)

    #END of simulation; just need to clean up: reset static elements
    #ai: mostly don't need to clean up in Python unless we intend to reuse objects
    # (and EHG don't even reuse the persons, so the rest of object creation is a trivial cost)
    # prepare classes for reuse
    schedule.clear_partnerships() # clears the `partnerships` and `transmissions` multimaps
    schedule.deaths.clear()

    tempfh.close()
    dt = time.clock() - t0
    #since the simulation ran to completion, we can use the output file
    outfilename = params['outfilename']
    shutil.move(tempfname, outfilename)
    msg = """
    {0} completed successfully in {1:d} minutes.
    {0} output written to {2}.
    """.format(sim_name, int(dt)//60, outfilename)
    logging.info(msg)





#### add some helper functions


def record_output(males, females, params):
    """Return None. Write data to disk.
    :note: EHG store these values in arrays; we write them to disk instead
    :note: we record the same values (in same order) as EHG,
      then we append some additions
    """
    counter = params['counters']  #count transmissions by stage
    fout = params['fout']  #the output file handle
    n_infected_males = sum(m.is_infected for m in males)
    n_infected_females = sum(f.is_infected for f in females)
    #build up this period's data as a list
    data = [n_infected_males, n_infected_females]
    data += [
            counter['malePrimaryTransToday'],
            counter['femalePrimaryTransToday'],
            counter['maleAsymptomaticTransToday'],
            counter['femaleAsymptomaticTransToday'],
            counter['maleSymptomaticTransToday'],
            counter['femaleSymptomaticTransToday'],
            ]
    #NOTE: bincount returns np.int64, not int!
    maleDistPartnerships = np.bincount(
        [male.n_partners for male in males],
        minlength=nOutGroups
        )
    assert sum(maleDistPartnerships) == len(males)
    femaleDistPartnerships = np.bincount(
        [female.n_partners for female in females],
        minlength=nOutGroups
        )
    assert sum(femaleDistPartnerships) == len(females)
    maleDistHIVp = np.bincount(
        [male.n_partners for male in males if male.is_infected],
        minlength=nOutGroups
        )
    femaleDistHIVp = np.bincount(
        [female.n_partners for female in females if female.is_infected],
        minlength=nOutGroups
        )
    #the following were provided as output arrays in the EHG code
    for dist in maleDistPartnerships, femaleDistPartnerships, maleDistHIVp, femaleDistHIVp:
        data.extend(dist[:nOutGroups])
        if (len(dist) > nOutGroups):
            logging.warn('discarding very high partnership counts')
    assert len(data) == 2 + 6 + 4*nOutGroups
    assert all(isinstance(item, numbers.Integral) for item in data), \
        "non-integer data!\n{}".format([(item,type(item)) for item in data]) #chkchkchkchk
    #above shd match EHG's output; below are additions
    nMprimary = sum(1 for male in males if male.has_primary())
    nFprimary = sum(1 for female in females if female.has_primary())
    data += [nMprimary, nFprimary]
    #finally, write the data to file
    fout.write('\n')
    data = ','.join(str(d) for d in data)
    fout.write(data)


def seedInfections(males, females, day, schedule, params):
    """Return diseases; seed initial infections.
    """
    assert not any(f.is_infected for f in females)
    assert not any(m.is_infected for m in males)
    logging.info('Seed infections.')
    Disease = params['Disease']
    nMales = len(males)
    nFemales = len(females)
    #note: using pHIVseed instead of pctHIVseed
    _pmSeedHIV, _pfSeedHIV = params['pSeedHIV']
    nMaleSeed = int( round(nMales * _pmSeedHIV) )
    nFemaleSeed = int( round(nFemales * _pfSeedHIV) )
    if not(nFemaleSeed < len(females) and nMaleSeed < len(males)):
        raise ValueError('choose smaller seeds')
    #why bother to randomize?? (but harmless, & matches EHG) chk
    _prng = params['prng']
    idxF = _prng.permutation(nFemales)[:nFemaleSeed]  #chkchk just use sample
    idxM = _prng.permutation(nMales)[:nMaleSeed]
    diseases = list()
    for idx in idxF:
        diseases.append( females[idx].HIVSeedInfect(Disease, day, schedule) )
    for idx in idxM:
        diseases.append( males[idx].HIVSeedInfect(Disease, day, schedule) )
    assert nFemaleSeed == sum(f.is_infected for f in females)
    assert nMaleSeed == sum(m.is_infected for m in males)
    return diseases

def tally_transmission(day, transmitter, counter):
    """Return None; tally the transmission.
    :note: EHG increment these in partnership.py; we do this instead
    """
    assert isinstance(day, int)
    assert isinstance(transmitter, Person)
    assert isinstance(counter, defaultdict)
    # get the infection stage of the transmitting partner
    # and tally the stage for output
    stage = transmitter.get_HIV_stage(day);
    if transmitter.sex == 'M':
        if(stage == 'primary'):
            counter['malePrimaryTransToday'] += 1
        elif(stage == 'asymptomatic'):
            counter['maleAsymptomaticTransToday'] += 1
        elif(stage == 'sympotmatic'):
            counter['maleSymptomaticTransToday'] += 1
    elif transmitter.sex == 'F':
        if(stage == 'primary'):
            counter['femalePrimaryTransToday'] += 1
        elif(stage == 'asymptomatic'):
            counter['femaleAsymptomaticTransToday'] += 1
        elif(stage == 'sympotmatic'):
            counter['femaleSymptomaticTransToday'] += 1
    else:
        raise ValueError('Unknown sex: ' + str(transmitter.sex))

