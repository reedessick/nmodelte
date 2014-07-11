#!python_alias
usage="""an executable to compute and report summary/state information about a network using integration data. "s" is for "summary" """

import sys
import nmode_state as nm_s
import nmode_utils as nm_u
import numpy as np
from optparse import *

####################################################################################################
#
#
#                        Parse input options
#
#
####################################################################################################

parser=OptionParser(usage=usage)

parser.add_option("-v", "--verbose", dest="verbose", action="store_true", help="print updates on the progress of the script")
parser.add_option("-s", "--stefilename", default=False, type="string", help="the file into which ste summary will be written. Default is sys.stdout")
parser.add_option("-p", "--pstefilename", default=False, type="string", help="the file into which ste summary information will be written. File is dumped into a pickle format")

parser.add_option("-l", "--logfilename", default=False, type="string", help="the location of the log-file for this data")
parser.add_option("", "--modes", default=False, type="string", help="list of modes to analyze.")

### time domain
parser.add_option("-F", "--filename", dest="filename", default=False, type="string", help="the name of the file you want to analyze")
parser.add_option("", "--tmin", default=False, type="float", help="minimum time in analysis.")
parser.add_option("", "--tmax", default=False, type="float", help="maximum time in analysis.")
parser.add_option("", "--tdownsample", default=False, type="int", help="consider only 1 point out of this many in the time domain")
parser.add_option("", "--tcurrent", default=False, type="string", help="the variable type for the time-domain data")

parser.add_option("", "--exponential-fit", action="store_true", help="return exponential fits for mode amplitudes")

parser.add_option("", "--amplitude", action="store_true", help="return the most recent amplitudes of each mode")
parser.add_option("", "--phase-lag", action="store_true", help="compute the phase lag. Returns these last value as well as average and standard deviation.")
parser.add_option("", "--yE", action="store_true", help="compute \dot{E} = 2*\gamma_i*|x_i|^2, summation implied.")

parser.add_option("", "--Hns", action="store_true", help="compute Hns for this system")
parser.add_option("", "--ddt_P-E", action="store_true", help="compute dE/dt*Porb for each mode in the system")
parser.add_option("", "--HintHns", action="store_true", help="compute Hint+Hns")
parser.add_option("", "--ddt_P-HintHns", action="store_true", help="compute d(Hint+Hns)/dt*Porb")

#parser.add_option("", "--adot-lin", action="store_true", help="compute \dot{a}, the time rate of change of the semimajor axis of the orbit using only LINEAR contributions to the interaction energy.")
#parser.add_option("", "--adot-lin-mean", action="store_true", help="compute an average of \dot{a} using all time data.")
#
#parser.add_option("", "--edot-lin", action="store_true", help="compute \dot{e}, the time rate of change of the eccentricity of the orbit using only LINEAR contributions to the inteaction energy.")
#parser.add_option("", "--edot-lin-mean", action="store_true", help="compute an average of \dot{e} using all time data.")
#
#parser.add_option("", "--Hint-dot-mean", action="store_true", help="compute an average of \dot{Hint} using all time data. H_int is the interaction Hamiltonian of the system, and is how the modal network talks to the orbit.")

### concentration fitting
parser.add_option("", "--conc-fit", default=False, action="store_true")
parser.add_option("", "--conc-fit-start", default=0, type="int")
parser.add_option("", "--conc-fit-stop", default=np.infty, type="int")
parser.add_option("", "--conc-fit-iters", default=100, type="int")
parser.add_option("", "--conc-fit-rtol", default=1e-8, type="float")

parser.add_option("", "--conc-verbose", default=False, action="store_true")

### freq domain
parser.add_option("", "--freqfilename", default=False, type="string", help="load frequency data from this file.")
parser.add_option("", "--fmin", default=False, type="float", help="the lower limit on the frequencies analyzed.")
parser.add_option("", "--fmax", default=False, type="float", help="the upper limit on the frequencies analyzed.")
parser.add_option("", "--fdownsample", default=False, type="int", help="consider only 1 point out of this many in the frequency domain")
parser.add_option("", "--fcurrent", default=False, type="string", help="the variable type for the frequency data")

parser.add_option("", "--freq-peaks", action="store_true", help="identify main frequency component by simple maximization.")
parser.add_option("", "--freq-peaks-fit", action="store_true", help="identify main frequency component by fitting lorentzians to the PSD of integration data.")
parser.add_option("", "--freq-peaks-verbose", action="store_true", help="print statements about progress of fits")
parser.add_option("", "--max-freqfit-iters", default=500, type="int", help="maximum number of iterations allowed when fitting PSD to lorentzian to identify frequencies.")
parser.add_option("", "--freqfit-rtol", default=1e-8, type="float", help="relative tolerance of fitting.")

### reporting options
parser.add_option("", "--dynamical-relevance", default=False, action="store_true")
parser.add_option("", "--analytics", default=False, action="store_true")

opts, args = parser.parse_args()

time_domain = ( opts.amplitude or opts.exponential_fit or opts.phase_lag or opts.yE or opts.Hns or opts.ddt_P_E or opts.HintHns or opts.ddt_P-HintHns )
freq_domain = ( opts.freq_peaks or opts.freq_peaks_fit )

if time_domain and (not opts.filename):
  opts.filename = raw_input("time-domain filename = ")
if time_domain and (not opts.tcurrent):
  opts.tcurrent = raw_input("variable type for time-domain = ")

if (freq_domain or opts.analytics) and (not opts.freqfilename):
  opts.freqfilename = raw_input("freq-domain filename = ")
if freq_domain and (not opts.fcurrent):
  opts.fcurrent = raw_input("variable type for freq-domain = ")

if (not opts.logfilename):
  opts.logfilename = raw_input("log filename = ")

####################################################################################################
#
#
#                     Load network and compute physical characteristics of the orbit
#
#
####################################################################################################
if opts.verbose: print "loading parameters from %s" % opts.logfilename
system = nm_u.load_log(opts.logfilename)
network = system.network
Mprim = system.Mprim
Mcomp = system.Mcomp
Rprim = system.Rprim
Oorb = system.Oorb

N_m = len(network) # set up a list of which modes to analyze
if opts.modes:
  modes = [int(l) for l in opts.modes.split()]
else:
  modes = range(N_m)

# set up termination condtion for amplitude/dissiaption conc fits
opts.conc_fit_stop = min(opts.conc_fit_stop, len(network))

# compute dimensionful parameters for the system
dimMprim = Mprim*nm_u.Msun
dimMcomp = Mcomp*nm_u.Mjup
dimRprim = Rprim*nm_u.Rsun

Porb = 2*np.pi/Oorb

# normalization energy for the eigenmodes
dimEo, these_units = nm_u.Eo(system)

# orbital semi-major axis
aorb, these_units = nm_u.aorb(system)

# orbital energy
Eorb, these_units = nm_u.Eorb(system)

####################################################################################################
#
#
#                                  time domain
#
#
####################################################################################################
if time_domain:
  if opts.verbose: print "loading time-domain data from %s" % opts.filename
  t_P, q, N_m = nm_u.load_out(opts.filename, tmin=opts.tmin, tmax=opts.tmax, downsample=opts.tdownsample)

  if not len(t_P):
    raise ValueError, "no data loaded from ",opts.filename

  if opts.amplitude:
    if opts.verbose: print "computing mode amplitudes"
    x = nm_s.compute_A(q, Eo=1.)
    mA = []
    sA = []
    lA = []
    for A in x:
      mA.append( nm_s.sample_mean(A) )
      sA.append( nm_s.sample_var(A, xo=mA[-1])**0.5 )
      lA.append( A[-1] )
    del x

    x = nm_s.compute_E(q, Eo=1.)
    mE = []
    sE = []
    lE = []
    totE = [sum([x[modeNo][i] for modeNo in range(N_m)]) for i in range(len(t_P))]
    mean_tot_E = nm_s.sample_mean(totE)
    var_tot_E = nm_s.sample_var(totE, xo=mean_tot_E)
    for E in x:
      mE.append( nm_s.sample_mean(E) )
      sE.append( nm_s.sample_var(E, xo=mA[-1])**0.5 )
      lE.append( E[-1] )
    del x, totE

    if opts.conc_fit:
      if opts.verbose: print "\tattempting to fit concentration diagram"
      mE_fp, mE_cvr, mE_rchi2 = nm_s.broken_PowLaw_fitter(range(opts.conc_fit_start+1, opts.conc_fit_stop+1), np.array(sorted(mE)[opts.conc_fit_start:opts.conc_fit_stop])/mean_tot_E, np.ones((opts.conc_fit_stop-opts.conc_fit_start)), max_iters=opts.conc_fit_iters, rtol=opts.conc_fit_rtol, verbose=opts.conc_verbose)

  if opts.exponential_fit:
    if opts.verbose: print "computing exponential fits for mode amplitudes"
    A_fit, y_fit, chi2_red = nm_s.A_exp_fit(t_P, q, P=Porb, Eo=1.)

  if opts.phase_lag:
    if opts.verbose: print "computing phase lags"
    x, xx = nm_s.phase_lag(t_P, q, network, Eo=1.)
    mAsin_lag = []
    sAsin_lag = []
    lAsin_lag = []
    for A in x:
      mAsin_lag.append( nm_s.sample_mean(A) )
      sAsin_lag.append( nm_s.sample_var(A, xo=mAsin_lag[-1])**0.5 )
      lAsin_lag.append( A[-1] )
    del x
    msin_lag = []
    ssin_lag = []
    lsin_lag = []
    for A in xx:
      msin_lag.append( nm_s.sample_mean(A) )
      ssin_lag.append( nm_s.sample_var(A, xo=msin_lag[-1])**0.5 )
      lsin_lag.append( A[-1] )
    del xx
#    Asin_lag, sin_lag = nm_s.phase_lag(t_P, q, network, Eo=1.)

  if opts.yE:
    if opts.verbose: print "computing (viscous) energy disipation rates"
    x, xx = nm_s.viscous_disp(q, network, Eo=1.)
    mean_tot_yE = nm_s.sample_mean(x)
    var_tot_yE = nm_s.sample_var(x, xo=mean_tot_yE)
    del x
    myE = []
    syE = []
    lyE = []
    for A in xx:
      myE.append( nm_s.sample_mean(A) )
      syE.append( nm_s.sample_var(A, xo=myE[-1])**0.5 )
      lyE.append( A[-1] )
    del xx

    if opts.conc_fit:
      if opts.verbose: print "\tattempting to fit concentration diagram"
      myE_fp, myE_cvr, myE_rchi2 = nm_s.broken_PowLaw_fitter(range(opts.conc_fit_start+1, opts.conc_fit_stop+1), np.array(sorted(myE)[opts.conc_fit_start:opts.conc_fit_stop])/mean_tot_yE, np.ones((opts.conc_fit_stop-opts.conc_fit_start)), max_iters=opts.conc_fit_iters, rtol=opts.conc_fit_rtol, verbose=opts.conc_verbose)

  if opts.ddt_P_E:
    if opts.verbose: print "computing rate of change of energy in each mode"
    if opts.tcurrent == "q":
      x, xx = nm_s.compute_ddt_P_E_q(t_P, q, system, Eo=1.)
    elif opts.tcurrent == "x":
      x, xx = nm_s.compute_ddt_P_E_x(t_P, q, system, Eo=1.)
    mean_tot_dE = nm_s.sample_mean(x)
    var_tot_dE = nm_s.sample_var(x, xo=mean_tot_dE)
    del x
    mdE = []
    sdE = []
    ldE = []
    for dE in xx:
      mdE.append( nm_s.sample_mean(dE) )
      sdE.append( nm_s.sample_var(dE, xo=mdE[-1])**0.5 )
      ldE.append( dE[-1] )
    del xx

  if opts.Hns:
    if opts.verbose: print "computing energy stored in the stellar hamiltonian"
    if opts.tcurrent == "q":
      Hns = nm_s.Hns_q(q, system.network, Eo=1.)
    elif opts.tcurrent == "x":
      Hns = nm_s.Hns_x(t_P, q, system, Eo=1.)
    mean_Hns = nm_s.sample_mean(Hns)
    var_Hns = nm_s.sample_var(Hns, xo=mean_Hns)
    del Hns

  if opts.HintHns:
    if opts.verbose: print "computing Hint+Hns"
    if opts.tcurrent == "q":
      HintHns = nm_s.Hint_plus_Hns_no_NLT_q(t_P, q, system.network, Eo=1.0)
    elif opts.tcurrent == "x":
      HintHns = nm_s.Hint_plus_Hns_no_NLT_x(t_P, q, system, Eo=1.0)
    mean_HintHns = nm_s.sample_mean(HintHns)
    var_HintHns = nm_s.sample_var(HintHns, xo=mean_HintHns)
    del HintHns

  if opts.ddt_P_HintHns:
    if opts.verbose: print "computing rate of change of Hint+Hns"
    if opts.tcurrent == "q":
      dHintHns = nm_s.ddt_P_Hint_plus_Hns_no_NLT_q(t_P, q, system, Eo=1.0)
    elif opts.tcurrent == "x":
      dHintHns = nm_s.ddt_P_Hint_plus_Hns_no_NLT_x(t_P, q, system, Eo=1.0)
    mean_dHintHns = nm_s.sample_mean(dHintHns)
    var_dHintHns = nm_s.sample_var(dHintHns, xo=mean_dHintHns)
    del dHintHns

####################################################################################################
#
#
#                                frequency domain
#
#
####################################################################################################
if freq_domain:
  if opts.verbose: print "loading freq-domain data from %s" % opts.freqfilename
  freq, fq, N_m = nm_u.load_out(opts.freqfilename, tmin=opts.fmin, tmax=opts.fmax, downsample=opts.fdownsample)

  if not len(freq):
    raise ValueError, "no data loaded from ", opts.freqfilename

  if opts.freq_peaks:
    if opts.verbose: print "finding peak frequencies by closest point in the FFT"
    x, _ = nm_s.find_freq_peaks(freq, fq, delta=0, __sort=True)
    freq_maxs=[]
    for A in x:
      try:
        freq_maxs.append( [A[0]] )
      except IndexError:
        freq_maxs.append( [] )

  if opts.freq_peaks_fit:
    if opts.verbose: print "fitting lorentizians to the peak frequencies in the FFT"
    fit_params, fit_params_covar = nm_s.fit_peaks_lorentzian(freq, fq, max_iters=opts.max_freqfit_iters, delta=0, rtol=opts.freqfit_rtol, verbose=opts.freq_peaks_verbose)

####################################################################################################
#
#
#                                 analytics
#
#
####################################################################################################
if opts.analytics:
  if not opts.freq_peaks:
    if opts.verbose: print "finding peak freuqencies by closest point in the FFT"
    freq_maxs, freq_mins = nm_s.find_freq_peaks(freq, fq, delta=0, __sort=True)

  if opts.verbose: print "computing the most unstable pairings using observed frequencies"
  freqs = []
  for f in freq_maxs:
    if len(f) == 0:
      freqs.append( 0 )
    else:
      freqs.append( f[0][0]*Oorb )
  Ethrs = nm_s.Ethrs(freqs, network)

  if opts.verbose: print "computing expected amplitudes from best pairings"
  Amps = []
  for modeNo in range(len(network)):
    try:
      m,i,j,kmij, Ethr = Ethrs[modeNo][0]
      Amps.append( nm_s.threeMode_equilib((m,i,j,kmij), freqs[m], network, verbose=opts.verbose) ) # assumes "m" is the parent mode
    except IndexError:
      Amps.append( "none" )

####################################################################################################
#
#
#                                   report data!
#
#
####################################################################################################
if opts.stefilename:
  if opts.verbose: print "writing summary to "+opts.stefilename
  stefile = open(opts.stefilename, "w")
else:
  stefile = sys.stdout
if opts.pstefilename:
  ### instantiate data structure for pickling
  sdata = {} # system totals
  sdata["unit_system"] = nm_u.unit_system
  sdata["stats"] = {}

  mdata = {} # mode specifics

### print parameters of the data set
if time_domain:
  print >> stefile, """TIMEdomain
	outfilename : %s
	var         : %s
        t_P[0]      : %.2f
        dt_P        : %.2f
        t_P[-1]     : %.2f""" % (opts.filename, opts.tcurrent, t_P[0], t_P[1]-t_P[0], t_P[-1]) 
  ### add to pickle object
  if opts.pstefilename:
    sdata["time_domain"] = {"outfilename":opts.filename, "var":opts.tcurrent, "t_P[0]":t_P[0], "dt_P":t_P[1]-t_P[0], "t_P[-1]":t_P[-1]}

if freq_domain:
  print >> stefile, """FREQdomain
	frqfilename : %s
        var         : %s
        f_F[0]      : %f
        df_F        : %f
        f_F[-1]     : %f""" % (opts.freqfilename, opts.fcurrent, freq[0], freq[1]-freq[0], freq[-1])
  ### add to pickle object
  if opts.pstefilename:
    sdata["freq_domain"] = {"frqfilename":opts.freqfilename, "var":opts.fcurrent, "f_F[0]":freq[0], "df_F":freq[1]-freq[0], "f_F[-1]":freq[-1]}

### print parameters of the system
print >> stefile, """SYSTEM
	logfilename : %s
        Nmodes      : %d
        Ngens       : %d
        Ntriples    : %d
	Mprim/Msun  : %f
        Rprim/Rsun  : %f
        Mcomp/Mjup  : %f""" % (opts.logfilename, N_m, len(network.gens()[0]), len(network.to_triples()), Mprim, Rprim, Mcomp)
_val, _valexp = nm_u.float_to_scientific(dimEo)
print >> stefile, "\tEo          : %fe%d %s" % (_val, _valexp, nm_u.units['energy'])
_val, _valexp = nm_u.float_to_scientific(Porb)
print >> stefile, "\tPorb        : %fe%d %s" % (_val, _valexp, nm_u.units['time'])
_val, _valexp = nm_u.float_to_scientific(Oorb)
print >> stefile, "\tOorb        : %fe%d rad/%s" % (_val, _valexp, nm_u.units['time'])
print >> stefile, "\tecc         : %f%d" % nm_u.float_to_scientific(system.eccentricity)
_val, _valexp = nm_u.float_to_scientific(aorb)
print >> stefile, "\taorb        : %fe%d %s" % (_val, _valexp, nm_u.units['length'])
_val, _valexp = nm_u.float_to_scientific(Eorb)
print >> stefile, "\tEorb        : %fe%d %s" % (_val, _valexp, nm_u.units['energy'])
### add to pickle object
if opts.pstefilename:
  sdata["system"] = {"logfilename":opts.logfilename, "Nmodes":N_m, "Ngens":len(network.gens()[0]), "Ntriples":len(network.to_triples()), "Mprim/Msun":Mprim, "Rprim/Rsun":Rprim, "Mcomp/Mjup":Mcomp, "Eo":dimEo, "Porb":Porb, "ecc":system.eccentricity, "Eorb":Eorb}

### print network totals
print >> stefile, "STATS"
if opts.Hns:
  print >> stefile, "\tmean{Hns/|Eorb|}        :  %fe%d" % nm_u.float_to_scientific(mean_Hns*dimEo/abs(Eorb))
  print >> stefile, "\tstdv{Hns/|Eorb|}        :  %fe%d" % nm_u.float_to_scientific(var_Hns**0.5*dimEo/abs(Eorb))
  ### add to pickle object
  if opts.pstefilename:
    sdata["stats"].update( {"mean{Hns/|Eorb|}":mean_Hns*dimEo/abs(Eorb), "stdv{Hns/|Eorb|}":var_Hns**0.5*dimEo/abs(Eorb)} )

if opts.HintHns:
  print >> stefile, "\tmean{Hint+Hns}/|Eorb|   :  %fe%d" % nm_u.float_to_scientific(mean_HintHns*dimEo/abs(Eorb))
  print >> stefile, "\tstdv{Hint+Hns}/|Eorb|   :  %fe%d" % nm_u.float_to_scientific(var_HintHns**0.5*dimEo/abs(Eorb))
  ### add to pickle object
  if opts.pstefilename:
    sdata["stats"].update( {"mean{Hint+Hns}/|Eorb|":mean_HintHns*dimEo/abs(Eorb), "stdv{Hint+Hns}/|Eorb|":var_HintHns**0.5*dimEo/abs(Eorb)} )

if opts.ddt_P_HintHns:
  print >> stefile, "\tmean{d(Hint+Hns)/dt}*Porb/|Eorb| : %fe%d" % nm_u.float_to_scientific(mean_dHintHns*dimEo/abs(Eorb))
  print >> stefile, "\tstdv{d(Hint+Hns)/dt}*Porb/|Eorb| : %fe%d" % nm_u.float_to_scientific(var_dHintHns**0.5*dimEo/abs(Eorb))
  ### add to pickle object
  if opts.pstefilename:
    sdata["stats"].update( {"mean{d(Hint+Hns)/dt}*Porb/|Eorb|": mean_dHintHns*dimEo/abs(Eorb), "stdv{d(Hint+Hns)/dt}*Porb/|Eorb|":var_dHintHns**0.5*dimEo/abs(Eorb)} )

if opts.amplitude:
  print >> stefile, "\tmean{sum{E}/|Eorb|}      : %fe%d" % nm_u.float_to_scientific(mean_tot_E*dimEo/abs(Eorb))
  print >> stefile, "\tstdv{sum{E}/|Eorb|}      : %fe%d" % nm_u.float_to_scientific(var_tot_E**0.5*dimEo/abs(Eorb))
  mE_gini = nm_s.gini_index([0.0]+sorted(mE))
  print >> stefile, "\tGini_index{mean{sum{E}}} : %.5f" % mE_gini
  ### add to pickle object
  if opts.pstefilename:
    sdata["stats"].update( {"mean{sum{E}/|Eorb|}":mean_tot_E*dimEo/abs(Eorb), "stdv{sum{E}/|Eorb|}":var_tot_E**0.5*dimEo/abs(Eorb), "Gini_index{mean{sum{E}}}":mE_gini} )

if opts.yE:
  print >> stefile, "\tmean{sum{Edot}*(Porb/|Eorb|)} : %fe%d" %  nm_u.float_to_scientific(mean_tot_yE*dimEo * (Porb/abs(Eorb)))
  print >> stefile, "\tstdv{sum{Edot}*(Porb/|Eorb|)} :  %fe%d" % nm_u.float_to_scientific(var_tot_yE**0.5*dimEo * (Porb/abs(Eorb)))
  myE_gini = nm_s.gini_index([0]+sorted(list(-np.array(myE))))
  print >> stefile, "\tGini_index{mean{sum{Edot}}}   :  %.5f" % myE_gini
  ### add to pickle object
  if opts.pstefilename:
    sdata["stats"].update( {"mean{|sum{Edot}|*(Porb/|Eorb|)}":abs(mean_tot_yE*dimEo*(Porb/Eorb)), "stdv{|sum{Edot}|*(Porb/|Eorb|)}":abs(var_tot_yE**0.5*dimEo*(Porb/Eorb)), "Gini_index{mean{sum{Edot}}}":myE_gini} )

if opts.ddt_P_E:
  print >> stefile, "\tmean{sum{dE/dt*Porb}/|Eorb|}  : %fe%d" % nm_u.float_to_scientific(mean_tot_dE*dimEo/abs(Eorb))
  print >> stefile, "\tstdv{sum{dE/dt*Porb}/|Eorb|}  : %fe%d" % nm_u.float_to_scientific(var_tot_dE**0.5*dimEo/abs(Eorb))
#  mdE_gini = nm_s.gini_index([0]+sorted(list(np.array(mdE)-min(mdE))))
#  print >> stefile, "\tGini_index{mean{dE/dt}}       : %.5f" % mdE_gini
  ### add to pickle object
  if opts.pstefilename:
#    sdata["stats"].update( {"mean{sum{dE/dt*Porb}/|Eorb|}": mean_tot_dE*dimEo/abs(Eorb), "stdv{sum{dE/dt*Porb}/|Eorb|}":var_tot_dE**0.5*dimEo/abs(Eorb), "Gini_index{mean{sum{dE/dt}}}": mdE_gini } )
    sdata["stats"].update( {"mean{sum{dE/dt*Porb}/|Eorb|}": mean_tot_dE*dimEo/abs(Eorb), "stdv{sum{dE/dt*Porb}/|Eorb|}":var_tot_dE**0.5*dimEo/abs(Eorb) } )

### fitted parameters
if opts.conc_fit:
  if opts.amplitude:
    print >> stefile, """FITS : mean{E_i}
	C        : %.6f +/- %.6f 
	alpha    : %.6f +/- %.6f
	beta     : %.6f +/- %.6f
	gamma    : %.6f +/- %.6f
	red_chi2 : %.6f""" % (mE_fp[0], mE_cvr[0]**0.5, mE_fp[1], mE_cvr[1]**0.5, mE_fp[2], mE_cvr[2]**0.5, mE_fp[3], mE_cvr[3]**0.5, mE_rchi2)
 
  if opts.yE:
    print >> stefile, """FITS : mean{-2\sum{y_i*E_i}}
        C        : %.6f +/- %.6f 
        alpha    : %.6f +/- %.6f
        beta     : %.6f +/- %.6f
        gamma    : %.6f +/- %.6f
        red_chi2 : %.6f""" % (myE_fp[0], myE_cvr[0]**0.5, myE_fp[1], myE_cvr[1]**0.5, myE_fp[2], myE_cvr[2]**0.5, myE_fp[3], myE_cvr[3]**0.5, myE_rchi2)

########## values for each mode
if opts.dynamical_relevance:
  if opts.verbose: print "sorting modes by mean{A_i**2}"
  if not opts.amplitude:
    mE = [ nm_s.mean(E) for E in nms.compute_E(q)]
  dyn_rel_map = [ (modeNo, mE[modeNo]) for modeNo in modes ]
  dyn_rel_map.sort(key=lambda l: l[1], reverse=True) # sort so biggest amplitudes are first
  modes = [ l[0] for l in dyn_rel_map ]

for modeNo in modes:
  print >> stefile, "mode %d" % modeNo
  w, y, U = network.modes[modeNo].get_wyU()
  n, l, m = network.modes[modeNo].get_nlm()
  print >> stefile, "\tn=%d, l=%d, m=%d" % (n, l, m)
  _val, _valexp = nm_u.float_to_scientific(w)
  print >> stefile, "\tw=%fe%d rad/%s" % (_val, _valexp, nm_u.units['time'])
  _val, _valexp = nm_u.float_to_scientific(y)
  print >> stefile, "\ty=%fe%d %s" % (_val, _valexp, nm_u.units['frequency'])

  ### mdata for pickle file
  if opts.pstefilename:
    mdata[modeNo] = {"n":n, "l":l, "m":m, "w":w, "y":y}

  ### analytics
  if opts.analytics:
    try:
      m,i,j, kmij, Ethr = Ethrs[modeNo][0]
      _val, _valexp = nm_u.float_to_scientific(Ethr**0.5)
      print >> stefile, "\tminimum 3mode Ethr from pairing: (%d, %d, %d)\n\t\t\tAthr = %fe%d\n\t\t\tk_%d,%d,%d = %f" % (m,i,j, _val, _valexp, m,i,j, kmij)
    except IndexError:
      print >> stefile, "\t\t--this mode is not coupled--"

    if Amps[modeNo] != "none":
      (m,i,j), (Am, Ai, Aj), _, _ = Amps[modeNo] # we neglect phase information
      _valm, _valmexp = nm_u.float_to_scientific(Am)
      _vali, _valiexp = nm_u.float_to_scientific(Ai)
      _valj, _valjexp = nm_u.float_to_scientific(Aj)
      print >> stefile, "\t\texpected equilbirium amplitudes are:\n\t\t\tA_%d = %fe%d\n\t\t\tA_%d = %fe%d\n\t\t\tA_%d = %fe%d" % (m,_valm, _valmexp, i, _vali, _valiexp, j, _valj, _valjexp)

  ### time domain
  if time_domain:
    print >> stefile, "\ttime-domain variable type: %s" % opts.tcurrent

    ### mdata
    if opts.pstefilename:
      mdata[modeNo].update( {"var":opts.tcurrent} )

  if opts.yE:
    _val, _valexp = nm_u.float_to_scientific(dimEo*lyE[modeNo])
    print >> stefile, "\t\t            -2*Eo*y_i*A_i**2 @ max{t/Porb} : %fe%d %s/%s" % (_val, _valexp, nm_u.units['energy'], nm_u.units['time'])
    _val, _valexp = nm_u.float_to_scientific(dimEo*myE[modeNo])
    print >> stefile, "\t\t                    mean{-2*Eo*y_i*A_i**2} : %fe%d %s/%s" % (_val, _valexp, nm_u.units['energy'], nm_u.units['time'])
    _val, _valexp = nm_u.float_to_scientific(dimEo*syE[modeNo])
    print >> stefile, "\t\t                    stdv{-2*Eo*y_i*A_i**2} : %fe%d %s/%s" % (_val, _valexp, nm_u.units['energy'], nm_u.units['time'])
    print >> stefile, "\t\t   fraction of mean{sum{-2*Eo*y_i*A_i**2}} : %fe%d" % nm_u.float_to_scientific(myE[modeNo]/mean_tot_yE)

    ### mdata
    if opts.pstefilename:
      mdata[modeNo].update( {"mean{Edot}":myE[modeNo], "stdv{Edot}":syE[modeNo], "frac_mean{sum{Edot}}":myE[modeNo]/mean_tot_yE} )

  if opts.amplitude:
    print >> stefile, "\t\t                         A_i @ max{t/Porb} : %fe%d" % nm_u.float_to_scientific(lA[modeNo])
    print >> stefile, "\t\t                                 mean{A_i} : %fe%d" % nm_u.float_to_scientific(mA[modeNo])
    print >> stefile, "\t\t                                 stdv{A_i} : %fe%d" % nm_u.float_to_scientific(sA[modeNo])
    print >> stefile, "\t\t                         E_i @ max{t/Porb} : %fe%d" % nm_u.float_to_scientific(lE[modeNo])
    print >> stefile, "\t\t                                 mean{E_i} : %fe%d" % nm_u.float_to_scientific(mE[modeNo])
    print >> stefile, "\t\t                                 stdv{E_i} : %fe%d" % nm_u.float_to_scientific(sE[modeNo])
    print >> stefile, "\t\t                fraction of mean{sum{E_i}} : %fe%d" % nm_u.float_to_scientific(mE[modeNo]/mean_tot_E)

    ### mdata
    if opts.pstefilename:
      mdata[modeNo].update( {"mean{A}":mA[modeNo], "stdv{A}":sA[modeNo], "mean{E}":mE[modeNo], "stdv{D}":sE[modeNo], "frac_mean{sum{E}}":mE[modeNo]/mean_tot_E} )

  if opts.ddt_P_E:
    _val, _valexp = nm_u.float_to_scientific(dimEo*ldE[modeNo]/system.Porb)
    print >> stefile, "\t\t                        dE_i/dt @ max{t_P} : %fe%d %s/%s" % (_val, _valexp, nm_u.units['energy'], nm_u.units['time'])
    _val, _valexp = nm_u.float_to_scientific(dimEo*mdE[modeNo]/system.Porb)
    print >> stefile, "\t\t                             mean{dE_i/dt} : %fe%d %s/%s" % (_val, _valexp, nm_u.units['energy'], nm_u.units['time'])
    _val, _valexp = nm_u.float_to_scientific(dimEo*sdE[modeNo]/system.Porb)
    print >> stefile, "\t\t                             stdv{dE_i/dt} : %fe%d %s/%s" % (_val, _valexp, nm_u.units['energy'], nm_u.units['time'])
    print >> stefile, "\t\t            fraction of mean{sum{dE/dt}} : %fe%d" % nm_u.float_to_scientific(mdE[modeNo]/mean_tot_dE)

    ### mdata
    if opts.pstefilename:
      mdata[modeNo].update( {"mean{dE_i/dt}":mdE[modeNo], "stdv{dE_i/dt}":sdE[modeNo], "frac_mean{sum{dE/dt}}":mdE[modeNo]/mean_tot_dE} )

  if opts.exponential_fit:
    print >> stefile, "\t\t                       A_fit @ max{t/Porb} : %fe%d" % nm_u.float_to_scientific(A_fit[modeNo])
    _val, _valexp = nm_u.float_to_scientific(y_fit[modeNo])
    print >> stefile, "\t\t                                     y_fit : %fe%d %s" % (_val, _valexp, nm_u.units['frequency'])
    print >> stefile, "\t\t                            reduced \chi^2 : %f " % chi2_red[modeNo]

  if opts.phase_lag:
    print >> stefile, "\t\t    sin(m_i*Oorb*t - \phi_i) @ max{t/Porb} : %f" % lsin_lag[modeNo]
    print >> stefile, "\t\t            mean{sin(m_i*Oorb*t - \phi_i)} : %f" % msin_lag[modeNo]
    print >> stefile, "\t\t            stdv{sin(m_i*Oorb*t - \phi_i)} : %f" % ssin_lag[modeNo]

    print >> stefile, "\t\tA_i*sin(m_i*Oorb*t - \phi_i) @ max{t/Porb} : %fe%d" % nm_u.float_to_scientific( lAsin_lag[modeNo] )
    print >> stefile, "\t\t        mean{A_i*sin(m_i*Oorb*t - \phi_i)} : %fe%d" % nm_u.float_to_scientific( mAsin_lag[modeNo] )
    print >> stefile, "\t\t        stdv{A_i*sin(m_i*Oorb*t - \phi_i)} : %fe%d" % nm_u.float_to_scientific( sAsin_lag[modeNo] )

    ### mdata
    if opts.pstefilename:
      mdata[modeNo].update( {"mean{sin_lag}":msin_lag[modeNo], "stdv{sin_lag}":ssin_lag[modeNo], "mean{A*sin_lag}":mAsin_lag[modeNo], "stdv{A*sin_lag}":sAsin_lag[modeNo]} )

  ### freq domain
  if freq_domain:
    print >> stefile, "\tfreq-domain variable type: %s" % opts.fcurrent

  if opts.freq_peaks:
    if len(freq_maxs[modeNo]) != 0:
      fo, PSD = freq_maxs[modeNo][0]
      print >> stefile, "\t\tmax PSD = %fe%d" % nm_u.float_to_scientific(PSD)
      print >> stefile, "\t\tf*Porb @ max PSD = %f" % fo
    else:
      print >> stefile, "\t\tno peak frequency detected!"

  if opts.freq_peaks_fit:
    _fit_params = fit_params[modeNo]
    _fit_params_covar = fit_params_covar[modeNo]
    if len(_fit_params) > 0:
      _fit_params_comb = zip(_fit_params, _fit_params_covar)
      _fit_params_comb.sort(key=lambda line: line[0][2], reverse=True) # sort so highest amplitude is first
      [ffit_fP, ffit_yP, ffit_A], [ffit_fP_var, ffit_yP_var, ffit_A_var] = _fit_params_comb[0]
      _val, _valexp = nm_u.float_to_scientific(ffit_A)
      _val1, _valexp1 = nm_u.float_to_scientific(ffit_A_var**0.5)
      print >> stefile, "\t\tmax fit A = %fe%d +/- %fe%d" % (_val, _valexp, _val1, _valexp1)
      _val, _valexp = nm_u.float_to_scientific(ffit_fP_var**0.5)
      print >> stefile, "\t\tfo*Porb for max fit A = %f +/- %fe%d" % (ffit_fP, _val, _valexp)
      _val, _valexp = nm_u.float_to_scientific(ffit_yP)
      _val1, _valexp1 = nm_u.float_to_scientific(ffit_yP_var**0.5)
      print >> stefile, "\t\ty*Porb for max fit A = %fe%d +/- %fe%d" % (_val, _valexp, _val1, _valexp1)
    else:
      print >> stefile, "\t\tno fitted frequency peaks detected!"

if opts.stefilename:
  stefile.close()

if opts.pstefilename:
  if opts.verbose: print "writing pickled summary info to ",opts.pstefilename
  import pickle
  pstefile = open(opts.pstefilename, "w")
  pickle.dump(sdata, pstefile)
  pickle.dump(mdata, pstefile)
  pstefile.close()
