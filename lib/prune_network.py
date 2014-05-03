usage="""written to add / remove modes from a network"""

import nmode_state as nm_s
import mode_selection as ms
import networks

####################################################################################################
#
#
#                           identify subsets of modes
#
#
####################################################################################################
def downselect_Alast(q, Athr, network, mode_nums=None):
  """
  applies an amplitude cut on all modes listed in q using the last Ai in q.
  returns: list of (n,l,m) to be removed
  """
  # separate modes based on amplitude threshold
  if not mode_nums:
    mode_nums = range(len(network))

  remove = []
  for modeNo in mode_nums:
    A = sum([l**2 for l in q[modeNo][-1]])**0.5
    if A <= Athr:
      remove.append( network.modes[modeNo].get_nlms() )
    else:
      pass

  return remove

##################################################
def downselect_Amean(q, Athr, network, mode_nums=None):
  """
  applies an amplitude cut on all modes listed in q using <Ai>.
  returns: list of (n,l,m) to be removed
  """
  if not mode_nums:
    mode_nums = range(len(network))

  Amean = [1.0*sum(Ai)/len(Ai) for Ai in nm_s.compute_A(q)] # average amplitudes
  remove = []
  for modeNo in mode_nums:
    if Amean[modeNo] <= Athr:
      remove.append( network.modes[modeNo].get_nlms() )
    else:
      pass

  return remove

##################################################
def downselect_Afit(q, Athr, network, t_P=None, mode_nums=None):
  """
  applies an amplituded cut using the fit amplitude at t_P[-1] assuming an exponential decay
    also requires the fit decay constant to be less than zero (decay, not growth)
  returns: list of (n,l,m) to be removed

  if not t_P:
    we compute an evenly spaced time series with the length of A. This should only introduce a constant prefactor in the computation of 'yfit', but shouldn't change the decision about whether to keep a mode.

  """
  if not t_P:
    t_P = range(len(q[0]))

  if not mode_nums:
    mode_nums = range(len(network))

  Afit, yfit, chi2_red = nm_s.A_exp_fit(t_P, q)
  remove = []
  for modeNo in mode_nums:
    if (Afit[modeNo] <= Athr) and (yfit[modeNo] < 0):
      remove.append( network.modes[modeNo].get_nlms() )

  return remove

##################################################
def within_bandwidth_analytic(min_w, max_w, system, mode_nums=False):
  """
  computes the expected frequencies assuming 3mode equilibrium solutions

  uses family-tree generation functions defined within networks.network object
  """
  if not mode_nums:
    mode_nums = range(len(system.network))
  freq_map = dict( (modeNo, freq) for freq, modeNo in system.compute_3mdoe_freqs() )
  ans = []
  for modeNo in mode_nums:
    if (min_w <= freq) and (freq <= max_w):
      ans.append( (freq, system.network.modes[modeNo]) )
  return ans

#  return [ (freq, system.network.modes[modeNo]) for freq, modeNo in system.compute_3mode_freqs() if ( (min_w <= freq ) and ( freq <= max_w ) ) ]

##################################################
def within_bandwidth_maxPSD(freq, fq, min_w, max_w, network, Oorb=False, mode_nums=None):
  """
  selects modes using the frequency data (freq, fq) and returns a list of modes

  modes = [(f, (n,l,m,w,y,U)), ...]

  if Oorb, we assume that freq contains (O/Oorb) data, and we convert it back to O [rad/sec] by multiplying by Oorb
  """
  if not mode_nums:
    mode_nums = range(len(network))

  maxs, mins = nm_s.find_freq_peaks(freq, fq, delta=0, __sort=True)

  modes = []
#  for ind, fit in enumerate(maxs):
  for modeNo in mode_nums:
    fit = maxs[modeNo]
    if len(fit) == 0:
      raise StandardError, "No frequency peak detected for modeNo=%d" % modeNo
    f,_ = fit[0]

    if Oorb:
      f = f*Oorb

    if (min_w <= f) and (f <= max_w):
      modes.append( (f, network.modes[modeNo]) )

  return modes

##################################################
def within_bandwidth_lorentzian(freq, fq, min_w, max_w, network, Oorb=False, rtol=1e-8, max_iters=100, bw=0, n_lowpass=False, mode_nums=None):
  """
  selects modes using the frequency data (freq, fq) and returns a list of modes

  modes = [(f, (n,l,m,w,y,U)), ...]

  if Oorb, we assume that freq contains (O/Oorb) data, and we convert it back to O [rad/sec] by multiplying by Oorb
  """
  if not mode_nums:
    mode_nums = range(len(network))

  fit_params, fit_params_covar = nm_s.fit_peaks_lorentzian(freq, fq, max_iters=50, verbose=True)

  modes = []
#  for ind, fit in enumerate(fit_params):
  for modeNo in mode_nums:
    fit = fit_params[modeNo]
    if len(fit) == 0:
      raise StandardError, "No frequency peak detected for modeNo=%d" % modeNo
    f,_,_ = fit[0]

    if Oorb:
      f = f*Oorb

    if (min_w <= f) and (f <= max_w):
      modes.append( (f, network.modes[modeNo]) )

  return modes

##################################################
def num_couplings_greater_than(num_k, network, mode_nums=None):
  """ returns a set of modes with len(K) >= num_k """
  if not mode_nums:
    mode_nums = range(len(network))

  modes = []
  for modeNo in mode_nums:
    if len(network.K[modeNo]) >= num_k:
      modes.append( network.modes[modeNo] )

  return modes

##################################################
def num_couplings_less_than(num_k, network, mode_nums=None):
  """ returns a set of modes with len(K) <= num_k """
  if not mode_nums:
    mode_nums = range(len(network))

  modes = []
  for modeNo in mode_nums:
    if len(network.K[modeNo]) <= num_k:
      modes.append( network.modes[modeNo] )

  return modes

##################################################
def heuristic_greater_than(heuristic, system, freqs=None, mode_nums=None):
  """ returns a set of modes with heuristic greater than heuristic. If freqs is not supplied, we compute freqs with syste.compute_3mode_freqs()"""
  if not freqs:
    freqs = system.compute_3mode_freqs()
    freqs.sort(key=lambda l: l[0])
    freqs = [l[1] for l in freqs]
  if not mode_nums:
    mode_nums = range(len(system.network))

  # ms.compute_heuristic(O, w1, w2, y1, y2)
  raise StandardError, "write me"

##################################################
def heuristic_less_than(heuristic, system, freqs=None, mode_nums=None):
  """ returns a set of modes with heuristic less than heuristic. If freqs is not supplied, we compute freqs with syste.compute_3mode_freqs()"""
  if not freqs:
    freqs = system.compute_3mode_freqs()
    freqs.sort(key=lambda l: l[0])
    freqs = [l[1] for l in freqs]
  if not mode_nums:
    mode_nums = range(len(system.network))

  # ms.compute_heuristic(O, w1, w2, y1, y2)
  raise StandardError, "write me"

##################################################
def Ethr_greater_than(Ethr, system, freqs=None, mode_nums=None):
  """ returns a set of modes with Ethr greater than Ethr. If freqs is not supplied, we compute freqs with syste.compute_3mode_freqs()"""
  if not freqs:
    freqs = system.compute_3mode_freqs()
    freqs.sort(key=lambda l: l[0])
    freqs = [l[1] for l in freqs]
  if not mode_nums:
    mode_nums = range(len(system.network))

  #ms.compute_Ethr(O, w1, w2, y1, y2, k)
  raise StandardError, "write me"

##################################################
def Ethr_less_than(Ethr, system, freqs=None, mode_nums=None):
  """ returns a set of modes with Ethr less than Ethr. If freqs is not supplied, we compute freqs with syste.compute_3mode_freqs()"""
  if not freqs:
    freqs = system.compute_3mode_freqs()
    freqs.sort(key=lambda l: l[0])
    freqs = [l[1] for l in freqs]
  if not mode_nums:
    mode_nums = range(len(system.network))

  #ms.compute_Ethr(O, w1, w2, y1, y2, k)
  raise StandardError, "write me"

##################################################
def NmodeE_greater_than(Ethr, system, freqs=None, mode_nums=None):
  """ returns a set of modes with NmodeE greater than NmodeE. If freqs is not supplied, we compute freqs with syste.compute_3mode_freqs()"""
  if not freqs:
    freqs = system.compute_3mode_freqs()
    freqs.sort(key=lambda l: l[0])
    freqs = [l[1] for l in freqs]
  if not mode_nums:
    mode_nums = range(len(system.network))

  raise StandardError, "write me"

##################################################
def NmodeE_less_than(Ethr, system, freqs=None, mode_nums=None):
  """ returns a set of modes with NmodeE less than NmodeE. If freqs is not supplied, we compute freqs with syste.compute_3mode_freqs()"""
  if not freqs:
    freqs = system.compute_3mode_freqs()
    freqs.sort(key=lambda l: l[0])
    freqs = [l[1] for l in freqs]
  if not mode_nums:
    mode_nums = range(len(system.network))

  raise StandardError, "write me"


