usage="""written to add / remove modes from a network"""

import pickle, sys 
import nmode_utils as nm_u
import nmode_state as nm_s
import networks
import numpy

####################################################################################################
#
#
#                           identify subsets of modes
#
#
####################################################################################################
def downselect_Alast(q, Athr, network):
  """
  applies an amplitude cut on all modes listed in q using the last Ai in q.
  returns: list of (n,l,m) to be removed
  """
  # separate modes based on amplitude threshold

  remove = []
  for nlms, ind in network.modeNoD.items():
    A = sum([l**2 for l in q[ind][-1]])**0.5
    if A <= Athr:
      remove.append(nlms)
    else:
      pass

  return remove

##################################################
def downselect_Amean(q, Athr, network):
  """
  applies an amplitude cut on all modes listed in q using <Ai>.
  returns: list of (n,l,m) to be removed
  """

  Amean = [1.0*sum(Ai)/len(Ai) for Ai in nm_s.compute_A(q)] # average amplitudes
  remove = []
  for nlms, ind in network.modeNoD.items():
    if Amean[ind] <= Athr:
      remove.append(nlms)
    else:
      pass

  return remove

##################################################
def downselect_Afit(q, Athr, network, t_P=False):
  """
  applies an amplituded cut using the fit amplitude at t_P[-1] assuming an exponential decay
    also requires the fit decay constant to be less than zero (decay, not growth)
  returns: list of (n,l,m) to be removed

  if not t_P:
    we compute an evenly spaced time series with the length of A. This should only introduce a constant prefactor in the computation of 'yfit', but shouldn't change the decision about whether to keep a mode.

  """
  if not t_P:
    t_P = range(len(q[0]))

  Afit, yfit, chi2_red = nm_s.A_exp_fit(t_P, q)
  remove = []
  for nlms, ind in network.modeNoD.items():
    if (Afit[ind] <= Athr) and (yfit[ind] < 0):
      remove.append(nlms)

  return remove

##################################################
def within_bandwidth_analytic(min_w, max_w, system):
  """
  computes the expected frequencies assuming 3mode equilibrium solutions

  uses family-tree generation functions defined within networks.network object

  returns 

  """
  return [ (freq, system.network.modes[modeNo]) for freq, modeNo in system.compute_3mode_freqs() if ( (min_w <= freq ) and ( freq <= max_w ) ) ]

##################################################
def within_bandwidth_maxPSD(freq, fq, min_w, max_w, network, Oorb=False):
  """
  selects modes using the frequency data (freq, fq) and returns a list of modes

  modes = [(f, (n,l,m,w,y,U)), ...]

  if Oorb, we assume that freq contains (O/Oorb) data, and we convert it back to O [rad/sec] by multiplying by Oorb
  """
  maxs, mins = nm_s.find_freq_peaks(freq, fq, delta=0, __sort=True)

  modes = []
  for ind, fit in enumerate(maxs):
    if len(fit) == 0:
      sys.exit("No frequency peak detected for modeNo=%d" % ind)
    f,_ = fit[0]

    if Oorb:
      f = f*Oorb

    if (min_w <= f) and (f <= max_w):
      modes.append( (f, network.modes[ind]) )

  return modes

##################################################
def within_bandwidth_lorentzian(freq, fq, min_w, max_w, network, Oorb=False, rtol=1e-8, max_iters=100, bw=0, n_lowpass=False):
  """
  selects modes using the frequency data (freq, fq) and returns a list of modes

  modes = [(f, (n,l,m,w,y,U)), ...]

  if Oorb, we assume that freq contains (O/Oorb) data, and we convert it back to O [rad/sec] by multiplying by Oorb
  """
  fit_params, fit_params_covar = nm_s.fit_peaks_lorentzian(freq, fq, max_iters=50, verbose=True)

  modes = []
  for ind, fit in enumerate(fit_params):
    if len(fit) == 0:
      sys.exit("No frequency peak detected for modeNo=%d" % ind)
    f,_,_ = fit[0]

    if Oorb:
      f = f*Oorb

    if (min_w <= f) and (f <= max_w):
      modes.append( (f, network.modes[ind]) )

  return modes


