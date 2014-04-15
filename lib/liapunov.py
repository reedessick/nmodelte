usage=""" a module containing useful methods for nmode_l.py, meant to compare two sets of integration data """

import numpy as np
from pygsl import fft

import nmode_state as nms

####################################################################################################
#
#  generate residuals
#
####################################################################################################
def compute_residuals(t1, t2, q1, q2, common_mode):
  """
  computes the residuals between the two integration data sets.

  common_modes must have the form [ ((n,l,m,s), modeNo1, modeNo2), ... ], such as provided by common_nlms

  Meant to be used with liapunov.py
  """
  start = max(t1[0], t2[0])
  stop = min(t1[-1], t2[-1])

  ### check time-steps for the networks
  dt1 = t1[1]-t1[0] # assumes equal spacing in time
  dt2 = t2[1]-t2[0]

  if dt1 < dt2:
    return compute_residuals(start, stop, t2, t1, q2, q1, [(m,j,i) for m,i,j in common_modes])

  residual_times = []
  residuals = [ [] for m in common_modes ]
  rq1 = [ [] for m in common_modes ]
  rq2 = [ [] for m in common_modes ]

  # set up integration data
  ind1 = 0
  ind2 = 0
  while t1[ind1] < start:
    ind1 += 1
  while t2[ind2] < start:
    ind2 += 1

  _t1 = t1[ind1]
  _t2 = t2[ind2]

  while _t1 <= stop:
    while _t2 < _t1: # ensure _t2 is bigger than _t1
      ind2 += 1
      _t2 = t2[ind2]

    residual_times.append( _t1 )

    if _t1 == _t2: # no interpolation needed
      for m, (nlms,i,j) in enumerate(common_modes):
        _q1 = np.array(q1[i][ind1])
        _q2 = np.array(q2[j][ind2])

        rq1[m].append( _q1 )
        rq2[m].append( _q2 )
        residuals[m].append( _q1 - _q2 )

    else: # linearly interpolate q2 to find values at this time
      for m, (nlms,i,j) in enumerate(common_modes):
        _q1 = np.array(q1[i][ind1])

        q2old = np.array(q2[j][ind2-1])
        q2new = np.array(q2[j][ind2])
        _q2 = q2old + (q2new-q2old)*(_t1-t2[ind2-1])/dt2 # linearly interpolate

        rq1[m].append( _q1 )
        rq2[m].append( _q2 )
        residuals[m].append( _q1 - _q2 )

    ind1 += 1
    _t1 = t1[ind1]

  return residual_times, residuals, rq1, rq2


####################################################################################################
#
#  compute statistics
#
####################################################################################################
def inner_product(x, y):
  """
  computes the inner product of two series. This is defined as 

  (x|y) = \sum_i { conjugate(x[i]) * y[i] }

  we assume x, y are single-mode nmode amplitude lists [ [R,I], [R,I], ...]
  """
  if len(x) != len(y):
    sys.exit("lengths do not agree in liapunov.inner_product")
  x = np.array([complex(r,-i) for r,i in x])
  y = np.array([complex(r,i) for r,i in y])

  return sum(x*y)

##################################################
def nmode_inner_product(q1, q2):
  """
  computes the inner product of every mode in q1, q2 (assumes modes are identical). returns a list of inner products as well as a total inner product.
  These are normalized by the self-inner product of each mode

  <x|y> = (x|y) * [(x|x)*(y|y)]**-0.5
  """
  N_m = len(q1)
  if N_m != len(q2):
    sys.exit("mode lists must have same length in liapunov.nmode_inner_product")

  norms1 = [inner_product(x,x)**0.5 for x in q1]
  norms2 = [inner_product(x,x)**0.5 for x in q2]

  inner_prods = []
  tot = 0.0
  tot_norm1 = 0.0
  tot_norm2 = 0.0
  for m in range(N_m):
    x = q1[m]
    y = q2[m]

    norm1 = inner_product(x,x) 
    norm2 = inner_product(y,y)
    cross = inner_product(x,y)

    inner_prods.append( cross/(norm1*norm2)**0.5 )
    tot += cross
    tot_norm1 += norm1
    tot_norm2 += norm2

  return tot/(tot_norm1*tot_norm2)**0.5, inner_prods

##################################################
def chi2(x,y):
  """ 
  computes a chi2 goodness of fit measure defined as follows
  
  chi2(x,y) = \sum_i { (x-y)**2 / ((x+y)*0.5) } / len(x)

  which is the difference in the two time series divided by their average, normalized by their length
  """
  if len(x) != len(y):
    sys.exit("lengths do not agree in liapunov.chi2")
  x = np.array([complex(r,i) for r,i in x])
  y = np.array([complex(r,i) for r,i in y])

  return 2*sum( (x-y)**2/(x+y) )/len(x)

##################################################
def nmode_chi2(x,y):
  """
  computes the chi2 of every mode in q1, q2 (assumes modes are identical). returns a list of chi2 statistics as well as a total chi2 statistic.
  """
  N_m = len(q1)
  if N_m != len(q2):
    sys.exit("mode lists must have same length in liapunov.nmode_inner_product")

  chi2s = []
  for m in range(N_m):
    chi2s.append( chi2(q1[m], q2[m]) )

  return sum(chi2s), chi2s

##################################################
def cross_correlation(x, y):
  """
  computes cross correlation between the two signals in the frequency domain. Returns a time series and frequency series
  """
  N = len(x)
  if N != len(y):
    sys.exit("lengths don't agree in liapunov.cross_correlation")
  cfx = fft.complex_forward_float([complex(r,-i) for r,i in x]) # fourier transform
  cfy = fft.complex_forward_float([complex(r,i) for r,i in y])

  fcc = cfy * cfx.conjugate()
  tcc = fft.complex_inverse_float( fcc ) # multiply and take inverse fourier transform

  return [[c.real, c.imag] for c in tcc], [[c.real, c.imag] for c in fcc] # return as a list of reals

##################################################
def nmode_cross_correlation(q1, q2):
  """
  computes the cross correlation of every moe in q1, q2 (assumes modes are identical).

  NOTE: this does not require a time argument and it is assumed that q1, q2 have are sampled at the same times, and this will correspond to the return values
    this also does not return a frequency argument because of similar reasons
  """
  N_m = len(q1)
  if N_m != len(q2):
    sys.exit("mode lists must have same length in liapunov.nmode_inner_product")

  tccs = []
  fccs = []
  for m in range(N_m):
    tcc, fcc = cross_correlation(q1[m], q2[m])
    tccs.append( tcc )
    fccs.append( fcc )

  return tccs, fccs
 
##################################################
def nmode_liapunov_exp(residual_times, residuals):
  """
  attempts to compute a liapunov exponent for the residuals (fit an exponential). We define the distance between two trajectories as the absolute value of the residual

  try to fit |r| = A*exp(y*t)

  computations are delegated to nmode_state.A_exp_fit
  """

  # compute network total
  TA_fit, Ty_fit, Tchi2_red = nms.A_exp_fit(residual_times, [ sum([nmu.amp(residuals[m][ind][0], residuals[m][ind][1])**2 for m in range(len(residuals))])**0.5 for ind in range(len(residuals[0])) ] )

  # compute for individual modes
  A_fit, y_fit, chi2_red = nms.A_exp_fit(residual_times, residuals)

  return (TA_fit[0], Ty_fit[0], Tchi2_red[0]), (A_fit, y_fit, chi2_red)

##################################################
def sum_sqr_err(residuals):
  """
  computes the sum-square-errors from residuals
  """
  SSE = []
  tot_SSE = 0.0
  Npts = False

  for r in residuals:
    if not Npts:
      Npts = len(r)
    elif Npts != len(r):
      sys.exit("inconsitent lengths in liapunov.sum_sqr_err")

    sse = sum( [nmu.amp(l1,l2)**2 for l1,l2 in r] )
    tot_SSE += sse
    SSE.append( sse/Npts )
    
  return tot_SSE/Npts, SSE


####################################################################################################
#   
#  plotting
#   
####################################################################################################
def amp_residuals_plot(residual_times, residuals, rq1, rq2, n_l_m=False, mode_nums=False):
  """
  a specialty projection which plots both sets of integration data as well as the amplitude of the residuals as a function of time.
  """ 

  # want to have amplitudes overlaid? This will be VERY messy with multiple modes
  #   an alternative is to have amplitudes on stacked axes on which I enforce equal limits/scaling
  #
  # regardless, we need to place residuals below and ensure the correct colors match with mode amplitudes (this should happend automatically?)

  print "WRITE ME"
  return False







