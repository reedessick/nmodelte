usage="""computes relevant parameters for g-modes"""

import numpy as np
import nmode_utils as nm_u

import mode_selection as ms

import networks

####################################################################################################
#
#
#                           gmode class
#
#
####################################################################################################
class gmode(networks.mode):
  """
  a class representing g-modes.

  automatically computes w, y based of assymptotic relations for high-order gmodes in sun-like stars
  this is done through delegation to gmodes.compute_w() and gmodes.compute_y()
    alpha, c, wo have the meanings defined therein

  data includes:
    n
    l
    m
    w = alpha*l/n
    y = c * wo**3 * L**2 * wo**-1 | L**2 = l*(l+1)
    U  = [(Ui, hi), (Uj, hj), ...] where i,j are harmonics of the tide (used in forcing frequency for integration)
    alpha
    c
    wo
  """

  def __init__(self, n=False, l=False, m=False, alpha=False, c=False, wo=False, U=[]):
    self.mode_type = "gmode"
    self.n = n
    self.l = l
    self.m = m
    if (n and l and alpha):
      self.w = compute_w(n,l,alpha)
    else:
      self.w = False
    if (n and l and c and wo and alpha):
      self.y = compute_y(n, l, c, wo, alpha)
    else:
      self.y = False
    self.U = U
    self.alpha = alpha
    self.c = c
    self.wo = wo
    self.check()

  ###
  def update(self):
    self.w = compute_w(self.n, self.l, self.alpha)
    self.y = compute_y(self.n, self.l, self.c, self.wo, self.alpha)    

  ### for reading from pickle 
  def from_tuple(self, nlmwyU, wo=1e-5):
    self.n, self.l, self.m, self.w, self.y, self.U = nlmwyU
    self.wo = wo
    self.alpha = compute_alpha(self.w, self.n, self.l)
    self.c = compute_c(self.y, wo, self.w, self.l)
    return self

  ### for reading from ASCii
  def from_str_nlmwy(self, string, wo=1e-5):
    n, l, m, w, y = string.strip().split()
    self.n = int(n)
    self.l = int(l)
    self.m = int(m)
    self.w = float(w)
    self.y = float(y)
    self.wo = wo
    self.alpha = compute_alpha(self.w, self.n, self.l)
    self.c = compute_c(self.y, wo, self.w, self.l)
    return self

  ###
  def compute_forcing(self, Porb, eccentricity, Mprim, Mcomp, Rprim, Ialm_hat=2.5e-3):
    """
    Assumes a simple analytic model for the linear forcing coefficient, taken from Weinberg's paper.

    This should be applicable to all eccentricities (proper scaling, etc).
    Without any hansen coefficient, we have:

      Ua(t) = [ Mcomp/(Mcomp+Mprim) * W_lm * I_alm * (Rprim/aorb)**(l+1) ] * [ (aorb/r(t))**(l+1) * exp(-i*m*phi(t)) ]

    The hansen coefficients expand the time-dependent factor into a fourier series in harmonics of the orbital period. The scaling is correct for the implemented solutions, but the prefactor may be slightly off.

    NOTE: if eccentricity is NOT 0.0, then we set the harmonic number to 'all' which will break most of your integrators, except for dxdphi_no_NLT
    """

    if self.l != 2:
      self.U = [] # no forcing for l!=2 modes with eccentricity==0
    else:
      if eccentricity == 0:
        self.U = [ (compute_U(self, Mprim, Mcomp, Rprim, Porb, Ialm_hat=Ialm_hat), 1) ]
      else:
        self.U = [ (compute_U(self, Mprim, Mcomp, Rprim, Porb, Ialm_hat=Ialm_hat), "all") ]

    return self

  ###
  def global_travel_time(self):
    """
    computes the global travel time of the mode
    """
    return 2*np.pi*self.n**2/(self.alpha*self.l)

  ###
  def is_local(self):
    """
    checks to see if the global travel time is longer than the damping rate. If it is, the wave is local!
      t*y >= 1 ==> local wave!
    """
    return self.global_travel_time()*self.y >= 1.0

  ###
  def is_global(self):
    """ not is_local() """
    return (not self.is_local())

  ###
  def pessimistic_breaking_thr(self, thr=1.0):
    """
    computes the amplitude the mode must reach before it begins breaking.
      k_r * xi_r >= thr
    assumes maximum is reached near inner turning point
    """
    prefact = 4729.502999063644 ### numerically evaluated prefactor using constants from Weinberg (2012)
    kr_xir = prefact * abs(self.wo / self.w )**3 * (self.l*(self.l+1))**0.5 # times mode amplitude
    return thr / kr_xir ### this should be the mode amplitude at which the mode will break.

  ###
  def breaking_thr(self, thr=1.0, rtol=1e-6, max_iters=100, verbose=False):
    """
    computes the amplitude the mode must reach before it begins breaking
      k_r * xi_r >= thr
    solves a transcendental equation for the radius at which k_r*xi_r is extremized.
      we use Newtons method and require the solution to agree to within rtol for convergence
    """
    prefact1 = 4729.502999063644 ### numerically evaluated prefactor using constants from Weinberg (2012)
    L2 = self.l*(self.l+1)
    prefact2 = (L2)**(1.5) * abs(self.wo/self.w)**3

    ### use Newton's method to find the zero
    ### we solve 0.5*(y+xo) = tan(y)
    ### where y = k_r*r - xo
    ###       xo= L2**0.5 - pi/4

    xo = L2**0.5 - np.pi/4

    ### bounds on the solution
    if self.l == 1: ### the first allowable solution occurs in a different well of tan(y)
      miny = 5*np.pi/4
      maxy = 3*np.pi/2
    else:
      miny = np.pi/4
      maxy = np.pi/2

    ### initial guess
    oldy = miny + np.pi/8
    ###                           f(y)                                df(y)
    y = oldy - ( np.tan(oldy) - 0.5*(oldy+xo) )/ (np.cos(oldy)**-2 - 0.5)

    if y < miny:
      if verbose: print "WARNING: y<miny"
      y = miny
#      raise ValueError, "y < miny = %f"%miny
    if y > maxy:
      if verbose: print "WARNING: y>maxy"
      y = maxy
#      raise ValueError, "y > maxy = %f"%maxy

    for _ in xrange(max_iters):
      if abs(y-oldy) < rtol*oldy:
        break
      oldy = y
      y = oldy - ( np.tan(oldy) - 0.5*(oldy+xo) )/ (np.cos(oldy)**-2 - 0.5)

      if y < miny:
        if verbose: print "WARNING: y<miny"
        y = miny
#        raise ValueError, "y < miny = %f"%miny
      if y > maxy:
        if verbose: print "WARNING: y>maxy"
        y = maxy
#        raise ValueError, "y > maxy = %f"%maxy
    else:
      raise StandardError, "max_iters=%d exceeded without convergence to within rtol=%.9f"%(max_iters, rtol)

    k_rr = y + xo ### k_r*r at max{k_r*xi_r}

    kr_xir = prefact1 * prefact2 * k_rr**-2 ## times the mode amplitude
    return thr / kr_xir ### this should be the mode amplitude at which the mode will break


####################################################################################################
#
#
#                            utility functions
#
#
####################################################################################################
def compute_U(mode, Mprim, Mcomp, Rprim, Porb, Ialm_hat=2.5e-3):
  """
  computes linear tidal forcing amplitude for this mode (an instance of networks.mode())
   only includes the dependence on this particular mode. All other dependencies are subsumed into U_hat

  U = U_hat * ((2*pi/w) / (1dy) )**(-11/6)
  where
  U_hat = 10**-8 * (Mcomp/(Mcomp+Mprime)) * (Porb / 10dy)**-2

  should be valid for sun-like stars with l=2
  """
  n, l, m, w, y = mode.get_nlmwy()
  if l == 2:
    ### DEPRECATED: old result that assumes constant Ialm_hat
    ### maintained for sanity checking, but the new implementation appears to be consistent
#    Uhat = 1.0e-8 * (Mcomp*nm_u.Mjup / (Mcomp*nm_u.Mjup + Mprim*nm_u.Msun) ) * (Porb / 864000.)**-2 ### old result 
#    return Uhat*(2*np.pi/(abs(w)*86400.))**(-11./6)

    wo = ms.compute_wo(Mprim, Rprim)
    Ylm = (3*np.pi/10)**0.5 ### from a tabulated value of Ylm(theta=pi/2, phi=0) for l=2, m=+/-2
    return = Mcomp*nm_u.Mjup*(Rprim*nm_u.Rsun)**3 / (nm_u.G*Mprim*nm_u.Msun*(Mprim*nm_u.Msun + Mcomp*nm_u.Mjup)) * (2*np.pi/Porb)**2 * Ylm * Ialm_hat*(abs(w)/wo)**(11./6)

  else:
    sys.exit("no analytic expression exists for gmode linear forcing coefficients with l=$d" % l)

##################################################
def compute_w(n, l, alpha):
  """ w = alpha*l/n """
  return alpha*l/n

##################################################
def compute_y(n, l, c, wo, alpha):
  """ y = c * wo**3 * L**2 * w**-2 | L**2 = l*(l+1) """
  return c*wo**3/alpha**2 * n**2 * (1+1./l)

##################################################
def compute_alpha(w, n, l):
  """ alpha = w*n/l """
  return w*n/l

##################################################
def compute_c(y, wo, w, l):
  """ c = (y * w**2) / (L**2 * wo**3) | L**2 = l*(l+1) """
  return (y * w**2) / (l*(l+1) * wo**3)

##################################################
def compute_cwo3a2(y, n, l):
  return y / (n**2 * (1+1./l))
