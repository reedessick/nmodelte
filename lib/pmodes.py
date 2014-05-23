usage="""computes relevant parameters for p-modes"""

import sympy, sys
import numpy as np
import sympy.physics.wigner as sympy_phys_wigner
import nmode_utils as nm_u

import networks

#print "\nWARNING! pmodes.py is not yet complete\n"

####################################################################################################
#
#
#                           pmode class
#
#
####################################################################################################
class pmode(networks.mode):
  """
  WRITE ME

  a class representing p-modes.

  data includes:
    n
    l
    m
    w 
    y
    U  = [(Ui, hi), (Uj, hj), ...] where i,j are harmonics of the tide (used in forcing frequency for integration)
  """

  def __init__(self, n, l, m, w, y, U=[]):
    self.mode_type = "pmode"
    self.n = n
    self.l = l
    self.m = m
    self.w = w
    self.y = y
    self.U = U
    self.check()

####################################################################################################
#
#
#                            utility functions
#
#
####################################################################################################
def compute_U(mode, U_hat=1e-12):
  """
  WRITE ME
  """
  return 0

##################################################
def compute_Uhat(Mprim, Mcomp, Porb):
  """
  WRITE ME
  """
  return 0

##################################################
def compute_w(n, l, alpha):
  """ WRITE ME """
  return 0

##################################################
def compute_y(n, l, c, wo, alpha):
  """ WRITE ME """
  return 0

##################################################
def compute_alpha(w, n, l):
  """ WRITE ME """
  return 0

##################################################
def compute_c(y, wo, w, l):
  """ WRITE ME """
  return 0

##################################################
def compute_cwo3a2(y, n, l):
  """ WRITE ME """
  return 0


