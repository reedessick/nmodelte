usage=""" written to provide basic plotting functions that are common to many investigations """

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
pi = np.pi
infty = np.infty
from numpy import linalg
#from scipy import linalg

import nmode_utils as nm_u

####################################################################################################
#
#
#                            standard plotting functions
#                              labels should be set by caller
#
####################################################################################################
def amp_plot(x, y, n_l_m=False, mode_nums=False):
  fig = plt.figure()
  ax = plt.subplot(1,1,1)
  if not mode_nums:
    mode_nums = range(len(y) )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])
    ax.plot(x, [nm_u.amp(l1,l2) for l1,l2 in y[m]], label=label)

  return fig, ax

##################################################
def phs_plot(x, y, n_l_m=False, mode_nums=False):
  fig = plt.figure()
  ax = plt.subplot(1,1,1)
  if not mode_nums:
    mode_nums = range(len(y) )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])
    ax.plot(x, [nm_u.phs(l1,l2)/(2*pi) for l1,l2 in y[m]], label=label)

  return fig, ax

##################################################
def amp_phs_plot(x, y, n_l_m=False, mode_nums=False):
  fig = plt.figure()
  ax1 = plt.subplot(2,1,1)
  ax2 = plt.subplot(2,1,2)

  plt.subplots_adjust(hspace=0.05, wspace=0.05)

  if not mode_nums:
    mode_nums = range(len(y) )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])
    ax1.plot(x, [nm_u.amp(l1,l2) for l1,l2 in y[m]], label=label)
    ax2.plot(x, [nm_u.phs(l1,l2)/(2*pi) for l1,l2 in y[m]], label=label)

  return fig, ax1, ax2

##################################################
def real_plot(x, y, n_l_m=False, mode_nums=False):
  fig = plt.figure()
  ax = plt.subplot(1,1,1)
  if not mode_nums:
    mode_nums = range(len(y) )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])
    ax.plot(x, [l[0] for l in y[m]], label=label)

  return fig, ax

##################################################
def imag_plot(x, y, n_l_m=False, mode_nums=False):
  fig = plt.figure()
  ax = plt.subplot(1,1,1)
  if not mode_nums:
    mode_nums = range(len(y) )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])
    ax.plot(x, [l[1] for l in y[m]], label=label)

  return fig, ax

##################################################
def real_imag_plot(x, y, n_l_m=False, mode_nums=False):
  fig = plt.figure()
  ax1 = plt.subplot(2,1,1)
  ax2 = plt.subplot(2,1,2)

  plt.subplots_adjust(hspace=0.05, wspace=0.05)

  if not mode_nums:
    mode_nums = range(len(y) )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])
    ax1.plot(x, [l[0] for l in y[m]], label=label)
    ax2.plot(x, [l[1] for l in y[m]], label=label)

  return fig, ax1, ax2

##################################################
def multi_gen_amp_plot(x, y, gens, n_l_m=False, mode_nums=False):
  """
  returns a stacked plot with one panel for each generation
  expects gens to have the form: [ [modeNo00, modeNo01, modeNo02,...], [modeNo10, modeNo11, modeNo12,...], ...]
  """
  fig = plt.figure()
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  axs = []
  for genNo, gen in enumerate(gens):
    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    for m in gen:
      if mode_nums and (m not in mode_nums):
        continue
      label = "mode %d" % m
      if n_l_m:
        label += ":%s" % str(n_l_m[m])
      ax.plot(x, [nm_u.amp(l1,l2) for l1,l2 in y[m]], label=label)
    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs

##################################################
def multi_gen_phs_plot(x, y, gens, n_l_m=False, mode_nums=False):
  """
  returns a stacked plot with one panel for each generation
  expects gens to have the form: [ [modeNo00, modeNo01, modeNo02,...], [modeNo10, modeNo11, modeNo12,...], ...]
  """
  fig = plt.figure()
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  axs = []
  for genNo, gen in enumerate(gens):
    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    for m in gen:
      if mode_nums and (m not in mode_nums):
        continue
      label = "mode %d" % m
      if n_l_m:
        label += ":%s" % str(n_l_m[m])
      ax.plot(x, [nm_u.phs(l1,l2)/(2*pi) for l1,l2 in y[m]], label=label)
    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs

##################################################
def multi_gen_real_plot(x, y, gens, n_l_m=False, mode_nums=False):
  """
  returns a stacked plot with one panel for each generation
  expects gens to have the form: [ [modeNo00, modeNo01, modeNo02,...], [modeNo10, modeNo11, modeNo12,...], ...]
  """
  fig = plt.figure()
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  axs = []
  for genNo, gen in enumerate(gens):
    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    for m in gen:
      if mode_nums and (m not in mode_nums):
        continue
      label = "mode %d" % m
      if n_l_m:
        label += ":%s" % str(n_l_m[m])
      ax.plot(x, [l[0] for l in y[m]], label=label)
    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs

##################################################
def multi_gen_imag_plot(x, y, gens, n_l_m=False, mode_nums=False):
  """
  returns a stacked plot with one panel for each generation
  expects gens to have the form: [ [modeNo00, modeNo01, modeNo02,...], [modeNo10, modeNo11, modeNo12,...], ...]
  """
  fig = plt.figure()
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  axs = []
  for genNo, gen in enumerate(gens):
    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    for m in gen:
      if mode_nums and (m not in mode_nums):
        continue
      label = "mode %d" % m
      if n_l_m:
        label += ":%s" % str(n_l_m[m])
      ax.plot(x, [l[1] for l in y[m]], label=label)
    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  return fig, axs

####################################################################################################
#
# phase portraits
#
####################################################################################################
def _compute_time_derivatives(t,q,system,current):
  import network_flow as nf
  Porb = system.Porb

  N_m = len(q)
  N_p = len(q[0])
  dqdt_P = [ [] for n in range(N_m) ]

  if current == "x":
    for ind in range(N_p): # all points in data
      this_q = []
      for m in range(N_m):
        this_q += q[m][ind] # add real and imagninary part for each mode to single list
      this_dqdt_P = nf.dxdt_no_NLT(t[ind], this_q, system) # compute derivatives
      for m in range(N_m):
        dqdt_P[m].append( [this_dqdt_P[2*m]*Porb, this_dqdt_P[2*m+1]*Porb] ) # re-format list into familiar form

  elif current == "q":
    for ind in range(N_p): # all points in data
      this_q = []
      for m in range(N_m):
        this_q += q[m][ind] # add real and imagninary part for each mode to single list
      this_dqdt_P = nf.dqdt_no_NLT(t[ind], this_q, system) # compute derivatives
      for m in range(N_m):
        dqdt_P[m].append( [this_dqdt_P[2*m]*Porb, this_dqdt_P[2*m+1]*Porb] ) # re-format list into familiar form

  else:
    sys.exit("unknown variable type: %s in nmode_plotting.phase_portrait" % current)

  return dqdt_P, N_m, N_p

##################################################
def amp_phase_portrait(t, q, system, current, n_l_m=False, mode_nums=False, verbose=False):
  """
  computes and overlays phase portraits for each mode. These are phase portraits for the amplitude of the mode
  """
  if verbose: print "\tcomputing time derivatives"
  dqdt_P, N_m, N_p = _compute_time_derivatives(t,q,system,current)
  
  if verbose: print "\tcomputing A, \dot{A} and plotting"
  fig = plt.figure()
  ax = plt.subplot(1,1,1)
  
  if not mode_nums:
    mode_nums = range( N_m )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])

    this_q = q[m]
    this_dqdt_P = dqdt_P[m]

    ax.plot([nm_u.amp(l1,l2) for l1,l2 in this_q], [(this_q[ind][0]*this_dqdt_P[ind][0] + this_q[ind][1]*this_dqdt_P[ind][1]) / (nm_u.amp(this_q[ind][0], this_q[ind][1])) for ind in range(N_p)], label=label)

  return fig, ax

##################################################
def real_phase_portrait(t, q, system, current, n_l_m=False, mode_nums=False, verbose=False):
  """
  computes and overlays phase portraits for each mode. These are phase portraits for the real part of the mode
  """
  if verbose: print "\tcomputing time derivatives"
  dqdt_P, N_m, N_p = _compute_time_derivatives(t,q,system,current)
 
  if verbose: print "\tcomputing A, \dot{A} and plotting" 
  fig = plt.figure()
  ax = plt.subplot(1,1,1)
  
  if not mode_nums:
    mode_nums = range( N_m )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])

    ax.plot([l for l,_ in q[m]], [l for l,_ in dqdt_P[m]], label=label)

  return fig, ax

##################################################
def imag_phase_portrait(t, q, system, current, n_l_m=False, mode_nums=False, verbose=False):
  """
  computes and overlays phase portraits for each mode. These are phase portraits for the imaginary part of the mode
  """
  if verbose: print "\tcomputing time derivatives"
  dqdt_P, N_m, N_p = _compute_time_derivatives(t,q,system,current)

  if verbose: print "\tcomputing A, \dot{A} and plotting"
  fig = plt.figure()
  ax = plt.subplot(1,1,1)

  if not mode_nums:
    mode_nums = range( N_m )

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])

    ax.plot([l for _,l in q[m]], [l for _,l in dqdt_P[m]], label=label)

  return fig, ax

##################################################
def real_imag_phase_portrait(q, n_l_m=False, mode_nums=False):
  """
  plots the Real part of each mode agains the Imaginary part of the same mode
    In polar coordinates, this corresponds to amplitude and phase
  """
  fig = plt.figure()
  ax = plt.subplot(1,1,1)

  if not mode_nums:
    mode_nums = range(len(q))

  for m in mode_nums:
    label = "mode %d" % m
    if n_l_m:
      label += ":%s" % str(n_l_m[m])

    ax.plot([l for l,_ in q[m]], [l for _,l in q[m]], label=label)

  return fig, ax

##################################################
#
#  multi-gen phase portraits
#
##################################################
def multi_gen_amp_phase_portrait(t, q, system, current, gens, n_l_m=False, mode_nums=False, verbose=False):
  """
  computes and overlays phase portraits for each mode. These are phase portraits for the amplitude of the mode
  automatically separates modes by generation (specified by gens)
  """
  if verbose: print "\tcomputing time derivatives"
  dqdt_P, N_m, N_p = _compute_time_derivatives(t,q,system,current)

  if verbose: print "\tcomputing A, \dot{A} and plotting"
  fig = plt.figure()
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  g_min = infty
  g_max = -infty

  axs = []
  for genNo, gen in enumerate(gens):
    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    for m in gen:
      if mode_nums and (m not in mode_nums):
        continue
      label = "mode %d" % m
      if n_l_m:
        label += ":%s" % str(n_l_m[m])

      this_q = q[m]
      this_dqdt_P = dqdt_P[m]

      ax.plot([nm_u.amp(l1,l2) for l1,l2 in this_q], [(this_q[ind][0]*this_dqdt_P[ind][0] + this_q[ind][1]*this_dqdt_P[ind][1]) / (nm_u.amp(this_q[ind][0], this_q[ind][1])) for ind in range(N_p)], label=label)

    # figure out minimum/maximum
    this_min, this_max = ax.get_xlim()
    if this_min < g_min:
      g_min = this_min
    if this_max > g_max:
      g_max = this_max

    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  for ax in axs:
    ax.set_xlim(xmin=g_min, xmax=g_max)

  return fig, axs

##################################################
def multi_gen_real_imag_phase_portrait(q, gens, n_l_m=False, mode_nums=False):
  """
  computes and overlays phase portraits for each mode. These are phase portraits for the amplitude of the mode
  automatically separates modes by generation (specified by gens)
  """
  fig = plt.figure()
  num_gen = len(gens)
  ax_height = 0.8/num_gen
  buff = 0.01*ax_height

  g_min = infty
  g_max = -infty

  axs = []
  for genNo, gen in enumerate(gens):
    ax = fig.add_axes( [0.15, 0.95 - (1+genNo)*ax_height, 0.8, ax_height-buff] )
    for m in gen:
      if mode_nums and (m not in mode_nums):
        continue
      label = "mode %d" % m
      if n_l_m:
        label += ":%s" % str(n_l_m[m])
      ax.plot([l for l,_ in q[m]], [l for _,l in q[m]], label=label)

    # figure out minimum/maximum
    this_min, this_max = ax.get_xlim()
    if this_min < g_min:
      g_min = this_min
    if this_max > g_max:
      g_max = this_max

    plt.setp(ax.get_xticklabels(), visible=False)
    axs.append(ax)

  plt.setp(ax.get_xticklabels(), visible=True)

  for ax in axs:
    ax.set_xlim(xmin=g_min, xmax=g_max)

  return fig, axs

####################################################################################################
#
#                          Instantaneous linearization
#                  stability of perturbations about current flow
####################################################################################################
def eig_q(q, network, just_eigval=False):
  """
  computes the eigenvalues of system for the single sample in q. 
    q must have the form [R0, I0, R1, I1, R2, I2, R3, I3, ...]
  returns a list of eigenvalues and eigenvectors

  this is a linearization of the evolution of two infinitesimally close neighboring trajectories, one of which is currently at the point {q}
  """
  Nm = len(network)
  M = np.zeros((2*Nm, 2*Nm)) # holder for matrix we diagonalize

  ### build matrix
  for modeNo in range(Nm):
    wo, yo, _ = network.wyU[modeNo]

    ro_ind = 2*modeNo
    io_ind = ro_ind+1
    # diagonal 2x2 matrix    
    M[ro_ind][ro_ind] = -yo
    M[ro_ind][io_ind] =  wo
    M[io_ind][ro_ind] = -wo
    M[io_ind][io_ind] = -yo

    # off-diagonal elements
    for i,j,k in network.K[modeNo]:
      if i == j:
        ri_ind = 2*i
        ii_ind = ri_ind+1
        ri = q[ri_ind]
        ii = q[ii_ind]

        M[ro_ind][ri_ind] =  2*wo*k*ii
        M[ro_ind][ii_ind] =  2*wo*k*ri
        M[io_ind][ri_ind] =  2*wo*k*ri
        M[io_ind][ii_ind] = -2*wo*k*ii

      else: # i != j
        ri_ind = 2*i
        ii_ind = ri_ind+1
        ri = q[ri_ind]
        ii = q[ii_ind]

        rj_ind = 2*j
        ij_ind = rj_ind+1
        rj = q[rj_ind]
        ij = q[ij_ind]

        M[ro_ind][ri_ind] =  2*wo*k*ij
        M[ro_ind][ii_ind] =  2*wo*k*rj
        M[io_ind][ri_ind] =  2*wo*k*rj
        M[io_ind][ii_ind] = -2*wo*k*ij

        M[ro_ind][rj_ind] =  2*wo*k*ii
        M[ro_ind][ij_ind] =  2*wo*k*ri
        M[io_ind][rj_ind] =  2*wo*k*ri
        M[io_ind][ij_ind] = -2*wo*k*ij

  ### diagonalize matrix. get eigenvalues, vectors
  if just_eigvals:
    return linalg.eigvals(M, overwrite_a=True)

  eigenvalues, eigenvectors = linalg.eig(M, overwrite_a=True) # overwrite_a=True may help performance, and I don't care about "a" at this point (a=M within the function)

  ### now, the eigenvectors may be complex, and we really want to return the corresponding real and imaginary parts of each mode (which are mixed into two complex numbers). We convert this by hand
  converted_eigenvectors = np.zeros((2*Nm, 2*Nm))
  for vecNo in range(2*Nm):
    eigenvect = eigenvectors[:,vecNo]
    for modeNo in range(Nm):
      r = eigenvect[2*modeNo]
      i = eigenvect[2*modeNo+1]
      converted_eigenvectors[2*modeNo][vecNo]   = np.real(r) - np.imag(i)
      converted_eigenvectors[2*modeNo+1][vecNo] = np.imag(r) + np.real(i)

    ### correct normalization
    converted_eigenvectors[:,vecNo] = converted_eigenvectors[:,vecNo]/sum(converted_eigenvectors[:,vecNo]**2)**0.5

  return eigenvalues, converted_eigenvectors

##################################################
def eig_x(t_P, x, system, just_eigvals=False):
  """
  computes the eigenvalues of system for the single sample in q. 
    x must have the form [R0, I0, R1, I1, R2, I2, R3, I3, ...]
  returns a list of eigenvalues and eigenvectors

  this is a linearization of the evolution of two infinitesimally close neighboring trajectories, one of which is currently at the point {q}
  """
  network = system.network
  Nm = len(network)
  M = np.zeros((2*Nm, 2*Nm)) # holder for matrix we diagonalize

  t = t_P*system.Porb

  ### build matrix
  for modeNo in range(Nm):
    wo, yo, _ = network.wyU[modeNo]

    ro_ind = 2*modeNo
    io_ind = ro_ind+1
    # diagonal 2x2 matrix    
    M[ro_ind][ro_ind] = -yo
    M[ro_ind][io_ind] =  wo
    M[io_ind][ro_ind] = -wo
    M[io_ind][io_ind] = -yo

    # off-diagonal elements
    for i,j,k in network.K[modeNo]:
      if i == j:
        wi = network.modes[i].w

        phs = (wo + 2*wi)*t
        cosW = np.cos(phs)
        sinW = np.sin(phs)

        ri_ind = 2*i
        ii_ind = ri_ind+1
        ri = x[ri_ind]
        ii = x[ii_ind]

        M[ro_ind][ri_ind] = 2*wo*k*(cosW*ii - sinW*ri)
        M[ro_ind][ii_ind] = 2*wo*k*(cosW*ri + sinW*ii)
        M[io_ind][ri_ind] = 2*wo*k*(cosW*ri + sinW*ii)
        M[io_ind][ii_ind] = 2*wo*k*(cosW*ii + sinW*ri)

      else: # i != j
        wi = network.modes[i].w
        wj = network.modes[j].w

        phs = (wo+wi+wj)*t
        cosW = np.cos(phs)
        sinW = np.sin(phs)

        ri_ind = 2*i
        ii_ind = ri_ind+1
        ri = x[ri_ind]
        ii = x[ii_ind]

        rj_ind = 2*j
        ij_ind = rj_ind+1
        rj = x[rj_ind]
        ij = x[ij_ind]

        M[ro_ind][ri_ind] = 2*wo*k*(cosW*ij - sinW*rj)
        M[ro_ind][ii_ind] = 2*wo*k*(cosW*rj + sinW*ij)
        M[io_ind][ri_ind] = 2*wo*k*(cosW*rj + sinW*ij)
        M[io_ind][ii_ind] = 2*wo*k*(cosW*ij + sinW*rj)

        M[ro_ind][rj_ind] = 2*wo*k*(cosW*ii - sinW*ri)
        M[ro_ind][ij_ind] = 2*wo*k*(cosW*ri + sinW*ii)
        M[io_ind][rj_ind] = 2*wo*k*(cosW*ri + sinW*ii)
        M[io_ind][ij_ind] = 2*wo*k*(cosW*ii + sinW*ri)

  ### diagonalize matrix. get eigenvalues, vectors
  if just_eigvals:
    return linalg.eigvals(M, overwrite_a=True)

  eigenvalues, eigenvectors = linalg.eig(M, overwrite_a=True) # overwrite_a=True may help performance, and I don't care about "a" at this point (a=M within the function)

  ### now, the eigenvectors may be complex, and we really want to return the corresponding real and imaginary parts of each mode (which are mixed into two complex numbers). We convert this by hand
  converted_eigenvectors = np.zeros((2*Nm, 2*Nm))
  for vecNo in range(2*Nm):
    eigenvect = eigenvectors[:,vecNo]
    for modeNo in range(Nm):
      r = eigenvect[2*modeNo]
      i = eigenvect[2*modeNo+1]
      converted_eigenvectors[2*modeNo][vecNo]   = np.real(r) - np.imag(i)
      converted_eigenvectors[2*modeNo+1][vecNo] = np.imag(r) + np.real(i)

    ### correct normalization
    converted_eigenvectors[:,vecNo] = converted_eigenvectors[:,vecNo]/sum(converted_eigenvectors[:,vecNo]**2)**0.5

  return eigenvalues, converted_eigenvectors


##################################################
def nmode_eig_q(q, system, max_eig=False, verbose=False):
  """
  computes the eigenvalues of the sytem for each sample in q.
    assumes "q" is the integration variable
  Returns a list of eigenvalues and eigenvectors

  if max_eig: returns ONLY the eigenvalue/vector with max{Re{eigenvalue}} at each sample
  """
  if verbose: import time

  eigvals = []
  eigvecs = []
  
  L = len(q[0])
  for ind in range(L):
    if verbose: to = time.time()

    this_q = []
    for Q in q:
      this_q += Q[ind] # add real and imagninary part for each mode to single list
    eigenvalues, eigenvectors = eig_q(this_q, system.network)
    if max_eig:
      ev = [(np.real(e), ev_ind) for ev_ind, e in enumerate(eigenvalues)]
      ev.sort(key=lambda l: l[0], reverse=True)
      max_ev_ind = ev[0][1]
      eigenvalues = eigenvalues[max_ev_ind]
      eigenvectors = eigenvectors[:,max_ev_ind]

    eigvals.append(eigenvalues)
    eigvecs.append(eigenvectors)

    if verbose: print "\t%d / %d\t%f" % (ind+1, L, time.time()-to)

  return eigvals, eigvecs

##################################################
def nmode_eig_x(t_P, x, system, max_eig=False, verbose=False):
  """
  computes the eigenvalues of the sytem for each sample in x.
    assumes "x" is the integration variable
  Returns a list of eigenvalues and eigenvectors

  if max_eig: returns ONLY the eigenvalue/vector with max{Re{eigenvalue}} at each sample
  """
  if verbose: import time

  Porb = system.Porb
  network = system.network
  Nm = len(system.network)
  eigvals = []
  eigvecs = []

  L = len(t_P)
  for ind in range(L):

    if verbose: to = time.time()

    t = t_P[ind]*Porb
    this_x = []
    for modeNo in range(Nm):
      this_x += x[modeNo][ind]

    eigenvalues, eigenvectors = eig_x(t_P[ind], this_x, system)
    
    if max_eig:
      ev = [(np.real(e), ev) for ev, e in enumerate(eigenvalues)]
      ev.sort(key=lambda l: l[0], reverse=True)
      max_ev_ind = ev[0][1]
      eigenvalues = eigenvalues[max_ev_ind]
      eigenvectors = eigenvectors[:,max_ev_ind]

    eigvals.append(eigenvalues)
    eigvecs.append(eigenvectors)

    if verbose: print "\t%d / %d\t%f" % (ind+1, L, time.time()-to)

  return eigvals, eigvecs

##################################################
def nmode_eigval_q(q, system, max_eig=False, verbose=False):
  """
  computes the eigenvalues of the sytem for each sample in q.
    assumes "q" is the integration variable
  Returns a list of eigenvalues
  
  if max_eig: returns ONLY the eigenvalue/vector with max{Re{eigenvalue}} at each sample
  """
  if verbose: import time

  eigvals = []

  L = len(q[0])
  for ind in range(L):
    if verbose: to = time.time()

    this_q = []
    for Q in q: 
      this_q += Q[ind]
    eigenvalues = eig_q(this_q, system.network, just_eigvals=True)
    if max_eig:
      ev = [(np.real(e), ev_ind) for ev_ind, e in enumerate(eigenvalues)]
      ev.sort(key=lambda l: l[0], reverse=True)
      max_ev_ind = ev[0][1]
      eigenvalues = eigenvalues[max_ev_ind]
    eigvals.append(eigenvalues)

    if verbose: print "\t%d / %d\t%f" % (ind+1, L, time.time()-to)

  return eigvals

##################################################
def nmode_eigval_x(t_P, x, system, max_eig=False, verbose=False):
  """
  computes the eigenvalues of the sytem for each sample in x.
    assumes "x" is the integration variable
  Returns a list of eigenvalues

  if max_eig: returns ONLY the eigenvalue/vector with max{Re{eigenvalue}} at each sample
  """
  if verbose: import time

  eigvals = []

  L = len(t_P)
  for ind in range(L):
    if verbose: to = time.time()

    this_x = []
    for X in x:
      this_x += X[ind]

    eigenvalues = eig_x(t_P[ind], this_x, system, just_eigvals=True)

    if max_eig:
      ev = [(np.real(e), ev_ind) for ev_ind, e in enumerate(eigenvalues)]
      ev.sort(key=lambda l: l[0], reverse=True)
      max_ev_ind = ev[0][1]
      eigenvalues = eigenvalues[max_ev_ind]

    eigvals.append(eigenvalues)

    if verbose: print "\t%d / %d\t%f" % (ind+1, L, time.time()-to)

  return eigvals

##################################################
def eig_plot(t_P, vals, vecs1, vecs2, system, cmap=plt.cm.jet):
  """
  generates a plot showing the evolution of a single eigenvalue and corresponding eigenvector over time (given as inputs)
  plots the amplitude and phase of each mode
    vals must be an array: size (len(t_P))
    vecs1 must be an array: size (len(t_P), N_m)
    vecs2 must be an array: size (len(t_P), N_m)

  """
  Nm = len(vecs1[0])

  fig = plt.figure()

  ax1 = fig.add_axes([0.15, 0.75, 0.70, 0.2])
  ax1p = ax1.twinx()

  ax2 = fig.add_axes([0.15, 0.45, 0.70, 0.275])
  ax2_cb = fig.add_axes([0.855, 0.45, 0.020, 0.275])

  ax3 = fig.add_axes([0.15, 0.15, 0.70, 0.275])
  ax3_cb = fig.add_axes([0.855, 0.15, 0.020, 0.275])

  # plot eigenvalue time series
  ax1.plot(t_P, np.real(vals)*system.Porb, color='b')
  ax1p.plot(t_P, np.imag(vals)/system.Oorb, color='g')
  plt.setp(ax1.get_xticklabels(), visible=False)

  # plot eigenvector surface plot
  ### vecs1
  cmap_min = np.floor(round(np.amin(vecs1[vecs1>-np.infty]),5))
  cmap_max = np.ceil(round(np.amax(vecs1[vecs1<np.infty]),5))
  cmap_norm=matplotlib.colors.Normalize(vmin=cmap_min, vmax=cmap_max)
  vec_plot(ax2, vecs1, t_P, cmap=cmap, cmap_norm=cmap_norm)
  cb2 = matplotlib.colorbar.ColorbarBase(ax2_cb, cmap=cmap, norm=cmap_norm, orientation='vertical')
  plt.setp(ax2.get_xticklabels(), visible=False)

  ### vecs2
  cmap_min = np.floor(round(np.amin(vecs2[vecs2>-np.infty]),5))
  cmap_max = np.ceil(round(np.amax(vecs2[vecs2<np.infty]),5))
  cmap_norm = matplotlib.colors.Normalize(vmin=cmap_min, vmax=cmap_max)
  vec_plot(ax3, vecs2, t_P, cmap=cmap, cmap_norm=cmap_norm)
  cb3 = matplotlib.colorbar.ColorbarBase(ax3_cb, cmap=cmap, norm=cmap_norm, orientation='vertical')

  ax1.set_ylabel(r"$\mathbb{R}\{s\} \cdot P_{\mathrm{orb}}$", color='b')
  ax1p.set_ylabel(r"$\mathbb{I}\{s\} / \Omega_{\mathrm{orb}}$", color='g')
  ax2.set_ylabel(r"mode No.")
  ax3.set_ylabel(r"mode No.")
  ax3.set_xlabel(r'$t/P_{\mathrm{orb}}$')

  return fig, [ax1,ax1p,ax2,ax3], [cb2, cb3]


#########################
def vec_plot(ax, vecs, t_P, cmap=plt.cm.jet, cmap_norm=None):
  """ helper function for eig_plot that generates the proper checkerboard for vectors over time """
  if not cmap_norm:
    cmap_norm = matplotlib.colors.Normalize(vmin=np.amin(vecs), vmax=np.amax(vecs))

  return ax.imshow(np.transpose(vecs), cmap=cmap, norm=cmap_norm, interpolation='nearest', origin='lower', aspect='auto', extent=(t_P[0], t_P[-1],-0.5,len(vecs[0])-0.5) )

##################################################
def eigval_plot(t_P, eigvals, system):
  """
  plots the trajectories of the eigenvalues in eigvals in the complex plane
    eigvals must be an array: size (len(t_P), N_m)

  eigvals will be sorted within this method
  """
  fig = plt.figure()
  ax1 = plt.subplot(2,1,1)
  ax2 = plt.subplot(2,1,2)

  plt.subplots_adjust(hspace=0.05, wspace=0.05)

  eigenvalues = [[] for i in range(len(eigvals[0]))] # place holder

  # sort list of eigenvalues and form plotting arrays
  for evals in eigvals:
    e = list(evals[:])
    e.sort(key=lambda l: np.imag(l), reverse=True)
    e.sort(key=lambda l: np.real(l), reverse=True)
    for ind, ee in enumerate(e):
      eigenvalues[ind].append(ee)

  for evals in eigenvalues:
    ax1.plot(t_P, np.real(evals)*system.Porb, alpha=0.5)
    ax2.plot(t_P, np.imag(evals)/system.Oorb, alpha=0.5)
  
  ax1.set_ylabel(r'$\mathbb{R}\{s\}\cdot P_{\mathrm{orb}}$')
  plt.setp(ax1.get_xticklabels(), visible=False)
  ax2.set_ylabel(r'$\mathbb{I}\{s\} / \Omega_{\mathrm{orb}}$')
  ax2.set_xlabel(r'$t/P_{\mathrm{orb}}$')
  
  return fig, ax1, ax2

##################################################
def growth_rates(t,q,system,current, verbose=False):
  """
  computes the growth rates of the system at all times in sample. Does this by computing
    (dx/dt)/(x) 
  or 
    (dq/dt)/(q)
  as appropriate. Returns a list of complex numbers with len(list) = N_m

  for plotting, this is used in conjunction with real_imag_plot()
  """
  if verbose: print "\tcomputing time derivatives"
  dqdt_P, N_m, N_p = _compute_time_derivatives(t,q,system,current)
  
  if verbose: print "\tcomputing growth rates"
  s = []
  for modeNo in range(N_m):
    _q = q[modeNo]
    _dqdt_P = dqdt_P[modeNo]

    _s = []

    for ind in range(N_p):
      r, i = _q[ind]
      dr, di = _dqdt_P[ind]
      A2 = r**2 + i**2
      if A2 > 0:
        _s.append( [(r*dr + i*di)/A2, (r*di - i*dr)/(2*np.pi*A2)] ) # WARNING: different normalizations for real and imag parts
      else: # special case of zero divisor
        l = []
        if (r*dr + i*di) > 0:
          l.append( abs(r*dr + i*di)/(r*dr + i*di)*np.infty )
        else:
          l.append( 0 )
        if (r*di - i*dr) > 0:
          l.append( abs(r*di - i*dr)/(r*di - i*dr)*np.inty )
        else:
          l.append( 0 )
        _s.append( l )

    s.append( _s )

  return s

