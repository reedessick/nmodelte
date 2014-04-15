usage="""written to provide integration functions"""

import numpy as np
import multiprocessing as mp
import sys
import math
import networks

####################################################################################################
#
#
#                   mode amplitude derivatives (dq/dt)
#                          SINGLE CORE
#
####################################################################################################
def dqdt_no_NLT(t, q, args):
  """ a dimension-ful version of the equations of motion that explicitly neglects non-linear tidal forcing

  expect system to be an instance of networks.system() 

  forcing function has the form: U(t) = \sum_{j} Uj*e^{-i*Oorb*m*harm_j*t} where m : azimuthal mode number
                                                                            harm_j : harmonic of the tide
  """

  # find out the number of modes and instantiate the derivative vector
  dimension, system = args[:2]
  N_m = dimension/2
  dqdt = np.zeros((dimension,), np.float)

  # pull out parameters
  W_Y_U = system.network.wyU
  NLM = system.network.nlm
  K = system.network.K
  Oorb = system.Oorb

  p = t*Oorb # the orbital phase

  # loop over modes, computing contributions to derivative vector for each mode
  for m in range(N_m):
    _m = NLM[m][-1]

    # pull out this mode
    qm_r = q[2*m]
    qm_i = q[2*m+1]

    # pull out parameters
    Wm, Ym, Um = W_Y_U[m]

    # LINEAR CONTRIBUTIONS
    real_lin = -Ym*qm_r + Wm*qm_i 
    imag_lin = -Ym*qm_i - Wm*qm_r 
    for _Um, harm_m in Um: # a list of forcings for different harmonics of the tide
      real_lin += Wm*_Um*math.sin(p*_m*harm_m)
      imag_lin += Wm*_Um*math.cos(p*_m*harm_m)

    # 3mode COUPLING CONTRIBUTIONS
    real_3md = 0.
    imag_3md = 0.
    for i, j, k_mij in K[m]: # loop over couplings for m-th mode
      if i!=m and j!=m:
        qi_r = q[2*i]
        qi_i = q[2*i+1]
        qj_r = q[2*j]
        qj_i = q[2*j+1]
      elif i==m and j!=m:
        qi_r = qm_r
        qi_i = qm_i
        qj_r = q[2*j]
        qj_i = q[2*j+1]
      elif i!=m and j==m:
        qi_r = q[2*i]
        qi_i = q[2*i+1]
        qj_r = qm_r
        qj_i = qm_i
      else:
        qi_r = qj_r = qm_r
        qi_i = qj_i = qm_i

      if i!=j:
        real_3md += 2*k_mij*(qi_r*qj_i + qi_i*qj_r)
        imag_3md += 2*k_mij*(qi_r*qj_r - qi_i*qj_i)
      else:
        real_3md += k_mij*(qi_r*qj_i + qi_i*qj_r)
        imag_3md += k_mij*(qi_r*qj_r - qi_i*qj_i)

    real_3md *= Wm
    imag_3md *= Wm

    dqdt[2*m]   = real_lin + real_3md
    dqdt[2*m+1] = imag_lin + imag_3md

  return dqdt

#########################
def dqdt_no_NLT_jac(t, q, args):
  """ jacobian function for dqdp_no_NLT() """

  dimension, system = args[:2]

  N_m = dimension/2
  jac = np.zeros((dimension, dimension), np.float)

  for modeNoo in range(N_m):
    r_indo = 2*modeNoo
    i_indo = r_indo+1
    wo, yo, _ = system.network.wyU[modeNoo]
    jac[r_indo][r_indo] -= yo
    jac[r_indo][i_indo] += wo
    jac[i_indo][r_indo] -= wo
    jac[i_indo][i_indo] -= yo

    for modeNoi, modeNoj, k in system.network.K[modeNoo]:
      if modeNoi == modeNoj:
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = q[r_indi]
        ii = q[i_indi]

        jac[r_indo][r_indi] += 2*wo*k*ii
        jac[r_indo][i_indi] += 2*wo*k*ri
        jac[i_indo][r_indi] += 2*wo*k*ri
        jac[i_indo][i_indi] -= 2*wo*k*ii
      else:
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = q[r_indi]
        ii = q[i_indi]

        r_indj = 2*modeNoj
        i_indj = r_indj+1
        rj = q[r_indj]
        ij = q[i_indj]

        jac[r_indo][r_indi] += 2*wo*k*ij
        jac[r_indo][i_indi] += 2*wo*k*rj
        jac[i_indo][r_indi] += 2*wo*k*rj
        jac[i_indo][i_indi] -= 2*wo*k*ij

        jac[r_indo][r_indj] += 2*wo*k*ii
        jac[r_indo][i_indj] += 2*wo*k*ri
        jac[i_indo][r_indj] += 2*wo*k*ri
        jac[i_indo][i_indj] -= 2*wo*k*ii

  return jac

####################################################################################################
#
#
#                  transformed mode amplitude derivatives (dx/dt)
#                          SINGLE CORE
#
####################################################################################################
def dxdt_no_NLT(t, x, args):
  """ 
  a dimension-ful version of the equations of motion that explicitly neglects non-linear tidal forcing

  implements transformed equations of motion:  q = x*e^{-i*w*t}  <==> x = q*e^{i*w*t}

  expects system to be an instance of networks.system

  forcing function has the form: U(t) = \sum_{j} Uj*e^{-i*Oorb*m*harm_j*t} where m : azimuthal mode number
                                                                            harm_j : harmonic of the tide
  """

#  print t, x

  # find number of modes and instantiate derivative vector
  dimension, system = args[:2]
  N_m = dimension/2
  dxdt = np.zeros((dimension,), np.float)

  # pull out paramters
  W_Y_U = system.network.wyU
  NLM = system.network.nlm
  K = system.network.K
  Oorb = system.Oorb

  p=Oorb*t

  # loop over modes, computing contributions to derivative vector for each mode
  for m in range(N_m):
    # pull out parameters
    Wm, Ym, Um = W_Y_U[m]
    _m = NLM[m][-1]

    # compute linear contributions. We only use the array calls once, so we don't create local variables
    real_lin = -Ym*x[2*m] 
    imag_lin = -Ym*x[2*m+1] 
    for _Um, harm_m in Um: # a list of forcings for different harmonics of the tide
      real_lin += Wm*_Um*math.sin(harm_m*_m*p-Wm*t)
      imag_lin += Wm*_Um*math.cos(harm_m*_m*p-Wm*t)


    # compute 3mode contributions
    real_3md = 0.
    imag_3md = 0.
    for i, j, k_mij in K[m]: # loop over couplings for the m-th mode
      if i == j:
        xi_r = xj_r = x[2*i]
        xi_i = xj_i = x[2*i+1]
        Wi = Wj = W_Y_U[i][0]
      else:
        xi_r = x[2*i]
        xi_i = x[2*i+1]
        xj_r = x[2*j]
        xj_i = x[2*j+1]
        Wi = W_Y_U[i][0]
        Wj = W_Y_U[j][0]

      rr_ii = xi_r*xj_r - xi_i*xj_i
      ri_ri = xi_r*xj_i + xi_i*xj_r

      sin_Wmij = np.sin((Wm+Wi+Wj)*t)
      cos_Wmij = np.cos((Wm+Wi+Wj)*t)

      if i!=j:
        real_3md += 2*k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md += 2*k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)
      else:
        real_3md +=   k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md +=   k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)

    real_3md *= Wm
    imag_3md *= Wm

    dxdt[2*m]   = real_lin + real_3md
    dxdt[2*m+1] = imag_lin + imag_3md

#  print t, dxdt

  return dxdt

#########################
def dxdt_no_NLT_jac(t, x, args):
  """ jacobian function for dxdp_no_NLT() """

  dimension, system = args[:2]

  N_m = dimension/2
  jac = np.zeros((dimension, dimension), np.float)

  for modeNoo in range(N_m):
    r_indo = 2*modeNoo
    i_indo = r_indo+1
    wo, yo, _ = system.network.wyU[modeNoo]
    jac[r_indo][r_indo] -= yo
    jac[i_indo][i_indo] -= yo

    for modeNoi, modeNoj, k in system.network.K[modeNoo]:
      if modeNoi == modeNoj:
        wi,_, _ = system.network.wyU[modeNoi]
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = x[r_indi]
        ii = x[i_indi]

        cosWT = math.cos((wo+wi+wi)*t)
        sinWT = math.sin((wo+wi+wi)*t)

        rs_ic_i = ri*sinWT - ii*cosWT
        rc_is_i = ri*cosWT + ii*sinWT

        jac[r_indo][r_indi] -= 2*wo*k*rs_ic_i
        jac[r_indo][i_indi] += 2*wo*k*rc_is_i
        jac[i_indo][r_indi] += 2*wo*k*rc_is_i
        jac[i_indo][i_indi] += 2*wo*k*rs_ic_i
      else:
        wi,_,_ = system.network.wyU[modeNoi]
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = x[r_indi]
        ii = x[i_indi]

        wj,_,_ = system.network.wyU[modeNoj]
        r_indj = 2*modeNoj
        i_indj = r_indj+1
        rj = x[r_indj]
        ij = x[i_indj]

        cosWT = math.cos((wo+wi+wj)*t)
        sinWT = np.sin((wo+wi+wj)*t)

        rs_ic_i = ri*sinWT - ii*cosWT
        rc_is_i = ri*cosWT + ii*sinWT

        rs_ic_j = rj*sinWT - ij*cosWT
        rc_is_j = rj*cosWT + ij*sinWT

        jac[r_indo][r_indi] -= 2*wo*k*rs_ic_j
        jac[r_indo][i_indi] += 2*wo*k*rc_is_j
        jac[i_indo][r_indi] += 2*wo*k*rc_is_j
        jac[i_indo][i_indi] += 2*wo*k*rs_ic_j

        jac[r_indo][r_indj] -= 2*wo*k*rs_ic_i
        jac[r_indo][i_indj] += 2*wo*k*rc_is_i
        jac[i_indo][r_indj] += 2*wo*k*rc_is_i
        jac[i_indo][i_indj] += 2*wo*k*rs_ic_i

  return jac

####################################################################################################
#
#
#                   mode amplitude derivatives (dx/dt)
#                          MULTI CORE
#
####################################################################################################
def __dxmdt_no_NLT_mp(M, system, connection):
  """
  computes dxm/dt for modes in M only!

  sends result through connection back to user
  if successful, returns True
  """
  N_m = len(M)
  dxmdt = np.zeros((2*N_m,), np.float)

  # pull out paramters
  W_Y_U = system.network.wyU
  NLM = system.network.nlm
  K = system.network.K
  Oorb = system.Oorb

  while True:
    t, x = connection.recv()
    p = Oorb*t

    for ind_m, m in enumerate(M):
      # pull out parameters
      Wm, Ym, Um = W_Y_U[m]
      _m = NLM[m][-1]

      # compute linear contributions. We only use the array calls once, so we don't create local variables
      real_lin = -Ym*x[2*m]
      imag_lin = -Ym*x[2*m+1]
      for _Um, harm_m in Um: # a list of forcings for different harmonics of the tide
        real_lin += Wm*_Um*math.sin(harm_m*_m*p-Wm*t)
        imag_lin += Wm*_Um*math.cos(harm_m*_m*p-Wm*t)

      # compute 3mode contributions
      real_3md = 0.
      imag_3md = 0.
      for i, j, k_mij in K[m]: # loop over couplings for the m-th mode
        if i == j:
          xi_r = xj_r = x[2*i]
          xi_i = xj_i = x[2*i+1]
          Wi = Wj = W_Y_U[i][0]
        else:
          xi_r = x[2*i]
          xi_i = x[2*i+1]
          xj_r = x[2*j]
          xj_i = x[2*j+1]
          Wi = W_Y_U[i][0]
          Wj = W_Y_U[j][0]

        rr_ii = xi_r*xj_r - xi_i*xj_i
        ri_ri = xi_r*xj_i + xi_i*xj_r

        sin_Wmij = np.sin((Wm+Wi+Wj)*t)
        cos_Wmij = np.cos((Wm+Wi+Wj)*t)

        if i!=j:
          real_3md += 2*k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
          imag_3md += 2*k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)
        else:
          real_3md +=   k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
          imag_3md +=   k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)

      real_3md *= Wm
      imag_3md *= Wm

      dxmdt[2*ind_m] = real_lin + real_3md
      dxmdt[2*ind_m+1] = imag_lin + imag_3md

    connection.send( dxmdt ) # communicate back to parent process

#########################
def dxdt_no_NLT_mp(t, x, args):
  """
  computes the transformed equations of motion with parallelization.

  implements transformed equations of motion:  q = x*e^{-i*w*t}  <==> x = q*e^{i*w*t}

  expects system to be an instance of networks.system()

  forcing function has the form: U(t) = \sum_{j} Uj*e^{-i*Oorb*m*harm_j*t} where m : azimuthal mode number
                                                                            harm_j : harmonic of the tide
  """
  dimension, Msets, conns = args[:3]

  dxdt = np.zeros((dimension,), np.float)

  # iterate over connections and send out data
  for conn in conns:
    conn.send( (t,x) )

  # iterate over connections and get data
  for conn, Mset in zip(conns, Msets):
    dxMdt = conn.recv() # recieve output from process

    ### assign data via np.array (faster than iteration?)
    dxdt[2*Mset[0]:2*Mset[-1]+2] = dxMdt

  return dxdt

####################################################################################################

def __dxmdt_no_NLT_p(M, t, x, system, connection):
  """
  computes dxm/dt for modes in M only!

  sends result through connection back to user
  if successful, returns True
  """
  N_m = len(M)
  dxmdt = np.zeros((2*N_m,), np.float)

  # pull out paramters
  W_Y_U = system.network.wyU
  NLM = system.network.nlm
  K = system.network.K
  Oorb = system.Oorb

  p=Oorb*t

  for ind_m, m in enumerate(M):
    # pull out parameters
    Wm, Ym, Um = W_Y_U[m]
    _m = NLM[m][-1]

    # compute linear contributions. We only use the array calls once, so we don't create local variables
    real_lin = -Ym*x[2*m]
    imag_lin = -Ym*x[2*m+1]
    for _Um, harm_m in Um: # a list of forcings for different harmonics of the tide
      real_lin += Wm*_Um*math.sin(harm_m*_m*p-Wm*t)
      imag_lin += Wm*_Um*math.cos(harm_m*_m*p-Wm*t)

    # compute 3mode contributions
    real_3md = 0.
    imag_3md = 0.
    for i, j, k_mij in K[m]: # loop over couplings for the m-th mode
      if i == j:
        xi_r = xj_r = x[2*i]
        xi_i = xj_i = x[2*i+1]
        Wi = Wj = W_Y_U[i][0]
      else:
        xi_r = x[2*i]
        xi_i = x[2*i+1]
        xj_r = x[2*j]
        xj_i = x[2*j+1]
        Wi = W_Y_U[i][0]
        Wj = W_Y_U[j][0]

      rr_ii = xi_r*xj_r - xi_i*xj_i
      ri_ri = xi_r*xj_i + xi_i*xj_r

      sin_Wmij = np.sin((Wm+Wi+Wj)*t)
      cos_Wmij = np.cos((Wm+Wi+Wj)*t)

      if i!=j:
        real_3md += 2*k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md += 2*k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)
      else:
        real_3md +=   k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md +=   k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)

    real_3md *= Wm
    imag_3md *= Wm

    dxmdt[2*ind_m] = real_lin + real_3md
    dxmdt[2*ind_m+1] = imag_lin + imag_3md

  connection.send( dxmdt ) # communicate back to parent process
  return True

#########################
def dxdt_no_NLT_p(t, x, args):
  """
  computes the transformed equations of motion with parallelization.

  implements transformed equations of motion:  q = x*e^{-i*w*t}  <==> x = q*e^{i*w*t}

  expects system to be an instance of networks.system()

  forcing function has the form: U(t) = \sum_{j} Uj*e^{-i*Oorb*m*harm_j*t} where m : azimuthal mode number
                                                                            harm_j : harmonic of the tide
  """
  dimension, system, Msets = args[:3]

  dxdt = np.zeros((dimension,), np.float)

  # holders for process identification
  connections = []
  processes = []

  if len(Msets) == 1:
    return dxdt_no_NLT(t, x, system)

  # iterate over modes and launch parallel jobs
  for Mset in Msets:
    con1, con2 = mp.Pipe()
    proc = mp.Process(target=__dxmdt_no_NLT_p, args=(Mset, t, x, system, con2))
    proc.start()
    connections.append(con1)
    processes.append(proc)

  # wait for all parallel jobs to finish and get data
  for ind, Mset in enumerate(Msets):
    processes[ind].join() # wait for process to finish
    dxMdt = connections[ind].recv() # recieve output from process
    
    ### assign data via np.array (faster than iteration?)
    dxdt[2*Mset[0]:2*Mset[-1]+2] = dxMdt

    ### assign data via iteration
#    for ind_m, m in enumerate(Mset):
#      dxdt[2*m] = dxMdt[2*ind_m] # map output into final array
#      dxdt[2*m+1] = dxMdt[2*ind_m+1]

  return dxdt

#########################
def __dxmdt_no_NLT_p_jac(M, t, x, args, connection):
  """
  computes jacobian for only a subset of modes

  sends result through connection back to user
  if successful, returns True
  """
  dimension, system = args[:2]
  jacM = np.zeros((2*len(M), dimension), np.float)

  for indo, modeNoo in enumerate(M):
    r_Mindo = 2*indo # row index for jacM
    i_Mindo = r_indo+1 # column index

    r_indo = 2*modeNoo
    i_indo = r_indo+1

    wo, yo, _ = system.network.wyU[modeNoo]
    jacM[r_Mindo][r_indo] -= yo
    jacM[i_Mindo][i_indo] -= yo

    for modeNoi, modeNoj, k in system.network.K[modeNoo]:
      if modeNoi == modeNoj:
        wi,_, _ = system.network.wyU[modeNoi]
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = x[r_indi]
        ii = x[i_indi]

        cosWT = math.cos((wo+wi+wi)*t)
        sinWT = math.sin((wo+wi+wi)*t)

        rs_ic_i = ri*sinWT - ii*cosWT
        rc_is_i = ri*cosWT + ii*sinWT

        jacM[r_Mindo][r_indi] -= 2*wo*k*rs_ic_i
        jacM[r_Mindo][i_indi] += 2*wo*k*rc_is_i
        jacM[i_Mindo][r_indi] += 2*wo*k*rc_is_i
        jacM[i_Mindo][i_indi] += 2*wo*k*rs_ic_i
      else:
        wi,_,_ = system.network.wyU[modeNoi]
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = x[r_indi]
        ii = x[i_indi]

        wj,_,_ = system.network.wyU[modeNoj]
        r_indj = 2*modeNoj
        i_indj = r_indj+1
        rj = x[r_indj]
        ij = x[i_indj]

        cosWT = math.cos((wo+wi+wj)*t)
        sinWT = math.sin((wo+wi+wj)*t)

        rs_ic_i = ri*sinWT - ii*cosWT
        rc_is_i = ri*cosWT + ii*sinWT

        rs_ic_j = rj*sinWT - ij*cosWT
        rc_is_j = rj*cosWT + ij*sinWT

        jacM[r_Mindo][r_indi] -= 2*wo*k*rs_ic_j
        jacM[r_Mindo][i_indi] += 2*wo*k*rc_is_j
        jacM[i_Mindo][r_indi] += 2*wo*k*rc_is_j
        jacM[i_Mindo][i_indi] += 2*wo*k*rs_ic_j

        jacM[r_Mindo][r_indj] -= 2*wo*k*rs_ic_i
        jacM[r_Mindo][i_indj] += 2*wo*k*rc_is_i
        jacM[i_Mindo][r_indj] += 2*wo*k*rc_is_i
        jacM[i_Mindo][i_indj] += 2*wo*k*rs_ic_i


  connection.send( jacM )
  return True

#########################
def dxdt_no_NLT_p_jac(t, x, args):
  """ jacobian function for dxdt_no_NLT_p() , also implemented through parallelization"""

  dimension, system, Msets = args[:3]

  jac = np.zeros((dimension,dimension), np.float)

  # holders for process identification
  connections = []
  processes = []

  if len(Msets) == 1:
    return dxdt_no_NLT_jac(t, x, system)

  # iterate over modes and launch parallel jobs
  for Mset in Msets:
    con1, con2 = mp.Pipe()
    proc = mp.Process(target=__dxmdt_no_NLT_p_jac, args=(Mset, t, x, (dimension,system), con2))
    proc.start()
    connections.append(con1)
    processes.append(proc)

  # wait for all parallel jobs to finish and get data
  for ind, Mset in enumerate(Msets):
    processes[ind].join() # wait for process to finish
    jacM = connections[ind].recv() # recieve output from process
    for ind_m, m in enumerate(Mset):
      jac[2*m,:] = jacM[2*ind_m,:] # map output into final array
      dxdt[2*m+1,:] = jacM[2*ind_m+1,:]

  return jac


##################################################
def __dxmdt_no_NLT_mpi(dxmdt, M, t, x, Oorb, W_Y_U, NLM, K):
  """
  computes dxm/dt for modes in M only!
  """
#  print M[0], t, x

  p=Oorb*t

  for m in M:
    # pull out parameters
    Wm, Ym, Um = W_Y_U[m]
    _m = NLM[m][-1]

    # compute linear contributions. We only use the array calls once, so we don't create local variables
    real_lin = -Ym*x[2*m]
    imag_lin = -Ym*x[2*m+1]
    for _Um, harm_m in Um: # a list of forcings for different harmonics of the tide
      real_lin += Wm*_Um*math.sin(harm_m*_m*p-Wm*t)
      imag_lin += Wm*_Um*math.cos(harm_m*_m*p-Wm*t)

    # compute 3mode contributions
    real_3md = 0.
    imag_3md = 0.
    for i, j, k_mij in K[m]: # loop over couplings for the m-th mode
      if i == j:
        xi_r = xj_r = x[2*i]
        xi_i = xj_i = x[2*i+1]
        Wi = Wj = W_Y_U[i][0]
      else:
        xi_r = x[2*i]
        xi_i = x[2*i+1]
        xj_r = x[2*j]
        xj_i = x[2*j+1]
        Wi = W_Y_U[i][0]
        Wj = W_Y_U[j][0]

      rr_ii = xi_r*xj_r - xi_i*xj_i
      ri_ri = xi_r*xj_i + xi_i*xj_r

      sin_Wmij = np.sin((Wm+Wi+Wj)*t)
      cos_Wmij = np.cos((Wm+Wi+Wj)*t)


      if i!=j:
        real_3md += 2*k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md += 2*k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)
      else:
        real_3md +=   k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md +=   k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)

    real_3md *= Wm
    imag_3md *= Wm

    dxmdt[2*m]   = real_lin + real_3md
    dxmdt[2*m+1] = imag_lin + imag_3md

#  print M[0], t, dxmdt

  return dxmdt

#########################
def dxdt_no_NLT_mpi(t, x, args):
  """
  an implementation of dxdt_no_NLT using parallelization via mpi4py. 
  """
  dimension, snds, rcvs, dxmdt = args[:4]

  ### final vector (initialized to zeros)
  dxdt = np.zeros(dimension, dtype="d")

  ### iterate over children and send them t, x
  ### IMPORTANT: may be able to speed things up by combining these into one send!
#  print "parent", x
  for sndt, sndx in snds:
    sndt.Start() # send t to child
    sndx.Start() # send x to child

  ### iterate over children and receive dxmdt from each child
  ### add that contribution to dxdt
  for rcv in rcvs:
    rcv.Start() # request dxmdt from child
    rcv.Wait()
    dxdt += dxmdt

  return dxdt

#########################
def __dxmdt_no_NLT_mpi_jac(dxmdt, M, t, x, Oorb, W_Y_U, NLM, K):
  """ computes jacobian for modes in M only! """
  print "write _dxdmdt_no_NLT_mpi_jac"
  sys.exit()

#########################  
def dxdt_no_NLT_mpi_jac(t, x, args):
  """ computes jacobian for dxdt_no_NLT_mpi_jac(), also implemented through mpi """

  dimension, snds, rcvs, dxdmt = args[:4]
  print "write dxdt_no_NLT_mpi_jac"
  sys.exit()

####################################################################################################
#
#
#                  transformed mode amplitude derivatives (dx/dt)
#                          SINGLE CORE
#
####################################################################################################
def dxdt_no_NLT_withPHI(t, x, args):
  """ 
  a dimension-ful version of the equations of motion that explicitly neglects non-linear tidal forcing

  implements transformed equations of motion:  q = x*e^{-i*w*t}  <==> x = q*e^{i*w*t}

  expects system to be an instance of networks.system

  we expect the forcing to be of the form U(t) = U (a/r)^(l+1) exp^(-i*phi)
  the orbital separation (r) is calculated from the system parameters and phi

  WARNING: we attempt to take the first element of the list of forcing coefficients (in the hansen expansion) as "U"
    when building such a list for use with this integrator, that should be kept in mind.
  """
  dimension, system = args[:2]
  dxdt = np.zeros((dimension,), np.float)

  # separate mode amplitudes from orbital phase
  phi = x[-1]
  x = x[:-1]

  # find number of modes and instantiate derivative vector
  N_m = dimension/2

  # pull out paramters
  W_Y_U = system.network.wyU
  NLM = system.network.nlm
  K = system.network.K
  ec = system.eccentricity

  ec_cosphi = ec*math.cos(phi)
  a_r = (1+ec_cosphi)/(1-ec**2) # orbital semi-major-axis divided by the orbital separation 

  ## dphi/dt
  dxdt[-1] = system.Oorb * (1+ec_cosphi)**2 / (1-ec**2)**3 # rate of change of phi with time

  ## dx_m/dt
  # loop over modes, computing contributions to derivative vector for each mode
  for m in range(N_m):
    # pull out parameters
    Wm, Ym, Um = W_Y_U[m]

    # compute linear contributions. We only use the array calls once, so we don't create local variables
    real_lin = -Ym*x[2*m]
    imag_lin = -Ym*x[2*m+1]
    if len(Um) > 0:
      _l, _m = NLM[m][1:]
      _Um, harm_m = Um[0]
      real_lin += Wm*_Um*a_r**(_l+1)*math.sin(_m*phi-Wm*t)
      imag_lin += Wm*_Um*a_r**(_l+1)*math.cos(_m*phi-Wm*t)

    # compute 3mode contributions
    real_3md = 0.
    imag_3md = 0.
    for i, j, k_mij in K[m]: # loop over couplings for the m-th mode
      if i == j:
        xi_r = xj_r = x[2*i]
        xi_i = xj_i = x[2*i+1]
        Wi = Wj = W_Y_U[i][0]
      else:
        xi_r = x[2*i]
        xi_i = x[2*i+1]
        xj_r = x[2*j]
        xj_i = x[2*j+1]
        Wi = W_Y_U[i][0]
        Wj = W_Y_U[j][0]

      rr_ii = xi_r*xj_r - xi_i*xj_i
      ri_ri = xi_r*xj_i + xi_i*xj_r

      sin_Wmij = np.sin((Wm+Wi+Wj)*t)
      cos_Wmij = np.cos((Wm+Wi+Wj)*t)

      if i!=j:
        real_3md += 2*k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md += 2*k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)
      else:
        real_3md +=   k_mij*(ri_ri*cos_Wmij - rr_ii*sin_Wmij)
        imag_3md +=   k_mij*(ri_ri*sin_Wmij + rr_ii*cos_Wmij)

    real_3md *= Wm
    imag_3md *= Wm

    dxdt[2*m]   = real_lin + real_3md
    dxdt[2*m+1] = imag_lin + imag_3md

  return dxdt

#########################
def dxdt_no_NLT_withPHI_jac(t, x, system):
  """ jacobian function for dxdp_no_NLT_withPHI() """

  dimension = len(x)+1
  N_m = dimension/2 # integer division takes care of the extra 1 (for phi)
  jac = np.zeros((dimension, dimension), np.float)

  ec = system.eccentricity
  ec_cosphi = ec*math.cos(phi)
  ec_sinphi =  ec*math.sin(phi)

  a_r = (1+ec_cosphi)/(1-ec**2) # orbital semi-major-axis divided by the orbital separation
  da_rdphi = -ec_sinphi/(1-ec**2) # partial of a_r with respect to phi

  for modeNoo in range(N_m):
    r_indo = 2*modeNoo
    i_indo = r_indo+1
    wo, yo, Uo = system.network.wyU[modeNoo]
    jac[r_indo][r_indo] -= yo
    jac[i_indo][i_indo] -= yo

    if len(Uo) > 0:
      _l, _m = system.network.nlm[modeNo][1:]
      _Uo,_ = Um[0]
      sin_mp_wt = math.sin(_m*phi-wo*t)
      cos_mp_wt = math.cos(_m*phi-wo*t)

      jac[r_indo][-1] += wo*_Uo*a_r**(_l+1)*_m*cos_mp_wt + wo*_Uo*(_l+1)*a_r**_l*sin_mp_wt*da_rdphi # dependence of forcing on phi
      jac[i_indo][-1] -= wo*_Uo*a_r**(_l+1)*_m*sin_mp_wt + wo*_Uo*(_l+1)*a_r**_l*cos_mp_wt*da_rdphi

    for modeNoi, modeNoj, k in system.network.K[modeNoo]:
      if modeNoi == modeNoj:
        wi,_, _ = system.network.wyU[modeNoi]
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = x[r_indi]
        ii = x[i_indi]

        cosWT = math.cos((wo+wi+wi)*t)
        sinWT = math.sin((wo+wi+wi)*t)

        rs_ic_i = ri*sinWT - ii*cosWT
        rc_is_i = ri*cosWT + ii*sinWT

        jac[r_indo][r_indi] -= 2*wo*k*rs_ic_i
        jac[r_indo][i_indi] += 2*wo*k*rc_is_i
        jac[i_indo][r_indi] += 2*wo*k*rc_is_i
        jac[i_indo][i_indi] += 2*wo*k*rs_ic_i
      else:
        wi,_,_ = system.network.wyU[modeNoi]
        r_indi = 2*modeNoi
        i_indi = r_indi+1
        ri = x[r_indi]
        ii = x[i_indi]

        wj,_,_ = system.network.wyU[modeNoj]
        r_indj = 2*modeNoj
        i_indj = r_indj+1
        rj = x[r_indj]
        ij = x[i_indj]

        cosWT = math.cos((wo+wi+wj)*t)
        sinWT = math.sin((wo+wi+wj)*t)

        rs_ic_i = ri*sinWT - ii*cosWT
        rc_is_i = ri*cosWT + ii*sinWT

        rs_ic_j = rj*sinWT - ij*cosWT
        rc_is_j = rj*cosWT + ij*sinWT

        jac[r_indo][r_indi] -= 2*wo*k*rs_ic_j
        jac[r_indo][i_indi] += 2*wo*k*rc_is_j
        jac[i_indo][r_indi] += 2*wo*k*rc_is_j
        jac[i_indo][i_indi] += 2*wo*k*rs_ic_j

        jac[r_indo][r_indj] -= 2*wo*k*rs_ic_i
        jac[r_indo][i_indj] += 2*wo*k*rc_is_i
        jac[i_indo][r_indj] += 2*wo*k*rc_is_i
        jac[i_indo][i_indj] += 2*wo*k*rs_ic_i

  jac[-1][-1] = -2*system.Oorb*ec_sinphi*(1-ec_cosphi) / (1-ec**2)**3 # d(dphi/dt)/dphi

  return jac

