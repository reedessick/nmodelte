usage=""" contains a class which represents a modal network. Includes basic mode classes """

import math
import numpy as np
from nmode_utils import float_to_scientific
from mode_selection import compute_Ethr, compute_Elin

####################################################################################################
#
#
#                                  basic classes for modal networks
#
#
####################################################################################################
class mode:
  """
  basic class that encapsulates a mode

  data contains: 
    n
    l 
    m 
    w 
    y
    U = [(Ui, hi), (Uj, hj), ...] where i,j are harmonics of the tide (used in forcing frequency for integration)
  """

  def __init__(self, n=False, l=False, m=False, w=False, y=False, U=[]):
    self.mode_type = "generic"
    self.n = n
    self.l = l
    self.m = m
    self.w = w
    self.y = y
    self.U = U
    self.check()

  ###
  def check(self):
    return True
    if self.n < 0:
      raise ValueError, "n=%d < 0"%self.n
    if self.l <= 0:
      raise ValueError, "l=%d <= 0"%self.l
    if abs(self.m) > self.l: 
      raise ValueError, "abs(m=%d) > l=%d"%(self.m, self.l)

  ###
  def get_nlm(self):
    return self.n, self.l, self.m

  ### 
  def get_nlms(self):
    
    if self.w != 0:
      sign = self.w/abs(self.w)
    else:
      sign = +1
    return self.n, self.l, self.m, sign

  ###
  def get_wyU(self):
    return self.w, self.y, self.U

  ###
  def get_nlmwy(self):
    return self.n, self.l, self.m, self.w, self.y

  ### for dumping into pickle
  def into_tuple(self):
    return (self.n, self.l, self.m, self.w, self.y, self.U)

  ### for reading from pickle 
  def from_tuple(self, nlmwyU):
    self.n, self.l, self.m, self.w, self.y, self.U = nlmwyU
    return self

  ###
  def to_str_nlm(self, tuple=False):
    if tuple:
      return "(%d, %d, %d)" % self.get_nlm()
    else:
      return "%d %d %d" % self.get_nlm()

  ###
  def to_str_nlms(self, tuple=False):
    if tuple:
      return "(%d, %d, %d, %s)" % self.get_nlms()
    else:
      return "%d %d %d %d" % self.get_nlms()

  ### for dumping to ASCii
  def to_str_nlmwy(self, tuple=False):
    w, wexp = float_to_scientific(self.w)
    y, yexp = float_to_scientific(self.y)
    if tuple:
      return "(%d, %d, %d, %.14fe%d, %.14fe%d)" % (self.n, self.l, self.m, w, wexp, y, yexp)
    else:
      return "%d %d %d %.14fe%d %.14fe%d" % (self.n, self.l, self.m, w, wexp, y, yexp)

  ### for reading from ASCii
  def from_str_nlmwy(self, string):
    n, l, m, w, y = string.strip().split()
    self.n = int(n)
    self.l = int(l)
    self.m = int(m)
    self.w = float(w)
    self.y = float(y)
    return self

##################################################
class network:
  """
  a class that encapsulates a modal network

  data contains:
    modes : list of mode objects
    K     : list of couplings (sparse matrix)
    modeNoD : dictionary mapping : nlms --> modeNo in modes (where s = +/-1, depeding on frequency sign)
    Oorb  : angular frequency of orbit
    nlm   : [_mode.get_nlm() for _mode in modes]
    wyU   : [_mode.get_wyU() for _mode in modes]
  """

  def __init__(self):
    self.modeNoD = {}
    self.modes = []
    self.K = []
    
    self.nlm = []
    self.wyU = []
   
  ###
  def __len__(self):
    return len(self.modes)
 
  ###
  def _update(self):
    self.modeNoD = dict( (m.get_nlms(), i) for i, m in enumerate(self.modes) )
    self.nlm = [m.get_nlm() for m in self.modes]
    self.wyU = [m.get_wyU() for m in self.modes]

  ###
  def find_G0(self):
    """
    finds the gen-0 set of modes defined as the modes which are parent modes in all their couplings

    NOTE: this includes un-coupled modes which is desirable
    """
    parents = []
    not_parents = []
    for o in range(len(self.modes)):
      if o in not_parents:
        continue
      wo = abs(self.modes[o].get_wyU()[0])
      is_parent=True
      for i,j,_ in self.K[o]:
        wi = abs(self.modes[i].get_wyU()[0])
        wj = abs(self.modes[j].get_wyU()[0])
        if wo < max(wi,wj):
          is_parent=False
          break
        else:
          not_parents.append(i)
          not_parents.append(j)
      if is_parent:
        parents.append( o )

    return parents # already sorted

  ###
  def find_Gip1(self, modeNos_Gi):
    """
    finds the gen-i+1 modes, defined as the modes which are daughter modes in the couplings associated with modeNos_Gi in which modeNos_Gi are parents

    expects modeNos_Gi to be a list of mode numbers
    """
    modesGip1 = []
    couplings = []
    for o in modeNos_Gi:
      wo = abs(self.modes[o].get_wyU()[0])
      for i,j,k in self.K[o]:
        wi = abs(self.modes[i].get_wyU()[0])
        wj = abs(self.modes[j].get_wyU()[0])
        if wo > max(wi,wj):
          modesGip1.append(i)
          modesGip1.append(j)
          couplings.append( (o,i,j,k) )

    return sorted(list(set(modesGip1))), couplings

  ###
  def gens(self):
    """
    decomposes the network into generations.
      returns gens, coups
    """
    gi = self.find_G0()
    gens = []
    coups=[]
    N_m = len(self.modes)
    Nm = len(gi)
    gens.append(gi)
    while N_m > Nm: # continue until all modes are included
      gi, coup = self.find_Gip1(gi)
      Nm += len(gi)
      gens.append(gi)
      coups.append(coup)

    return gens, coups

  ###
  def add_modes(self, add_parents, verbose=True):
    """
    by assumption, these modes have no couplings and should just be fed into the params structure 
    add_parents = [mode1, mode2, ...]

    where mode1, etc are instances of networks.mode() objects
    """
    for mode_o in add_parents:
      nlms_o = mode_o.get_nlms()
      if self.modeNoD.has_key( nlms_o ):
        modeNo_o = self.modeNoD[nlms_o]
        if verbose: print "\t"+str(nlms_o)+" already in system! updating parameters:\n\t\tmodeNo. %d: %s --> %s" % (modeNo_o, self.modes[modeNo_o].get_wyU(), mode_o.get_wyU())
        self.modes[modeNo_o] = mode_o
      else:
        modeNo_o = len(self.modes)
        self.modeNoD[nlms_o] = modeNo_o
        self.modes.append( mode_o )
        self.K.append( [] )
        if verbose: print "\t"+str(mode_o.get_nlms())+" added to system!\n\t\tmode No. %d: %s" % (modeNo_o, mode_o.get_wyU())

    self._update()
    return self

  ###
  def to_triples(self):
    """
    converts the network to a list of triples compatible with network.add_couplings()
    meant to be used when combining networks
      eg:  new_network = copy.deepcopy(network1).add_couplings( network2.to_triples() )
    """
    triples_modeNo = []
    for o, coups in enumerate(self.K):
      for i,j,kappa in coups:
        triples_modeNo.append( (tuple(sorted([o,i,j])), kappa) )

    return [ (self.modes[o], self.modes[i], self.modes[j], kappa) for (o,i,j), kappa in set(triples_modeNo) ]

  ###
  def add_couplings(self, add_triples, verbose=True):
    """ 
    add mode (nlm) to network described by params
    add_triples must have the form [(mode_o, mode_a, mode_b, koab),...]
  
    where all modes are instances of networks.mode() objects
    """
    # ensure that all modes are in the network
    for Mo, Ma, Mb, koab in add_triples:
      # ensure that all modes are in the network
      self.add_modes([Mo, Ma, Mb], verbose=verbose)

      modeNoo = self.modeNoD[Mo.get_nlms()]
      modeNoa = self.modeNoD[Ma.get_nlms()]
      modeNob = self.modeNoD[Mb.get_nlms()]
      # add couplings for each triple
      for indo, coupo in enumerate(self.K[modeNoo]):
        i, j, koij = coupo
        if ( (modeNoa, modeNob) == (i, j) ) or ( (modeNoa, modeNob) == (j, i) ):
          if verbose: print "\ttriple: "+str(modeNoo)+","+str(modeNoa)+","+str(modeNob)+" : coupling already exists, updating k = "+str(koij)+" --> "+str(koab)
          self.K[modeNoo][indo] = (modeNoa, modeNob, koab)

          if (modeNoo == modeNoa) and (modeNoo == modeNob): # only one mode
            pass # everything already updated

          elif (modeNoa == modeNob): # only two modes
            for inda, coupa in enumerate(self.K[modeNoa]):
              i, j, koij = coupa
              if ((modeNoo, modeNob) == (i,j)) or ((modeNoo, modeNob) == (j,i)):
                self.K[modeNoa][inda] = (modeNoo, modeNob, koab)
                break

          elif (modeNoo == modeNoa): # only two modes
            for indb, coupb in enumerate(self.K[modeNob]):
              i, j, koij = coupb
              if ((modeNoo, modeNoa) == (i,j)):
                self.K[modeNob][indb] = (modeNoo, modeNoa, koab)
                break

          elif (modeNoo == modeNob): # only two modes
            for inda, coupa in enumerate(self.K[modeNoa]):
              i, j, koij = coupa
              if ((modeNoo, modeNob) == (i, j)):
                self.K[modeNoa][inda] = (modeNoo, modeNob, koab)
                break

          else: # three distinct modes
            for inda, coupa in enumerate(self.K[modeNoa]):
              i, j, koij = coupa
              if ((modeNoo, modeNob) == (i, j)) or ((modeNoo, modeNob) == (j, i)):
                self.K[modeNoa][inda] = (modeNoo, modeNob, koab)
                break
            for indb, coupb in enumerate(self.K[modeNob]):
              i, j, koij = coupb
              if ((modeNoo, modeNoa) == (i, j)) or ((modeNoo, modeNoa) == (j, i)):
                self.K[modeNob][indb] = (modeNoo, modeNoa, koab)
                break

          break
      else:
        if verbose: print "\ttriple: "+str(modeNoo)+","+str(modeNoa)+","+str(modeNob)+" : adding couping k --> "+str(koab)
        if (modeNoo == modeNoa) and (modeNoo == modeNob): # all three modes identical
          self.K[modeNoo].append( (modeNoa, modeNob, koab))

        elif (modeNoa == modeNob): 
          self.K[modeNoo].append((modeNoa, modeNob, koab))
          self.K[modeNoa].append((modeNoo, modeNob, koab))

        elif (modeNoo == modeNoa):
          self.K[modeNoo].append((modeNoa, modeNob, koab))
          self.K[modeNob].append((modeNoo, modeNoa, koab))

        elif (modeNoo == modeNob):
          self.K[modeNoo].append((modeNoa, modeNob, koab))
          self.K[modeNoa].append((modeNoo, modeNob, koab))

        else:
          self.K[modeNoo].append((modeNoa, modeNob, koab))
          self.K[modeNoa].append((modeNoo, modeNob, koab))
          self.K[modeNob].append((modeNoo, modeNoa, koab))

    self._update()
    return self

  ###
  def remove_modes(self, remove_modes, verbose=True):
    """
    remove modes from this network
    remove_modes must be of the form [(n,l,m,s),(n,l,m,s),...]
    """
    removed_modeNos = []

    for nlms in remove_modes:
      if verbose: print "removing "+str(nlms)+" from network."
      if self.modeNoD.has_key(nlms):
        modeNo = self.modeNoD[nlms]
      else:
        if verbose: print "\t"+str(nlms)+" not in network. skipping..."
        continue

      removed_modeNos.append( modeNo ) # remember which numbers we've removed

      self.modes[modeNo] = None

    # build new lists with only remaining modes
    new_modes = []
    new_K = []
    mode_map = {}
    new_modeNo = 0 
    for modeNo, mode in enumerate(self.modes):
      if mode: # if mode has not been removed, we add it to the new list
        new_modes.append( mode )
        new_K.append( self.K[modeNo] )
        mode_map[modeNo] = new_modeNo # define dictionary for updating K
        new_modeNo += 1 # increment new mode number
      else:
        pass

    # update modeNos in remaining list
    for ind, coup in enumerate(new_K): # update couplings now that all new mode numbers are known
      new_coup = []
      for (i,j,k) in coup:
        if (i not in removed_modeNos) and (j not in removed_modeNos): # only keep the coupling if none of the involved modes have been removed
          new_coup.append( (mode_map[i], mode_map[j], k) )
      new_K[ind] = new_coup

    self.modes = new_modes
    self.K = new_K

    self._update()
    return self

##################################################
class system:
  """
  a class that stores the mode network (networks.network) as well as binary orbital parameters. Used for integration and standard reading/loading from pickle file.
  """

  def __init__(self, Mprim, Mcomp, Rprim, Porb, eccentricity, net=False):
    self.Mprim = Mprim
    self.Mcomp = Mcomp
    self.Rprim = Rprim
    self.Porb = Porb
    self.Oorb = 2*math.pi/Porb
    self.eccentricity = eccentricity
    if net:
      self.network = net
    else:
      self.network = network()

  ###
  def compute_linear_freqs(self, modeNos):
    """
    returns the linear oscillation frequency of the modes in modeNos
      if multichromatic forcing is present (len(U)>1), we select the frequency with highest linear energy
    """
    network = self.network

    freqs = []
    for modeNo in modeNos:
      wo, yo, Uo = network.wyU[modeNo]
      mo = network.nlm[modeNo][-1]
      Elin = -np.infty
      O = wo
      for u, k in Uo:
        o = mo*k*self.Oorb # parent oscillates at -m*k*Oorb
        elin = compute_Elin(o, wo, yo, u)
        if elin > Elin:
          Elin = elin
          O = o
      freqs.append( O )

    return freqs

  ###
  def compute_3mode_freqs(self):
    """
    compute all analytic frequencies for this network
    """
    network = self.network
    modes = {}

    g0 = network.find_G0() # find G0
    for modeNo, O in zip(g0, self.compute_linear_freqs(g0)): # compute linear frequencies and assign them
      modes[modeNo] = O

    gi = g0
    while len(modes.keys()) < len(network):
      gip1, couplings = network.find_Gip1(gi) # get daughters
      Ethrs = []
      for o, i, j, k in couplings: # compute all Ethrs for these daughters
        Oo = modes[o] # parent frequency
        wi, yi, _ = network.wyU[i]
        wj, yj, _ = network.wyU[j]
        Ethrs.append( ( compute_Ethr(Oo, wi, wj, yi, yj, k), (o, i, j, k), (wi,wj,yi,yj) ) ) # compute_Ethr imported from mode_selection
    
      Ethrs.sort(key=lambda l: l[0]) # sort the couplings by dynamical relevance so we pick the most important terms first
    
      for Ethr, (o, i, j, k), (wi,wj,yi,yj) in Ethrs: # iterate through and define each daughter's frequency
        if modes.has_key(i) and modes.has_key(j): # both freqs already defined
          pass
        elif modes.has_key(i): # only one freq defined. The other frequency oscillates to cancel the time-dependence
          modes[j] = -(modes[o] + modes[i])
        elif modes.has_key(j):
          modes[i] = -(modes[o] + modes[j])
        else: # modes oscillate at 3 mode stability solutions (our best guess)
          D = modes[o] + wi + wj
          modes[i] = wi - D/(1+yj/yi)
          modes[j] = wj - D/(1+yi/yj)

      gi = gip1

    return [modes[modeNo] for modeNo in sorted(modes.keys())]
#    return modes # (modeNo, freq)
#    return [(value, key) for key, value in modes.items()]

  ###
  def compute_lin_eq(self, t=0.0, default=1e-20, tcurrent="x"):
    """
    computes the linear equlibrium state at t by choosing the most resonant forcing for each mode and computing the linear amplitude.
    for modes without linear forcing, we set the real and imaginary parts to rand*default separately.
    for multi-chromatic forcing, we choose the most resonant frequency.
    """
    network = self.network
    Nm = len(network)
    q = np.zeros((2*Nm,))
    for modeNo in range(Nm):
      wo, yo, Uo = network.wyU[modeNo]
      if len(Uo): # it is linearly forced
        mo = network.nlm[modeNo][-1]
        Elin = -np.infty
        O = wo
        U = 0.0
        for u, k in Uo:
          o = mo*k*self.Oorb
          elin = compute_Elin(o, wo, yo, u)
          if elin > Elin:
            Elin = elin
            O = o
            U = u

        cd = np.cos((wo-O)*wo*U/((wo-O)**2+yo**2))
        sd = np.sin(yo*wo*U/((wo-O)**2+yo**2))

        if tcurrent == "q":
          p = o*t
          cp = np.cos(p)
          sp = np.sin(p)
          q[2*modeNo:2*modeNo+2] = Elin**0.5 * np.array([cp*cd + sp*sd, -sp*cd + sd*cp])
        elif tcurrent == "x":
          p = (O-wo)*t
          cp = np.cos(p)
          sp = np.sin(p)
          q[2*modeNo:2*modeNo+2] = Elin**0.5 * np.array([cp*cd + sp*sd, -sp*cd + sd*cp])
        else:
          raise ValueError, "unknown tcurrent = %s"%tcurrent

      else: # not linearly forced, so we set it to rand*default
        q[2*modeNo:2*modeNo+2] = np.random.rand(2)*default

    return q

  ###
  def compute_3mode_eq(self, t=0.0, default=1e-20, tcurrent="x"):
    """
    computes the 3mode equilibrium state at t by choosing the G0-G1 triple that minimizes Ethr. Ignores all modes with genNo >1
    for modes that do not participate in the 3mode equilib, the real and imaginary parts are set to rand*default, with rand chosen separately for the real and imag parts

    if there are no couplings between G0-G1, we return the linear equilibrium state 
    """
    network = self.network
    gens, coups = network.gens()

    if not len(coups[0]): # no couplings from parents to next generation -> only one generation
      return self.compute_lin_eq(t=t, default=default)

    from nmode_state import threeMode_equilib # we import this here because it avoids conflicts when importing the module as a whole

    freqs = {}
    for modeNo, O in zip(gens[0], self.compute_linear_freqs(gens[0])): # compute linear frequencies and assign them
      freqs[modeNo] = O
      
    Ethr = np.infty
    for o,i,j,k in coups[0]: # these are couplings with G0 as parent (o) and G1 as children (i,j)
      Oo = freqs[o]
      wi, yi, _ = network.wyU[i]
      wj, yj, _ = network.wyU[j]
      ethr = compute_Ethr(Oo, wi, wj, yi, yj, k)
      if ethr < Ethr: # choose tuple with lowest Ethr
        Ethr = ethr
        best_tuple = (o,i,j,k)

    ### for best tuple
    (o,i,j), (Ao, Ai, Aj), ((so, co), (si, ci), (sj, cj)), (do, di, dj) = threeMode_equilib(best_tuple, freqs[best_tuple[0]], network, default_A=default) # get amplitudes, phases and detunings appropriate for tcurrent == "q"
 
    ### check whether the daughters are non-zero
    if (Ai==0) and (Aj==0):
      return self.compute_lin_eq(self, t=t, default=default, tcurrent=tcurrent)
 
    ### instantiate IC vector
    q = np.random.rand(2*len(network))*default # the vast majority are set to rand*default
    
    ### fill in best triple
    if tcurrent == "q":
      p = (wo-do)*t
      cp = np.cos(p)
      sp = np.sin(p)
      q[2*o:2*o+2] = Ao*np.array([cp*co + sp*so, -sp*co + cp*so])

      p = (wi-di)*t
      cp = np.cos(p)
      sp = np.sin(p)
      q[2*i:2*i+2] = Ai*np.array([cp*ci + sp*si, -sp*ci + cp*si])

      p = (wj-dj)*t
      cp = np.cos(p)
      sp = np.sin(p)
      q[2*j:2*j+2] = Aj*np.array([cp*cj + sp*sj, -sp*cj + cp*sj])

    elif tcurrent == "x":
      p = do*t
      cp = np.cos(p)
      sp = np.sin(p)
      q[2*o:2*o+2] = Ao*np.array([cp*co - sp*so, sp*co + cp*so])

      p = di*t
      cp = np.cos(p)
      sp = np.sin(p)
      q[2*i:2*i+2] = Ai*np.array([cp*ci - sp*si, sp*ci + cp*si])

      p = dj*t
      cp = np.cos(p)
      sp = np.sin(p)
      q[2*j:2*j+2] = Aj*np.array([cp*cj - sp*sj, sp*cj + cp*sj])
    else:
      raise ValueError, "unknown tcurrent = %s"%tcurrent

    return q
