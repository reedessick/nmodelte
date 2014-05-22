usage=""" written to encapsulate ggg couplings in sun-like stars """

import mode_selection as ms
import gmodes as gm
import copy, math, sys, os, sympy
import numpy as np
import multiprocessing as mp
import sympy.physics.wigner as sympy_phys_wigner
from sympy.core import sympify
import glob

####################################################################################################
#
#
#                       automatic file name generation
#
#
####################################################################################################
def generate_filename(metric, no, lo, mo, absO, min_l, max_l, min_w, max_w):
  """
  generates a standard catalog filename for ggg couplings
  """
  filename = "%s__no=%d_lo=%d_mo=%d__absO=%fe-05__min-l=%d_max-l=%d__min-w=%.6fe-05_max-w=%.6fe-05.ctg" % (metric, no, lo, mo, absO*1e5, min_l, max_l, min_w*1e5, max_w*1e5)

  return filename

##################################################
def parse_filename(filename):
  """
  interprets a standard catalog filename for ggg couplings
  """
  metric, pro_nlm, pro_absO, pro_l, pro_w = filename.split("/")[-1].split("__")

  no, lo, mo = [ int(f.split("=")[-1]) for f in pro_nlm.split("_") ]
  absO = float(pro_absO.split("=")[-1])
  min_l, max_l = [ int(f.split("=")[-1]) for f in pro_l.split("_") ]
  min_w, max_w = [ float(f.strip(".ctg").split("=")[-1]) for f in pro_w.split("_") ]

  return metric, no, lo, mo, absO, min_l, max_l, min_w, max_w

##################################################
def check_filename(filename, metric, no, lo, mo, absO, min_l, max_l, min_w, max_w, verbose=False):
  """
  checks to see whether we can use filename for this particular data set. Requires
  
  abs(mo) == abs(_mo)
  absO == _absO (to within 1e-6)
  lo == _lo
  there must be some overlap between [min_l, max_l] and [_min_l, _max_l]
  there must be some overlap between [min_w, max_w] and [_min_w, _max_w]

  no != _no is fine (we re-normalize elsewhere).
  """
  absO = float( "%fe-5" % (absO*1e5) ) # done so that it will be compatible with the string in the filename
  min_w = float( "%fe-5" % (min_w*1e5) )
  max_w = float( "%fe-5" % (max_w*1e5) )
  try:
    _metric, _no, _lo, _mo, _absO, _min_l, _max_l, _min_w, _max_w = parse_filename(filename)

    metric_test = (metric == _metric)
    absO_test = (absO == _absO)
    lo_test = (lo == _lo)
    absm_test = (abs(mo) == abs(_mo))
   
    lmin_test = (min_l <= _min_l) and (_min_l <= max_l)
    lmax_test = (min_l <= _max_l) and (_max_l <= max_l)

    wmin_test = (min_w <= _min_w) and (_min_w <= max_w)
    wmax_test = (min_w <= _max_w) and (_max_w <= max_w)

    if verbose:
      print "metric_test : " + str(metric_test)
      if not metric_test:
        print "\t metric : %s\n\t_metric : %s" % (metric, _metric)
      print "absO_test : " + str(absO_test)
      if not absO_test:
        print "\t absO : %f\n\t_absO : %f" % (absO, _absO)
      print "lo_test : " + str(lo_test)
      if not lo_test:
        print "\t lo : %d\n\t_lo : %d" % (lo, _lo)
      print "absm_test : " + str(absm_test)
      if not absm_test:
        print "\t absm : %d\n\t_absm : %d" % (absm, _absm)
      print "lmin_test : " + str(lmin_test)
      if not lmin_test:
        print "\t min_l : %d\n\t_min_l : %d" % (min_l, _min_l)
      print "lmax_test : " + str(lmax_test)
      if not lmax_test:
        print "\t max_l : %d\n\t_max_l : %d" % (max_l, _max_l)
      print "wmin_test : " + str(wmin_test)
      if not wmin_test:
        print "\t min_w : %f\n\t_min_w : %f" % (min_w, _min_w)
      print "wmax_test : " + str(wmax_test)
      if not wmax_test:
        print "\t max_w : %f\n\t_max_w : %f" % (max_w, _max_w)

    return metric_test and absO_test and lo_test and absm_test and lmin_test and lmax_test and wmin_test and wmax_test
  except:
    print "could not understand catalog filename: %s" % filename
    return False


####################################################################################################
#
#
#                               utility functions for ggg couplings
#
#
####################################################################################################
def compute_kabc(la, ma, lb, mb, lc, mc, k_hat=5e4, P=10*86400):
  """ 
  compute 3mode coupling coefficient 
    assume kabc = k_hat * (P/Po)^2 * (T/0.2)   (Weinberg, A67)
    where
      T = ( (2*la+1)*(2*lb+2)*(2*lc+1) / (4*np.pi) )**0.5 * Wigner3j(la,lb,lc;ma,mb;mc) * Wigner3j(la,lb,lc;0,0,0)
      P = long-wavelength mode's period period
      Po = normalization for period dependence
  """
  Po=10*86400 # normalization for period dependence
  # check selection rules:
  if ((la+lb+lc)%2 == 0) and (abs(lb-lc) <= la <= lb+lc) and (abs(la-lc) <= lb <= la+lc) and (abs(lb-la) <= lc <= lb+la) and (ma+mb+mc == 0):
    return k_hat * (P/Po)**2 * ( compute_T(la,ma,lb,mb,lc,mc) /0.2 ) 
  else: # selection rules make kabc vanish
    return 0

#########################
def compute_T(la,ma,lb,mb,lc,mc):
  return float( ( (2*la+1)*(2*lb+2)*(2*lc+1) / (4*np.pi) )**0.5 * sympy.N( sympy_phys_wigner.wigner_3j(la,lb,lc,ma,mb,mc) * sympy_phys_wigner.wigner_3j(la,lb,lc,0,0,0) ) )

##################################################
def renormalize_kabc(kabc, Pold, Pnew):
  """
  renormalizes kabc with new period
  """
  return kabc * (Pnew/Pold)**2

##################################################
def all_possible_partners(mode1, mode2, min_w, max_w, alpha, c, wo):
  """
  computes all possible modes that could be coupled to mode1 and mode2 that satisfy the bounds given.

  this can be used with intercouple_network() to include all possible couplings for a given mode within a specified bandwidth
  """
  l1 = mode1.l
  l2 = mode2.l

  # "m" is determined exactly for mode3
  m3 = -(mode1.m + mode2.m)
  
  # compute bounds on l3 from triangle inequalities
  max_l3 = l1+l2
  l3 = max( abs(m3), abs(l1-l2) )
  if (l3 + l1 + l2) % 2 != 0: # sum must be mod 2
    l3 += 1

  k1_fact = mode1.n*(1+1./l1)**0.5
  k2_fact = mode1.n*(1+1./l2)**0.5
  
  # iterate and build a list of modes
  modes = []
  while l3 <= max_l3: # iterate over all allowed l3

    l3_1 = (1+1./l3)**0.5

    n3 = int(math.ceil( max(1, alpha*l3/max_w, (k2_fact-k1_fact)/l3_1) )) # lower bound on n3, requirement that k1 >= | k2 - k3| and min_w <= w3 <= max_w
    max_n3 = math.floor( min((k2_fact+k1_fact)/l3_1, alpha*l3/min_w) ) # upper bound on n3

    while n3 <= max_n3:
      modes.append( gm.gmode(n=n3, l=l3, m=m3, alpha=alpha, c=c, wo=wo, U=[]) )
      n3 += 1
    l3 += 2 # l3+l1+l2 must be mod 2

  return modes

##################################################
def multiple_collective_instabilities(parent, O, Eo, maxp=1, Nmin=0, Nmax=10000, alpha=4e-3, c=2e-10, wo=1e-5, k_hat=5e4, verbose=False, min_l=False, max_l=False, min_n=False, max_n=False, min_absw=False, max_absw=False):
  """
  systematically searches for all allowed collectively unstable modes coupled to parent
  """
  if verbose:
    import time

  n,l,m,w,_ = parent.get_nlmwy()
  P = 2*np.pi/w
  sgnO = abs(O)/O

  ### look for all allowed l1,m1,l2,m2 combinations
  if verbose: 
    print "looking for all allowed l1,m1,l2,m2 combinations"
  starting_triples = []
  l1 = min_l
  while l1 <= max_l:
    l2 = max(l1, abs(l-l1)) # start no lower than l1 to avoid redundant work
    if (l+l1+l2)%2:
      l2 += 1
    max_l2 = min(max_l, l+l1)
    while l2 <= max_l2:
      ### we approximate the detuning as perfect, in which case the 3mode Ethr is minimized when w1 = w2 = O/2
      n1 = int(math.floor(abs(2*alpha*l1/O)))
      n2 = int(math.ceil(abs(2*alpha*l2/O)))
      if n1**3 * n2**3 < (n1+1)**3 * (n2-1)**3: # check whether we've got the floor/ceil correct
        n1 += 1
        n2 -= 1

      m1 = max(-l1, -(m+l2)) 
      if (l1 == l2): # avoid duplicate work
        if (m == 0):
          max_m1 = 0
        elif (m > 0):
          max_m1 = -1
        else: # mo < 0
          max_m1 = min(l1, l2-m)
          m1 = 1 # avoid duplicate pairs
      else:
        max_m1 = min(l1, l2-m)
      while m1 <= max_m1:
        m2 = -(m+m1)
        k = compute_kabc(l, m, l1, m1, l2, m2, k_hat=k_hat, P=P) # compute coupling

        d1 = gm.gmode(n1, l1, m1, alpha=alpha, c=c, wo=wo) # set up daughter modes
        d1.w *= -sgnO
        d2 = gm.gmode(n2, l2, m2, alpha=alpha, c=c, wo=wo)
        d2.w *= -sgnO

        starting_triples.append( (parent, d1, d2, k) ) # add starting point to list

        m1 += 1
      l2 += 2 # increment by 2 because l+l1+l2 must be event
    l1 += 1

  starting_triples.sort(key=lambda l: l[1].n + l[2].n) # sort by the sum of the daughter n's, which should capture the rough size of the parameter space searched?
  if verbose:
    print "\tfound %d starting points"%len(starting_triples)
    tos = []
  triples = []
#  Nmodes = 0
  procs = []
  for triple_ind, triple in enumerate(starting_triples):
    con1, con2 = mp.Pipe()
    new_filename = "single_collective_instability-%d.txt" % triple_ind
    args = (triple, O, Nmin, Nmax, alpha, c, wo, k_hat, Eo, new_filename, min_l, max_l, min_n, max_n, min_absw, max_absw, con2)
    if verbose:
      print "launching search for collective instability : %s\n\tNmin\t: %d\n\tNmax\t: %d\n\tparent\t\t: %s\n\tdaughter\t: %s\n\tdaughter\t: %s" % (new_filename, Nmin, Nmax, triple[0].to_str_nlmwy(), triple[1].to_str_nlmwy(), triple[2].to_str_nlmwy())
      tos.append( time.time() )

    p = mp.Process(target=single_collective_instability, args=args)
    p.start()
    con2.close() # this way the process is the only thing that can write to con2
    procs.append((p, args, con1))

    while len(procs) >= maxp: # wait
      for ind, (p, _, _) in enumerate(procs):
        if not p.is_alive():
          break
      else:
        continue
      p, p_args, con1 = procs.pop(ind)
      if verbose:
        triple = p_args[0]
        print "receiving from collective instability around\n\tparent\t\t: %s\n\tdaughter\t: %s\n\tdaughter\t: %s" % (triple[0].to_str_nlmwy(), triple[1].to_str_nlmwy(), triple[2].to_str_nlmwy())

      if p.exitcode < 0:
        sys.exit("\njob Failed! with status "+p.exitcode() )

      new_triples, _ = con1.recv()
      triples += new_triples
      if verbose:
        print "done: %f seconds\n\tfound %d triples" % (time.time()-tos.pop(ind), len(new_triples))

  while len(procs):
    for ind, (p, _, _) in enumerate(procs):
      if not p.is_alive():
        break
    else:
      continue
    p, p_args, con1 = procs.pop(ind)
    if verbose:
      triple=p_args[0]
      print "receiving from collective instability around\n\tparent\t\t: %s\n\tdaughter\t: %s\n\tdaughter\t: %s" % (triple[0].to_str_nlmwy(), triple[1].to_str_nlmwy(), triple[2].to_str_nlmwy())

    if p.exitcode < 0:
      sys.exit("\njob Failed! with status "+p.exitcode() )

    new_triples, _ = con1.recv()
    triples += new_triples
    if verbose:
      print "done: %f seconds\n\tfound %d triples" % (time.time()-tos.pop(ind), len(new_triples))

  return triples

##################################################
def single_collective_instability(triple, O, Nmin=0, Nmax=10000, alpha=4e-3, c=2e-10, wo=1e-5, k_hat=5e4, Eo=1e-16, verbose=False, min_l=False, max_l=False, min_n=False, max_n=False, min_absw=False, max_absw=False, connection=None):
  """
  looks for and returns a list of triples that constitute a collectively unstable set of modes. We use as an approximate metric
    N**2 * Eo >= Ethr \\forall pairs with k_o12 != 0.0
  Does this by searching neighboring modes to child_mode1, child_mode2 and incrementally adding them to the set

  returns list of triples compatible with network.add_couplings
  """
  verbose_file = isinstance(verbose, str)
  if verbose_file:
    stdout = open(verbose, "w")
  else:
    import sys
    stdout = sys.stdout

  from collections import defaultdict

  absO = abs(O)
  sgnO = O/absO

  # define triples, which will be returned
  triples = []

  # determine which mode in triple is the parent
  if verbose: print >> stdout, "determining which mode is the parent, setting up included objects, computing initial metric"
  modes = [(abs(mode.w), mode) for mode in triple[:3]]
  modes.sort(key=lambda l: l[0], reverse=True) # largest w first
  parent, child1, child2 = [mode for _, mode in modes]

  nlmo = parent.get_nlm()
  no,lo,mo = nlmo
  P = abs(2*np.pi/parent.w)
  Po=10*86400 # normalization for period dependence
  k_hat = k_hat * (P/Po)**2 / 0.2

  starting_metric = ms.compute_Ethr(-absO, abs(child1.w), abs(child2.w), child1.y, child2.y, triple[-1]) # this is the metric we use. 

  # define included children
  included_nlm = set([child1.get_nlm(), child2.get_nlm()]) # the total set of modes in the collective instability.
  included1 = defaultdict( set ) 
  included1[child1.l].add( child1.n )
  included_modes1 = [child1]
  included_nlm1 = [child1.get_nlm()]
  included2 = defaultdict( set ) 
  included2[child2.l].add( child2.n )
  included_modes2 = [child2]
  included_nlm2 = [child2.get_nlm()]

  # define computed T values
  Ts = {} # [(l,m),(l,m),(l,m)] in a particular order

  # define boarder sets
  ### boarder around child 1
  if verbose: print >> stdout, "setting up boarder1"
  boarder1 = defaultdict( set )
  boarder_modes1 = []
  __extend_boarder(child1, boarder1, boarder_modes1, included1, included_nlm2, included_modes2, alpha, c, wo, k_hat, sgnO, absO, nlmo, Ts, min_l=min_l, max_l=max_l, min_n=min_n, max_n=max_n, min_absw=min_absw, max_absw=max_absw) # delegate extension of boarder

  ### boarder around child2
  if verbose: print >> stdout, "setting up boarder2"
  boarder2 = defaultdict( set )
  boarder_modes2 = []
  __extend_boarder(child2, boarder2, boarder_modes2, included2, included_nlm1, included_modes1, alpha, c, wo, k_hat, sgnO, absO, nlmo, Ts, min_l=min_l, max_l=max_l, min_n=min_n, max_n=max_n, min_absw=min_absw, max_absw=max_absw) 

  ### iterate and include modes into collectively unstable set
  if verbose: print >> stdout, "iterating to fill out included"
  Nmodes = len(included_nlm) # we always start off with 2 modes in the network
  Nboarder1 = len(boarder_modes1)
  Nboarder2 = len(boarder_modes2)

  if (Nboarder1==0) and (Nboarder2==0):
    metric = np.infty
  elif Nboarder1==0:
    metric = boarder_modes2[0][0]
  elif Nboarder2==0:
    metric = boarder_modes1[0][0]
  else:
    metric = min(boarder_modes1[0][0], boarder_modes2[0][0])
  include_first_triple = (starting_metric <= Eo)
  if include_first_triple:
    metric = starting_metric
    triples.append(triple)

  while (Nmodes < Nmax) and ((Nboarder1 > 0) or (Nboarder2 > 0)): 
    
    # determine which set of boarder modes to use
    if (len(boarder_modes1) == 0): # use boarder 2
      this_included = included2
      this_included_modes = included_modes2
      this_included_nlm = included_nlm2
      this_boarder = boarder2
      this_boarder_modes = boarder_modes2
      other_included = included1
      other_included_modes = included_modes1
      other_included_nlm = included_nlm1
      other_boarder = boarder1
      other_boarder_modes = boarder_modes1

    elif (len(boarder_modes2) == 0): # use boarder 1
      this_included = included1
      this_included_modes = included_modes1
      this_included_nlm = included_nlm1
      this_boarder = boarder1
      this_boarder_modes = boarder_modes1
      other_included = included2
      other_included_modes = included_modes2
      other_included_nlm = included_nlm2
      other_boarder = boarder2
      other_boarder_modes = boarder_modes2

    elif boarder_modes1[0][0] < boarder_modes2[0][0]: # boarder1 has the lower metric -> use that
      this_included = included1
      this_included_modes = included_modes1
      this_included_nlm = included_nlm1
      this_boarder = boarder1
      this_boarder_modes = boarder_modes1
      other_included = included2
      other_included_modes = included_modes2
      other_included_nlm = included_nlm2
      other_boarder = boarder2
      other_boarder_modes = boarder_modes2

    else:
      this_included = included2
      this_included_modes = included_modes2
      this_included_nlm = included_nlm2
      this_boarder = boarder2
      this_boarder_modes = boarder_modes2
      other_included = included1
      other_included_modes = included_modes1
      other_included_nlm = included_nlm1
      other_boarder = boarder1
      other_boarder_modes = boarder_modes1

    this_metric, [this_mode, this_pk] = this_boarder_modes[0]

    if (not include_first_triple) and (starting_metric <= this_metric):
      triples.append( triple )
      include_first_triple = True

#    judge = (Nmodes+1)**2 * Eo
#    judge = len(this_pk)**2 * Eo
    judge = Eo

    # determine whether we need to add a mode      
    if (judge >= this_metric): # condtions under which we add the mode

      if verbose: print >> stdout, "adding %s" % this_mode.to_str_nlmwy(tuple=True), "\n\tEo = ", judge, " > ", this_metric, " = this_metric"

      metric = max(this_metric, metric) # we've included this mode so we update the set's metric

      ### we're including this mode, so we remove it from the boarder
      this_boarder_modes.pop(0) # remove the first element
      this_boarder[this_mode.l].remove( this_mode.n ) # remove it from the list

      ### update included objects
      included_nlm.add( this_mode.get_nlm() )
      this_included[this_mode.l].add( this_mode.n )
      this_included_modes.append( this_mode )
      this_included_nlm.append( this_mode.get_nlm() )

      ### update triples
      for Ethr, (partner_ind, k) in this_pk:
        triples.append( (parent, this_mode, other_included_modes[partner_ind], k) )

      ### update other boarder
      __update_boarder_modes(this_mode, len(this_included_modes)-1, other_boarder_modes, k_hat, nlmo, absO, Ts)

      ### update this boarder
      __extend_boarder(this_mode, this_boarder, this_boarder_modes, this_included, other_included_nlm, other_included_modes, alpha, c, wo, k_hat, sgnO, absO, nlmo, Ts, min_l=min_l, max_l=max_l, min_n=min_n, max_n=max_n, min_absw=min_absw, max_absw=max_absw)

    else: # we don't add the mode. So, we extend both boarders in hope of finding more modes
      if verbose: print >> stdout, "no satisfactory mode found.\nextending boarder1. Nboarder1 = ", Nboarder1
      for ind in range(Nboarder1):
        __extend_boarder(boarder_modes1[ind][1][0], boarder1, boarder_modes1, included1, included_nlm2, included_modes2, alpha, c, wo, k_hat, sgnO, absO, nlmo, Ts, min_l=min_l, max_l=max_l, min_n=min_n, max_n=max_n, min_absw=min_absw, max_absw=max_absw)
      if verbose: print >> stdout, "extending boarder2. Nboarder2 = ", Nboarder2
      for ind in range(Nboarder2):
        __extend_boarder(boarder_modes2[ind][1][0], boarder2, boarder_modes2, included2, included_nlm1, included_modes1, alpha, c, wo, k_hat, sgnO, absO, nlmo, Ts, min_l=min_l, max_l=max_l, min_n=min_n, max_n=max_n, min_absw=min_absw, max_absw=max_absw)

      if (Nboarder1 == len(boarder_modes1)) and (Nboarder2 == len(boarder_modes2)):
        if (Nmodes < Nmin):
          if verbose: print >> stdout, "adding %s" % this_mode.to_str_nlmwy(tuple=True), "\n\tNmodes = %d < %d = Nmin" % (Nmodes, Nmin)
          
          metric = max(this_metric, metric) # we've included this mode so we update the set's metric

          ### we're including this mode, so we remove it from the boarder
          this_boarder_modes.pop(0) # remove the first element
          this_boarder[this_mode.l].remove( this_mode.n ) # remove it from the list

          ### update included objects
          included_nlm.add( this_mode.get_nlm() )
          this_included[this_mode.l].add( this_mode.n )
          this_included_modes.append( this_mode )
          this_included_nlm.append( this_mode.get_nlm() )

          ### update triples
          for Ethr, (partner_ind, k) in this_pk:
            triples.append( (parent, this_mode, other_included_modes[partner_ind], k) )

          ### update other boarder
          __update_boarder_modes(this_mode, len(this_included_modes)-1, other_boarder_modes, k_hat, nlmo, absO, Ts)

          ### update this boarder
          __extend_boarder(this_mode, this_boarder, this_boarder_modes, this_included, other_included_nlm, other_included_modes, alpha, c, wo, k_hat, sgnO, absO, nlmo, Ts, min_l=min_l, max_l=max_l, min_n=min_n, max_n=max_n, min_absw=min_absw, max_absw=max_absw)

        else:
          if verbose: print >> stdout, "no new boarder modes found. exiting loop"
          break

    Nmodes = len(included_nlm)
    Nboarder1 = len(boarder_modes1)
    Nboarder2 = len(boarder_modes2)

    if verbose_file:
      stdout.flush()

  if verbose:
    print >> stdout, "Ntriples = ", len(triples)
    print >> stdout, "Nmodes = ", Nmodes, " > ", Nmin, "= Nmin"
    if Nmodes >= Nmax:
      print >> stdout, "Nmodes = ", Nmodes, " >= ", Nmax, "= Nmax"
    elif (Nboarder1 == 0) and (Nboarder2 == 0):
      print >> stdout, "could not find any more allowed modes: Nboarder1 = Nboarder2 = 0"

#  print >> stdout, Nmin, Nmodes, Nmax
#  print >> stdout, (Nmodes+1)**2 * Eo, metric
#  print >> stdout, Nboarder1
#  print >> stdout, Nboarder2

  if connection:
    if verbose:
      print >> stdout, "sending triples"

    connection.send( (triples, Nmodes) )

    if verbose: 
      print >> stdout, "triples sent"
    if verbose_file:
      stdout.close()

    return True

  else:
    if verbose_file:
      stdout.close()

    return triples, Nmodes

#########################
def __collective_metric( pk ):
  """
  computes the ranking metric for an element of a boarder_modes list

  pk is assumed to have the following form
    pk = [ [Ethr, (partner_ind, k)], ... ] and is sorted in order of increasing Ethr

  """
#  return max([Ethr for Ethr, _ in pk])
#  return 1.0*max([Ethr for Ethr, _ in pk])/len(pk)**2
  return np.min( np.array([Ethr for Ethr, _ in pk]) / np.arange(1,len(pk)+1)**2 )

#########################
def __update_boarder_modes(this_mode, this_mode_ind, other_boarder_modes, k_hat, nlmo, absO, Ts):
  """
  updates boarder_modes with all possible couplings to this_mode
  """     
  no, lo, mo = nlmo
  ko_fact = no*(1+1./lo)**0.5

  tn, tl, tm = this_mode.get_nlm() 
  tk_fact = tn*(1+1./tl)**0.5

  lo_minus_tl = abs(lo-tl)
  lo_plus_tl = lo+tl
  tabsw = abs(this_mode.w)

  for ind, [other_metric, (other_mode, other_pk)] in enumerate(other_boarder_modes):
    on,ol,om = other_mode.get_nlm()

    if ((lo+ol+tl)%2 == 0) and (lo_minus_tl <= ol) and (ol <= lo_plus_tl): # check to see whether this_mode can couple to other_mode
      if (ko_fact >= abs(tk_fact-on*(1+1./ol)**0.5) ): # constraint on n

        # compute k 
        key = [(lo,mo), (other_mode.l,other_mode.m), (tl,tm)] ; key.sort(key=lambda L: L[1]) ; key.sort(key=lambda L: L[0]) ; key = tuple(key)
        if Ts.has_key(key):
          T = Ts[key]
        else:
          T = compute_T(lo,mo,other_mode.l,other_mode.m,tl,tm)
          Ts[key] = T
        k = k_hat * T

        Ethr =  ms.compute_Ethr(-absO, tabsw, abs(other_mode.w), this_mode.y, other_mode.y, k)

        ### add index to other_pk
        other_pk = ms.insert_in_place( [Ethr, (this_mode_ind, k)], other_pk ) # we've appended it to _included_modes

        ### update metric
#        other_boarder_modes[ind][0] = __collective_metric( other_pk )
        other_boarder_modes[ind][0] = ms.compute_collE( [ethr for ethr, _ in other_pk] )

  other_boarder_modes.sort(key=lambda L: L[0]) # can't get rid of this sort, but hopefully quicksort will still run efficiently

  return True

#########################
def __extend_boarder(this_mode, this_boarder, this_boarder_modes, this_included, other_included_nlm, other_included_modes, alpha, c, wo, k_hat, sgnO, absO, nlmo, Ts, min_l=False, max_l=False, min_n=False, max_n=False, min_absw=False, max_absw=False):
  """
  extends the boarder for collective_instability(). modifies lists in place
  """
  no,lo,mo = nlmo
  nlm_this = this_mode.get_nlm()
  ### iterate over all possible new neighboring modes
  for n,l,m, partner_inds in __get_neighbors(nlmo, other_included_nlm, nlm_this, this_boarder, this_included, min_l=min_l, max_l=max_l, min_n=min_n, max_n=max_n, min_absw=min_absw, max_absw=max_absw, alpha=alpha):
    new_mode = gm.gmode(n=n, l=l, m=m, alpha=alpha, c=c, wo=wo, U=[])
    new_absw = abs(new_mode.w)
    new_mode.w = -sgnO*new_absw

    ### iterate over all possible partner modes (we know there is at least one)
    metric = - np.infty 
    pk = []
    for partner_ind in partner_inds:
      partner_mode = other_included_modes[partner_ind]

      # compute k 
      key = [(lo,mo), (partner_mode.l,partner_mode.m), (l,m)] ; key.sort(key=lambda L: L[1]) ; key.sort(key=lambda L: L[0]) ; key = tuple(key)
      if Ts.has_key(key):
        T = Ts[key]
      else:
        T = compute_T(lo,mo,partner_mode.l,partner_mode.m,l,m)
        Ts[key] = T
      k = k_hat * T

      # compute Ethr
      Ethr = ms.compute_Ethr(-absO, new_absw, abs(partner_mode.w), new_mode.y, partner_mode.y, k)

      ### add to list of partner_inds and couplings
      pk = ms.insert_in_place( [Ethr, (partner_ind, k)], pk )

      # compute Ethr
      _metric = 1.0*Ethr/len(partner_inds)**2
#      _metric = Ethr
      if _metric > metric:
        metric = _metric

    ### add mode to boarders1
    this_boarder[new_mode.l].add( new_mode.n )
#    this_boarder_modes = ms.insert_in_place( [__collective_metric( pk ) , (new_mode, pk)], this_boarder_modes )
    this_boarder_modes = ms.insert_in_place( [ms.compute_collE( [ethr for ethr, _ in pk] ) , (new_mode, pk)], this_boarder_modes ) # reduces the overhead for sorting during this step. 

  return True

#########################
def __get_neighbors(nlmo, nlm1_list, nlm2, boarders2, included2, min_l=False, max_l=False, min_n=False, max_n=False, min_absw=False, max_absw=False, alpha=4e-3):
  """
  finds the neighbors to n2,l2,m2 that can couple successfully to no,lo,mo and n1,l1,m1
    assumes the initial coupling is valid (does not check)

  makes sure the new mode is not already included in either boarders2 or included2
  """
  no,lo,mo=nlmo
  n2,l2,m2=nlm2

  nlmk1_list = [(n1,l1,m1,n1*(1+1./l1)**0.5) for n1,l1,m1 in nlm1_list]

  new_modes = []
  ko_fact = no*(1+1./lo)**0.5
  
  for new_l in [l2+2, l2-2]:
    if (included2.has_key( new_l ) and (n2 in included2[new_l])) or (boarders2.has_key( new_l ) and (n2 in boarders2[new_l]) ):
      continue
    if (min_l and (new_l < min_l)) or (max_l and (new_l > max_l)):
      continue
    if (min_absw and (1.0*alpha*new_l/n2 < min_absw)) or (max_absw and (1.0*alpha*new_l/n2 > max_absw)):
      continue

    k_fact = n2*(1+1./new_l)**0.5
    include=[]
    for ind, (n1,l1,m1,k1_fact) in enumerate(nlmk1_list):
      if (new_l >= abs(m2)) and (abs(lo-l1) <= new_l) and (new_l <= lo+l1): # check angular momentum constraints
        if abs(k1_fact - k_fact ) <= ko_fact: # make sure the n values are close enough
          include.append( ind )

    if include != []: new_modes.append( (n2, new_l, m2, include) )

  for new_n in [n2-1, n2+1]: # we already know the angular momentum constraints work
    if (included2.has_key( l2 ) and (new_n in included2[l2])) or (boarders2.has_key( l2 ) and (new_n in boarders2[l2]) ):
      continue
    if (min_n and (new_n < min_n)) or (max_n and (new_n > max_n)):
      continue
    if (min_absw and (1.0*alpha*l2/new_n < min_absw)) or (max_absw and (1.0*alpha*l2/new_n > max_absw)):
      continue

    k_fact = new_n*(1+1./l2)**0.5
    include = []
    for ind, (n1,l1,m1,k1_fact) in enumerate(nlmk1_list):
      if (abs(lo-l1) <= l2) and (l2 <= lo+l1): # make sure we have ang momentum conservation (nlm2 may not couple to everything in nlm1_list)
        if abs(k1_fact - k_fact ) <= ko_fact: # make sure the n values are close enough
          include.append( ind )
        
    if include != []: new_modes.append( (new_n, l2, m2, include) )

  return new_modes

##################################################
def intercouple_network(network, k_hat=5e4, verbose=False):
  """
  decomposes network into generations and computes all possible couplings between generations.

  WARNING! adds the couplings in place!
           
  """
  gens, coups = network.gens() # decompose into generations
  N_g = len(gens)

  triples = []
  for genNo in range(N_g-1): # iterate and compute all possible parent -> child,child couplings (these are the ones that make physical sense)
    gi = gens[genNo]  # current parents
    giP1 = gens[genNo+1] # current daughters

    existing_coups = [tuple(sorted(coup[:3])) for coup in coups[genNo] ]

    for p in gi: # parent modeNo
      parent = network.modes[p]
      no, lo, mo = parent.get_nlm()
      wo, _, _ = parent.get_wyU()

      ko_fact = no*(1+1./lo)**0.5

      for cind, c1 in enumerate(giP1): # child1 modeNo
        child1 = network.modes[c1]
        nc1, lc1, mc1 = child1.get_nlm()

        kc1_fact = nc1*(1+1./lc1)**0.5

        for c2 in giP1[cind:]: # child2 modeNo, start at 'cind' to avoid duplicate work
          if tuple(sorted([p, c1, c2])) in existing_coups:
            if verbose: print "\tcoupling (%d,%d,%d) already exists!" % (p,c1,c2)

          else:
            child2 = network.modes[c2]
            nc2, lc2, mc2 = child2.get_nlm()
            kc2_fact = nc2*(1+1./lc2)**0.5

            ### check triangle inequality, etc, check that wavenumber constraint is satisfied: k_p >= | k_d2 - k_d1|
            if (abs(lc1-lc2) <= lo) and (lo <= lc1+lc2) and ( ((lo+lc1+lc2)%2) == 0 ) and ((mo+mc1+mc2) == 0) and (ko_fact >= abs(kc1_fact - kc2_fact)):
              if verbose: print "\tcomputing coupling (%d,%d,%d)" % (p,c1,c2)
              triples.append( (parent, child1, child2, compute_kabc(lo, mo, lc1, mc1, lc2, mc2, k_hat=k_hat, P=2*np.pi/wo ) ) )

  network.add_couplings(triples, verbose=verbose)

  return network

####################################################################################################
#
#
#                                      ggg coupling lists
#
#
####################################################################################################
class ggg_coupling_list(ms.coupling_list):
  """
  this class handles coupling lists, reading and writing them to/from disk, etc.

  contains basic functionality, like sorting and combining different lists

  data contains:
    filenames : list of filenames associated with this list
    parent_mode : an instance of networks.mode()
    O : forcing freq against which the daughters couple

    num_pairs : the maximum number of pairings read from each file in filenames
      we should not trust this list if we request more couplings than num_pairs

    alpha, c, wo, k_hat : used to define the parameters of the primary star. These are assumed constant for all mode lists added, but this is not checked.
  """


  def __init__(self, alpha, c, wo, k_hat, parent_mode=None):
    self.filenames = []
    self.couplings = []
    self.metric = None
    self.alpha = alpha
    self.c = c
    self.wo = wo
    self.k_hat = k_hat
    self.O = None
    self.parent_mode = parent_mode
    self.num_pairs = -1 # an unphysical number so we don't get confused.

  ### 
  def load_mode_lists(self, metric, filenames, num_pairs=-1, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """
    loads data from sorted files into a sorted mode list of g-modes. Looks for the first num_pairs modes in each list and loads them. 

    the internal list is kept sorted.

    a subset of modes can be specified with the optional arguments

    *** overrides the existing method to use ggg_coupling_list_element
    """

    # check for basic compatibility
    if self.num_pairs == -1:
      self.num_pairs = num_pairs
    elif self.num_pairs < num_pairs:
      sys.exit("already loaded lists with num_pairs=%d. Cannot load more lists with num_pairs=%d>%d" %(self.num_pairs, num_pairs, self.num_pairs))
    if self.metric == None:
      self.metric = metric
    elif self.metric != metric:
      sys.exit("conflict of metric. Mode list already has metric=%s, cannot load lists with metric $s" % (self.metric, metric))

    # set up iteration
    if isinstance(filenames, str):
      filenames = [filenames]
    couplings = [ (c.metric_value, c) for c in self.couplings]

    for filename in filenames:

      f = open(filename, 'r')
      absO = float(f.readline().strip())
      Mo = gm.gmode().from_str_nlmwy(f.readlin().strip())

      # check for forcing freq compatibility
      if self.O == None:
        self.O = absO
      elif self.O != absO:
        sys.exit("incompatible forcing frequencies!")
      # check for parent mode compatibility
      if self.parent_mode == None:
        self.parent_mode = Mo
      elif (self.parent_mode.l != Mo.l) or (abs(self.parent_mode.m) != abs(Mo.m)):
        sys.exit("incompatible parent modes!\n\texisting: %s\n\tnew: %s" % (str(self.parent_mode.gen_nlmwy()), str(Mo.get_nlmwy())))

      self.filenames.append(filename)

      # read in modes
      ind == 0
      _couplings = []
      for line in f:
        if (num_pair != -1) and (ind >= num_pair): # only load num_pair couplings
          break

        if line[0] != "#":
          coupling = ggg_coupling_list_element().from_str( line.strip(), self.alpha, self.c, self.wo ).setup( self.parent_mode ).renorm_koab( self.parent_mode, Mo )
          if ms.check_mode(coupling.dmode1, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w) and ms.check_mode(coupling.dmode2, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w):
            _couplings = ms.insert_in_place( (coupling.metric_value, coupling), _couplings )

      f.close()
      couplings = ms.insert_sortedList_in_place( _couplings, couplings, sorted=True)

    self.couplings = [c for metric_value, c in couplings]

    return self

  ###
  def load_unsorted_mode_lists(self, metric, filenames, num_pairs=-1, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """
    loads data from unsorted files into a sorted mode list. Looks for the first num_pairs modes in each list and loads them. 

    the internal list is kept sorted.

    a subset of modes can be specified with the optional arguments

    *** overrides the existing method to use ggg_coupling_list_element
    """

    # check for basic compatibility
    if self.num_pairs == -1:
      self.num_pairs = num_pairs
    elif self.num_pairs < num_pairs:
      sys.exit("already loaded lists with num_pairs=%d. Cannot load more lists with num_pairs=%d>%d" %(self.num_pairs, num_pairs, self.num_pairs))
    if self.metric == None:
      self.metric = metric
    elif self.metric != metric:
      sys.exit("conflict of metric. Mode list already has metric=%s, cannot load lists with metric $s" % (self.metric, metric))

    # set up iteration
    if isinstance(filenames, str):
      filenames = [filenames]
    couplings = [ (c.metric_value, c) for c in self.couplings]

    for filename in filenames:

      f = open(filename, 'r')
      absO = float(f.readline().strip())
      Mo = gm.gmode().from_str_nlmwy(f.readline().strip())

      # check for forcing freq compatibility
      if self.O == None:
        self.O = absO
      elif self.O != absO:
        print "incompatible forcing frequencies!: \n\texisting: %.9f\n\t new: %.9f\nskipping %s" % (self.O, absO, filename)
        continue
#        sys.exit("incompatible forcing frequencies!")
      # check for parent mode compatibility
      if self.parent_mode == None:
        self.parent_mode = Mo
      elif (self.parent_mode.l != Mo.l) or (abs(self.parent_mode.m) != abs(Mo.m)):
        sys.exit("incompatible parent modes!\n\texisting: %s\n\tnew: %s" % (str(self.parent_mode.get_nlmwy()), str(Mo.get_nlmwy())))

      self.filenames.append(filename)

      # read in modes
      _couplings = []
      for line in f:
        if line[0] != "#":
          coupling = ggg_coupling_list_element().from_str( line.strip(), self.alpha, self.c, self.wo ).setup( self.parent_mode ).renorm_koab( self.parent_mode, Mo )
          if ms.check_mode(coupling.dmode1, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w) and ms.check_mode(coupling.dmode2, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w):
            if num_pairs == -1: # load everything 
              _couplings = ms.insert_in_place( (coupling.metric_value, coupling), couplings )
            else:
              if len(_couplings) < num_pairs: # add the mode
                _couplings = ms.insert_in_place( (coupling.metric_value, coupling), _couplings)
              else:
                if coupling.metric_value < _couplings[-1][0]:
                  _couplings.pop()
                  _couplings = ms.insert_in_place( (coupling.metric_value, coupling), _couplings)

      f.close()
#      for c in _couplings:
#        couplings = ms.insert_in_place( c, couplings )
      couplings = ms.insert_sortedList_in_place( _couplings, couplings, sorted=False )

    self.couplings = [c for metric_value, c in couplings]

    return self

  ###
  def to_triples(self, num_pairs, parent_forcing=False, daughter_forcing=False, Mprim=False, Mcomp=False, Porb=False, eccentricity="none", Rprim=1.):
    """
    converts the mode list to a triples format compatible with networks.network.add_couplings()
    """
    if (self.num_pairs != -1) and (num_pairs > self.num_pairs):
      sys.exit("cannot report more pairs than the maximum read from any given file")

    if parent_forcing or daughter_forcing:
      if not (Mprim and Mcomp and Porb and (eccentricity!="none")):
        sys.exit("must supply Mprim, Mcomp, Porb, eccentricity to compute forcing")
    if parent_forcing:
      self.parent_mode.compute_forcing(Porb, eccentricity, Mprim, Mcomp, Rprim)
    else:
      self.parent_mode.U = []

    triples = []
    for coupling in self.couplings[:num_pairs]:
      if daughter_forcing:
        coupling.set_forcings(Porb, eccentricity, Mprim, Mcomp, Rprim=Rprim)
      else:
        coupling.remove_forcings()

      triples.append( (self.parent_mode, coupling.dmode1, coupling.dmode2, coupling.k) )

    return triples

##################################################
class ggg_minima_coupling_list(ms.coupling_list):
  """
  this class handles minima coupling lists, reading and writing them to/from disk, etc. 
  this list loads in only local minima from file. The best "num_pairs" modes are only computed upon the call "to_triples"

  contains basic functionality, like sorting and combining different lists

  data contains:
    filenames : list of filenames associated with this list
    parent_mode : an instance of networks.mode()
    O : forcing freq against which the daughters couple

    num_pairs : the maximum number of pairings read from each file in filenames
      we should not trust this list if we request more couplings than num_pairs

    alpha, c, wo, k_hat : used to define the parameters of the primary star. These are assumed constant for all mode lists added, but this is 
not checked.
  """

  def __init__(self, alpha, c, wo, k_hat, parent_mode=None):
    self.filenames = []
    self.couplings = []
    self.metric = None
    self.alpha = alpha
    self.c = c
    self.wo = wo
    self.k_hat = k_hat
    self.O = None
    self.parent_mode = parent_mode

  ###
  def load_mode_lists(self, metric, filenames, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """
    delegates to load_unsorted_mode_lists
    """
    return self.load_unsorted_mode_lists(metric, filenames, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)

  ###
  def load_unsorted_mode_lists(self, metric, filenames, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """ 
    Read in local minima recorded in filenames. Reads in all elements of the files (as it must) 
    assumes the alpha, c, wo, k_hat values are consistent
    
    actual "best triples" are only computed when we call "to_triples()"

    metric : [min_Ethr, min_heuristic] (used in _load_mode_minima_list_update)
    """
    
    # check for basic compatibility
    if self.metric == None:
      self.metric = metric
    elif self.metric != metric:
      sys.exit("conflict of metric. Mode list already has metric=%s, cannot load lists with metric $s" % (self.metric, metric))

    # set up iteration
    if isinstance(filenames, str):
      filenames = [filenames]
    couplings = [ (c.metric_value, c) for c in self.couplings]

    for filename in filenames:

      f = open(filename, 'r')
      absO = float(f.readline().strip())
      Mo = gm.gmode().from_str_nlmwy(f.readline().strip())

      # check for forcing freq compatibility
      if self.O == None:
        self.O = absO
      elif self.O != absO:
        print "incompatible forcing frequencies!: \n\texisting: %.9f\n\t new: %.9f\nskipping %s" % (self.O, absO, filename)
        continue
#        sys.exit("incompatible forcing frequencies!")
      # check for parent mode compatibility
      if self.parent_mode == None:
        self.parent_mode = Mo
      elif (self.parent_mode.l != Mo.l) or (abs(self.parent_mode.m) != abs(Mo.m)):
        sys.exit("incompatible parent modes!\n\texisting: %s\n\tnew: %s" % (str(self.parent_mode.get_nlmwy()), str(Mo.get_nlmwy())))

      self.filenames.append(filename)

      # set up iteration
#      _couplings = copy.deepcopy( couplings )
      _couplings = []
      ind = 0
      for line in f:
        if line[0] != "#":
          ind += 1
          coupling = ggg_coupling_list_element().from_str( line.strip(), self.alpha, self.c, self.wo ).setup( self.parent_mode ).renorm_koab( self.parent_mode, Mo )
          if ms.check_mode(coupling.dmode1, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w) and ms.check_mode(coupling.dmode2, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w):
            _couplings = ms.insert_in_place( (coupling.metric_value, coupling), _couplings )
      f.close()

#      for c in _couplings:
#        couplings = ms.insert_in_place( c, couplings )
      couplings = ms.insert_sortedList_in_place( _couplings, couplings, sorted=True )

    self.couplings = [c for metric_value, c in couplings]
    return self

  ###
  def load_minima_list_update(self, metric, _lo, _mo, _na, _la, _ma, _nb, _lb, _mb, cpairsD, couplings, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """
    updates the list of minima. checks both modes a and b agains specified bounds

    metric : [min_Ethr, min_heuristic]
    """
    modea = gm.gmode(_na, _la, _ma, self.alpha, self.c, self.wo)
    modeb = gm.gmode(_nb, _lb, _mb, self.alpha, self.c, self.wo)

    _nlma = modea.get_nlm()
    _nlmb = modeb.get_nlm()
    if not cpairsD.contains( modea, modeb):
      if ms.check_mode(modea, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w) and ms.check_mode(modeb, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w):
        # generate new coupling
        coupling = ggg_coupling_list_element()
        coupling.dmode1 = modea
        coupling.dmode1 = modeb
        # compute koab
        coupling.k = compute_kabc(_lo, _mo, _la, _ma, _lb, _mb, k_hat=self.k_hat, P=2*np.pi/self.parent_mode.w) 
        # compute rank
        if metric == "min_Ethr":
          coupling.metric_value = ms.compute_Ethr(-abs(self.O), modea.w, modeb.w, modea.y, modeb.y, coupling.k)
        elif metric == "min_heuristic":
          coupling.metric_value = ms.compute_heuristic(-abs(self.O), modea.w, modeb.w, modea.y, modeb.y)
        else:
          sys.exit("rank_func=%s not recognized in ggg._load_mode_minima_list_update")

        cpairsD.add( modea, modeb )
        couplings = ms.insert_in_place( (coupling.metric_value, coupling), couplings)
 
    return True

  ###
  def to_triples(self, num_pairs, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False, parent_forcing=False, daughter_forcing=False, Mprim=False, Mcomp=False, Porb=False, eccentricity="none", Rprim=1.):
    """ this is where we build around local minima 
    if supplied, parent_forcing must be accompanied by Mprim, Mcomp, Porb, eccentricity
    same for daughter forcing
    """

    if parent_forcing or daughter_forcing:
      if not (Mprim and Mcomp and Porb and (eccentricity!="none")):
        sys.exit("must supply Mprim, Mcomp, Porb, eccentricity to compute forcing")
    if parent_forcing:
      self.parent_mode.compute_forcing(Porb, eccentricity, Mprim, Mcomp, Rprim)

    no, lo, mo = self.parent_mode.get_nlm()
    cwo3a2 = self.c * self.wo**3 * self.alpha**-2

    cpairsD = _pairsD()
    for coupling in self.couplings:
      cpairsD.add( coupling.dmode1, coupling.dmode2 )

    # pull out the best pairs
    triples = []
    couplings = copy.deepcopy(self.couplings)
    _couplings = [(c.metric_value, c) for c in couplings]
    while len(triples) < num_pairs:
      best_coupling = couplings.pop(0)
      best_coupling.setup( self.parent_mode ) # we don't need to renormalize k because we compute it from the parent mode anyway
      E = best_coupling.metric_value
      mode1 = best_coupling.dmode1
      mode2 = best_coupling.dmode2
      if daughter_forcing:
        best_coupling.compute_forcings(Porb, eccentricity, Mprim, Mcomp, Rprim=Rprim)
      else:
        best_coupling.remove_forcings()
      koab = best_coupling.k
        
      triples.append( (self.parent_mode, mode1, mode2, koab) )

      na, la, ma = mode1.get_nlm()
      nb, lb, mb = mode2.get_nlm()

      ### na+1
      self.load_minima_list_update(self.metric, lo , mo , na+1, la , ma , nb , lb , mb , cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### na-1
      if na >= 2:
        self.load_minima_list_update(self.metric, lo, mo, na-1, la, ma, nb, lb, mb, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### nb+1
      self.load_minima_list_update(self.metric, lo, mo, na, la, ma, nb+1, lb, mb, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### na-1
      if nb >= 2:
        self.load_minima_list_update(self.metric, lo, mo, na+1, la, ma, nb-1, lb, mb, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### la+1
      if (ma > -la) and (la+1 <= lb+lo):
        self.load_minima_list_update(self.metric, lo, mo, na, la+1, ma, nb, lb, mb, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### la-1
      if (la >= 2) and (ma < la) and (abs(lb-lo) <= la-1):
        self.load_minima_list_update(self.metric, lo, mo, na, la-1, ma, nb, lb, mb, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### lb+1
      if (mb > -lb) and (lb+1 <= la+lo):
        self.load_minima_list_update(self.metric, lo, mo, na, la, ma, nb, lb+1, mb, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ## lb-1
      if (lb >= 2) and (mb < lb) and (abs(la-lo) <= lb-1):
        self.load_minima_list_update(self.metric, lo, mo, na, la, ma, nb, lb-1, mb, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### ma+1 ( and mb )
      if (ma < la) and (mb > -lb):
        self.load_minima_list_update(self.metric, lo, mo, na, la, ma+1, nb, lb, mb-1, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
      ### ma-1 ( and mb )
      if (ma > -la) and (mb < lb):
        self.load_minima_list_update(self.metric, lo, mo, na, la, ma-1, nb, lb, mb+1, cpairsD, _couplings, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w)
  
    return triples

##################################################
class _pairsD():
  ###
  def __init__(self):
    self.pairsD = {}

  ###
  def add(self, mode1, mode2 ):
    na, la, ma = mode1.get_nlm()
    nb, lb, mb = mode2.get_nlm()
    add = self.add_which(na,la,ma,nb,lb,mb)
    if add == 'a':
      if self.pairsD.has_key((na,la,ma)): self.pairsD[(na,la,ma)].append( (nb, lb, mb) )
      else: self.pairsD[(na,la,ma)] = [(nb, lb, mb)]
    else:
      if self.pairsD.has_key((nb,lb,mb)): self.pairsD[(nb,lb,mb)].append( (na,la,ma) )
      else: self.pairsD[(nb,lb,mb)] = [(na,la,ma)]
  ###
  def contains(self, mode1, mode2 ):
    na, la, ma = mode1.get_nlm()
    nb, lb, mb = mode2.get_nlm()
    add = self.add_which(na,la,ma,nb,lb,mb)
    if add == 'a':
      if self.pairsD.has_key((na,la,ma)) and ((nb,lb,mb) in self.pairsD[(na,la,ma)]):
        return True
    if add == 'b':
      if self.pairsD.has_key((nb,lb,mb)) and ((na,la,ma) in self.pairsD[(nb,lb,mb)]):
        return True
    return False

  ###
  def add_which(self, na,la,ma,nb,lb,mb):
    if   la < lb: add = 'a'
    elif lb < la: add = 'b'
    else:
      if   na < nb: add = 'a'
      elif nb < na: add = 'b'
      else:
        if   ma < mb: add = 'a'
        else: add = 'b'
    return add

##################################################
class ggg_coupling_list_element(ms.coupling_list_element):
  """
  a class representing the coupling_list elements for ggg couplings. This inherits from mode_selection.coupling_list_element, but ensures that these modes are instances of gmodes.gmode()
  also generates forcing coefficients as modes are read in.

  data contains:
    metric_value : used to sort
    dmode1 : first daughter mode. must be an instance of networks.mode()
    dmode2 : second daughter mode. must be an instance of networks.mode()
    k : coupling coefficient

  coupling_list_element does not know about the parent mode. That's stored by the coupling_list as a whole.
  """

  def from_str(self, string, alpha, c, wo):
    """
    parse dpair from catalog ASCII files
    """
    string = string.strip().split("(")
    string = [item.strip().strip(")") for item in string[:-1] + string[-1].split(")")]

    metric_value, dmode1str, dmode2str, k = string
    
    self.metric_value = float(metric_value)
    try:
      self.k = float(k)
    except:
      print string
      raise ValueError

    n, l, m, w, y = dmode1str.split(",")
    self.dmode1 = gm.gmode(int(n), int(l), int(m), alpha, c, wo)

    n, l, m, w, y = dmode2str.split(",")
    self.dmode2 = gm.gmode(int(n), int(l), int(m), alpha, c, wo)

    return self

  ###
  def renorm_koab(self, new_parent_mode, old_parent_mode):
    """
    re-normalizes the coupling constant for ggg couplings
      old --> new
    this accounts for the freq dependence of the parent mode
    """
    self.k = renormalize_kabc(self.k, 2*math.pi/old_parent_mode.w, 2*math.pi/new_parent_mode.w)
    return self

  ###
  def set_forcings(self, Porb, eccentricity, Mprim, Mcomp, Rprim=1.):
    """
    computes forcings for each daughter mode
    """
    self.dmode1.compute_forcing(Porb, eccentricity, Mprim, Mcomp, Rprim)
    self.dmode1.compute_forcing(Porb, eccentricity, Mprim, Mcomp, Rprim)
    return self

  ###
  def remove_forcings(self):
    """
    removes forcings for each daughter mode
    """
    self.dmode1.U = []
    self.dmode2.U = []

####################################################################################################
#
#
#                             parent generation algorithms
#           returns a lists of modes (gmodes.gmode()) that can be added to a networks.network()
#
####################################################################################################
def compute_parents_detuning(Oorb, bounds, N=1, min_w=0.0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, Porb=False, Mprim=False, Mcomp=False, eccentricity="none", Rprim=1., forcing=False):
  """ 
  generates modes on the fly and returns the first N that obey bounds and fall within the specified bandwidth

  bounds = {l1:[m1_1, m1_2, ..], ...}
  where li is an allowed value for l and mi_k are the allowed values of m for that l
  """
  absO = abs(Oorb)
  signO = absO/Oorb

  if forcing:
    if not (Mprim and Mcomp and Porb and (eccentricity != "none")):
      sys.exit("must supply Mprim, Mcomp, Porb, eccentricity to compute forcing")

  parents = []
  for lo in bounds.keys():
    no = int(max(1, math.ceil(alpha*lo/max_w))) # lower bound on na
    max_no = math.floor(alpha*lo/min_w) # upper bound on na

    while no <= max_no:

      for mo in bounds[lo]:
        tmp_mode = gm.gmode(no, lo, mo, alpha, c, wo)
        signmo = mo/abs(mo)
        D = abs(mo*absO-signmo*tmp_mode.w)

        if len(parents) < N: # always add a mode
          tmp_mode.w *= signO*signmo
          if forcing:
            tmp_mode = tmp_mode.compute_forcing(Porb, eccentricity, Mprim, Mcomp, Rprim)
          parents = ms.insert_in_place( (D, tmp_mode), parents )
          min_D = parents[-1][0]

        elif D < min_D: # mode is better than others. add it
          parents.pop()
          tmp_mode.w *= signO*signmo
          if forcing:
            tmp_mode = tmp_mode.compute_forcing(Porb, eccentricity, Mprim, Mcomp, Rprim)
          parents = ms.insert_in_place( (D, tmp_mode), parents )
          min_D = parents[-1][0]

      no += 1

  return [l[1] for l in parents]

##################################################
def compute_parents_Elin(Oorb, bounds, N=1, min_w=0.0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, Porb=False, Mprim=False, Mcomp=False, eccentricity="none", Rprim=1.):
  """ 
  generates modes on the fly and returns the first N that obey bounds and fall within the specified bandwidth

  *** N is really only a lower bound because we can't distinguish between degenerate parents ***

  bounds = {l1:[m1_1, m1_2, ..], ...}
  where li is an allowed value for l and mi_k are the allowed values of m for that l
  """
  absO = abs(Oorb)
  signO = absO/Oorb

  if not (Mprim and Mcomp and Porb and (eccentricity != "none")):
    sys.exit("must supply Mprim, Mcomp, Porb, eccentricity to compute forcing")

  parents = []
  for lo in bounds.keys():
    no = int(max(1, math.ceil(alpha*lo/max_w))) # lower bound on na
    max_no = math.floor(alpha*lo/min_w) # upper bound on na

    while no <= max_no:

      for mo in bounds[lo]:
        tmp_mode = gm.gmode(no, lo, mo, alpha, c, wo)
        tmp_mode.compute_forcing(Porb, eccentricity, Mprim, Mcomp, Rprim)
        signmo = mo/abs(mo)
        D = -ms.compute_Elin(mo*absO, signmo*tmp_mode.w, tmp_mode.y, tmp_mode.U[0][0]) # we want maximum Elin, which means minimize -Elin

        if len(parents) < N: # always add a mode
          tmp_mode.w *= signO*signmo
          parents = ms.insert_in_place( (D, tmp_mode), parents )
          min_D = parents[-1][0]

        elif D < min_D: # mode is better than others. add it
          parents.pop() # remove worst mode
          tmp_mode.w *= signO*signmo
          parents = ms.insert_in_place( (D, tmp_mode), parents )
          min_D = parents[-1][0]

      no += 1

  return [l[1] for l in parents]

####################################################################################################
#
#
#                     brute force mode list generation methods
#             writes coupling lists to file (quickly), from which they can be read
#           formatting should be compatible with ggg_coupling_list_element.from_str()
####################################################################################################
def compute_pairs_Ethr(parent_mode, O, min_l=1, max_l=100, min_w=0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, filename=False, Emax=np.infty):
  """
  compute 3mode Ethr for daughter pairs. Daughter pairs are generated on the fly, with 
    min_l <= la <= max_l
    min_w <= wa <= max_w

  and all other parameters are limited by the following constraints
    |lo - la| <= lb <= lo + la
    lo + la + lb = even
    mo + ma + mb = 0
    ko >= |ka - kb|

  WARNING! mode parameters are computed internally without reference to compute_w() or compute_y()

  """
  no, lo, mo = parent_mode.get_nlm()
  Po = 2*np.pi/abs(parent_mode.w)
  cwo3a2 = c*wo**3*alpha**-2

  ko_fact = no * (1+1./lo)**0.5

  absO=abs(O)

  if filename:
    f=open(filename, "w")
    print >>f, absO
    print >>f, parent_mode.to_str_nlmwy()
  else:
    print absO
    print parent_mode.to_str_nlmwy()

  la = max(1, min_l) # instantiate la
  while la <= max_l:
    min_na = max(1, int(math.ceil(alpha*la/max_w))) # lower bound on na
    max_na = math.floor(alpha*la/min_w) # upper bound on na
    # demand that lb >= la (computational reasons)
    # demand that (la+lb+lo) is even
    if lo%2 == 0:
      lb = la
    else:
      lb = la + 1
    max_lb = la + lo # enforce triangle inequality
    while lb <= max_lb:
      lb_1 = (1+1./lb)**0.5

      # compute koab first
      kappas = []
      # require that (ma+mb+mo)=0
      ma = max(-la, -lb-mo) # bounds on ma and mb=-(mo+ma)
      if (la == lb): # avoid duplicate work
        if (mo == 0):
          max_ma = 0
        elif (mo > 0):
          max_ma = -1
        else: # mo < 0
          max_ma = min(la, lb-mo)
          ma = 1 # avoid duplicate pairs
      else:
        max_ma = min(la, lb-mo) # bounds on ma and mb=-(mo+ma)
      while ma <= max_ma:
        mb = -mo-ma
        ### exact computation
        kappas.append( (ma, mb, compute_kabc(lo, mo, la, ma, lb, mb, P=Po, k_hat=k_hat)) )
        ### assymptotic limit?
#        koab = gm.compute_kabc_assymptotic(lo, mo, la, ma, lb, -mo-ma, P=Po, k_hat=k_hat)
        ma += 1

      # compute N dependence once for each la, lb pair
      na = min_na
      while na <= max_na:
        ka_fact = na * (1+1./la)**0.5
        wa = alpha*la/na
        ya = cwo3a2 * ka_fact**2
        for ma, mb, koab in kappas:
          nb = int(math.ceil(max([1, (ka_fact - ko_fact)/lb_1, alpha*lb/max_w]))) # requirement that ko >= | ka - kb| and min_w <= wb <= max_w
          if (la == lb) and (ma == mb):
            nb = max(na, nb) # avoid duplicate pairs
          max_nb = math.floor(min( [(ka_fact + ko_fact)/lb_1, alpha*lb/min_w] )) # 
          while nb <= max_nb:
            wb = alpha*lb/nb
            yb = cwo3a2 * (nb*lb_1)**2
            _E = ms.compute_Ethr(-absO, wa, wb, ya, yb, koab)
            if _E <= Emax:
              if filename:
                print >>f, _E, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab
                f.flush()
              else:
                print      _E, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab

            nb += 1
        na += 1
      lb += 2
    la += 1

  if filename:
    f.close()
    return filename
  else:
    return True

##################################################
def compute_pairs_heuristic(parent_mode, O, min_l=1, max_l=100, min_w=0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, filename=False, Emax=np.infty):
  """
  compute a heuristic for daughter pairs. heuristic = (ya+yb)**2 + (O+wa+wb)**2. Daughter pairs are generated on the fly, with 
    min_l <= la <= max_l
    min_w <= wa <= max_w

  and all other parameters are limited by the following constraints
    |lo - la| <= lb <= lo + la
    lo + la + lb = even
    mo + ma + mb = 0
    ko >= |ka - kb|

  WARNING! mode parameters are computed internally without reference to compute_w() or compute_y()

  """
  no, lo, mo = parent_mode.get_nlm()
  Po = 2*np.pi/abs(parent_mode.w)
  cwo3a2 = c*wo**3*alpha**-2

  ko_fact = no*(1+1./lo)**0.5

  absO=abs(O)

  if filename:
    f=open(filename, "w")
    print >>f, absO
    print >>f, parent_mode.to_str_nlmwy()
  else:
    print absO
    print parent_mode.to_str_nlmwy()

  la = max(1, min_l) # instantiate la
  while la <= max_l:
    min_na = max(1, int(math.ceil(alpha*la/max_w))) # lower bound on na
    max_na = math.floor(alpha*la/min_w) # upper bound on na
    # demand that lb >= la (computational reasons)
    # demand that (la+lb+lo) is even
    if lo%2 == 0:
      lb = la
    else:
      lb = la + 1
    max_lb = la + lo # enforce triangle inequality
    while lb <= max_lb:
      lb_1 = (1+1./lb)**0.5

      # compute koab first
      kappas = []
      # require that (ma+mb+mo)=0
      ma = max(-la, -lb-mo) # bounds on ma and mb=-(mo+ma)
      if (la == lb): # avoid duplicate work
        if (mo == 0):
          max_ma = 0
        elif (mo > 0):
          max_ma = -1
        else: # mo < 0
          max_ma = min(la, lb-mo)
          ma = 1 # avoid duplicate pairs
      else:
        max_ma = min(la, lb-mo) # bounds on ma and mb=-(mo+ma)
      while ma <= max_ma:
        mb = -mo-ma
        ### exact computation
        kappas.append( (ma, mb, compute_kabc(lo, mo, la, ma, lb, mb, P=Po, k_hat=k_hat)) )
        ### assymptotic limit?
#        koab = gm.compute_kabc_assymptotic(lo, mo, la, ma, lb, -mo-ma, P=Po, k_hat=k_hat)
        ma += 1

      # compute N dependence once for each la, lb pair
      na = min_na
      while na <= max_na:
        ka_fact = na*(1+1./la)**0.5
        wa = alpha*la/na
        ya = cwo3a2 * ka_fact**2
        for ma, mb, koab in kappas:
          nb = int(math.ceil(max([1, (ka_fact - ko_fact)/lb_1, alpha*lb/max_w]))) # requirement that ko >= | ka - kb| and min_w <= wb <= max_w
          if (la == lb) and (ma == mb):
            nb = max(na, nb) # avoid duplicate pairs
          max_nb = math.floor(min( [(ka_fact + ko_fact)/lb_1, alpha*lb/min_w] )) # 
          while nb <= max_nb:
            wb = alpha*lb/nb
            yb = cwo3a2 * (nb*lb_1)**2
            _E = ms.compute_heuristic(-absO, wa, wb, ya, yb)
            if _E <= Emax:
              if filename:
                print >>f, _E, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab
                f.flush()
              else:
                print      _E, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab

            nb += 1
        na += 1

      lb += 2
    la += 1

  if filename:
    f.close()
    return filename
  else:
    return True

####################################################################################################
#
#
#                    local minima mode list generation methods
#             writes coupling lists to file (quickly), from which they can be read
#           formatting should be compatible with ggg_coupling_list_element.from_str()
####################################################################################################
def compute_min_pairs_Ethr(parent_mode, O, min_l=1, max_l=100, min_w=0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, filename=False, Emax=np.infty):
  """
  compute 3mode Ethr for daughter pairs. Daughter pairs are generated on the fly, with 
    min_l <= la <= max_l
    min_w <= wa <= max_w

  and all other parameters are limited by the following constraints
    |lo - la| <= lb <= lo + la
    lo + la + lb = even
    mo + ma + mb = 0
    ko >= |ka - kb|

  WARNING! mode parameters are computed internally without reference to compute_w() or compute_y()

  """
  no, lo, mo = parent_mode.get_nlm()
  Po = 2*np.pi/abs(parent_mode.w)
  cwo3a2 = c*wo**3*alpha**-2

  ko_fact = no*(1+1./lo)**0.5

  absO=abs(O)

  if filename:
    f=open(filename, "w")
    print >>f, absO
    print >>f, parent_mode.to_str_nlmwy()
  else:
    print absO
    print parent_mode.to_str_nlmwy()

  # we solve numerically for nb | na, la, lb
  eta = absO/alpha
  sigma = c*wo**3/alpha**3

  la = max(1, min_l) # instantiate la
  while la <= max_l:
    La = (1+1./la)
    la_1 = La**0.5
    min_na = max(1, int(math.ceil(alpha*la/max_w))) # lower bound on na
    max_na = math.floor(alpha*la/min_w) # upper bound on na

    # demand that lb >= la (computational reasons)
    # demand that (la+lb+lo) is even
    if lo%2 == 0:
      lb = la
    else:
      lb = la + 1
    max_lb = la + lo # enforce triangle inequality

    while lb <= max_lb:
      Lb = (1+1./lb)
      lb_1 = Lb**0.5
      bandwidth_min_nb = alpha*lb/max_w
      bandwidth_max_nb = alpha*lb/min_w

      # maximize koab over ma, mb

      # require that (ma+mb+mo)=0
      ma = max(-la, -lb-mo) # bounds on ma and mb=-(mo+ma)
      if (la == lb): # avoid duplicate work
        if (mo == 0):
          max_ma = 0
        elif (mo > 0):
          max_ma = -1
        else: # mo < 0
          max_ma = min(la, lb-mo)
          ma = 1 # avoid duplicate pairs
      else:
        max_ma = min(la, lb-mo) # bounds on ma and mb=-(mo+ma)

      best_ma = ma
      best_koab = 0.
      absbest_koab = 0.
      while ma <= max_ma:
        mb = -mo-ma
        ### exact computation
        koab = compute_kabc(lo, mo, la, ma, lb, mb, P=Po, k_hat=k_hat)
        if abs(koab) > absbest_koab:
          best_ma = ma
          best_koab = koab
          absbest_koab = abs(best_koab)
        ### assymptotic limit?
#        koab = gm.compute_kabc_assymptotic(lo, mo, la, ma, lb, -mo-ma, P=Po, k_hat=k_hat)
        ma += 1

      ma = best_ma
      mb = -mo-ma
      koab = best_koab

      # compute N dependence once for each la, lb pair
      na = min_na
      while na <= max_na:
        ka_fact = na*la_1

        min_nb = math.ceil(max([1, (ka_fact - ko_fact)/lb_1, bandwidth_min_nb])) # requirement that ko >= | ka - kb| and min_w <= wb <= max_w
        if (la == lb) and (ma == mb):
          min_nb = max(na, min_nb) # avoid duplicate pairs
        max_nb = math.floor(min( [(ka_fact + ko_fact)/lb_1, bandwidth_max_nb] ))

        if min_nb > max_nb: # avoid stupid work
          na += 1
          continue

        wa = gm.compute_w(na, la, alpha)
        ya = cwo3a2 * ka_fact**2 # faster than delegation?

        ### solve for nb(na,la,lb)
        lana_eta = 1.0*la/na-eta

        nb_roots = np.roots( [3*sigma**2*Lb**3 , 0 , 9*La*(sigma*Lb*na)**2, 0, Lb*((4+3*sigma**2)*lana_eta**2 + 9*(sigma*La*na**2)**2) ,    6*lana_eta*Lb*lb*(sigma**2+1) , 3*sigma**2*(na**6*La**3 + lana_eta*na**2*La + lb**2*Lb) + 2*lb**2*Lb, 2*lana_eta*lb*La*na**2*(6*sigma**2-1), La*na**2*lb*(3*sigma**2-2)] )
        nb_roots = nb_roots[nb_roots.imag == 0]
        nb_roots = nb_roots[(nb_roots >= min_nb)*(nb_roots <= max_nb)] # only keep positive definite roots

        for nb in nb_roots: ### iterate over all allowed roots and pick the best neighbors
          nb_f = int(math.floor(nb.real))
          wb_f = gm.compute_w(nb_f, lb, alpha)
          yb_f = cwo3a2 * nb_f**2 * Lb # faster than delegating?

          nb_c = int(math.ceil(nb.real))
          wb_c = gm.compute_w(nb_c, lb, alpha)
          yb_c = cwo3a2 * nb_c**2 * Lb

          E_f = ms.compute_Ethr(-absO, wa, wb_f, ya, yb_f, koab)
          E_c = ms.compute_Ethr(-absO, wa, wb_c, ya, yb_c, koab)
          if E_f <= E_c:
            E = E_f
            nb, wb, yb = nb_f, wb_f, yb_f
          else:
            E = E_c
            nb, wb, yb = nb_c, wb_c, yb_c

          if (E < Emax): # even if Emax=np.infty, we handle this correctly
            if filename:
              print >>f, E, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab
              f.flush()
            else:
              print E, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab

        na += 1
      lb += 2
    la += 1

  if filename:
    f.close()
    return filename
  else:
    return True


##################################################
def deprecated_compute_min_pairs_Ethr(parent_mode, O, min_l=1, max_l=100, min_w=0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, filename=False, Emax=np.infty):
  """
  compute 3mode Ethr for daughter pairs. Daughter pairs are generated on the fly, with 
    min_l <= la <= max_l
    min_w <= wa <= max_w

  and all other parameters are limited by the following constraints
    |lo - la| <= lb <= lo + la
    lo + la + lb = even
    mo + ma + mb = 0
    ko >= |ka - kb|

  WARNING! mode parameters are computed internally without reference to compute_w() or compute_y()

  """
  no, lo, mo = parent_mode.get_nlm()
  Po = 2*np.pi/abs(parent_mode.w)
  cwo3a2 = c*wo**3*alpha**-2

  ko_fact = no*(1+1./lo)**0.5

  absO=abs(O)

  if filename:
    f=open(filename, "w")
    print >>f, absO
    print >>f, parent_mode.to_str_nlmwy()
  else:
    print absO
    print parent_mode.to_str_nlmwy()

  # we solve numerically for nb | na, la, lb
  nb, _na, _la, _lb = sympy.symbols("nb na la lb")
  La = sympy.sqrt(_la*(_la+1))
  Lb = sympy.sqrt(_lb*(_lb+1))
  eta = absO/alpha
  sigma = c*wo**3/alpha**3

  func = (sigma**2 * La**2 * Lb**2 )/(_la **3 * _lb**3) * nb**3 * _na**3 * ( 1 + (eta - _la/_na -_la/nb)**2 / (sigma**2 * (_na**2 * La**2 / _la**2 + nb**2 * Lb**2/_lb**2)**2 ) )
  min_func = 3*sigma**2 * nb**2*_na**3 * ( (_na*La/_la)**2 + (nb*Lb/_lb)**2 )**3 + (3*nb**2*_na**3*(eta-_la/_na - _lb/nb)**2 + _na**3*nb**3*(eta - _la/_na - _lb/nb)*(2*_lb/nb**2) )*( (_na*La/_la)**2 + (nb*Lb/_lb)**2 ) - 2*_na**3*nb**3*(eta-_la/_na - _lb/nb)**2*(2*nb*(Lb/_lb)**2)

  la = max(1, min_l) # instantiate la
  while la <= max_l:
    func_la = func.subs(_la, la)
    min_func_la = min_func.subs(_la, la)

    min_na = max(1, int(math.ceil(alpha*la/max_w))) # lower bound on na
    max_na = math.floor(alpha*la/min_w) # upper bound on na
    print min_na, max_na

    # demand that lb >= la (computational reasons)
    # demand that (la+lb+lo) is even
    if lo%2 == 0:
      lb = la
    else:
      lb = la + 1
    max_lb = la + lo # enforce triangle inequality

    while lb <= max_lb:
      func_lalb = func_la.subs(_lb, lb)
      min_func_lalb = min_func_la.subs(_lb, lb)

      lb_1 = (1+1./lb)**0.5
      bandwidth_min_nb = alpha*lb/max_w
      bandwidth_max_nb = alpha*lb/min_w

      # maximize koab over ma, mb

      # require that (ma+mb+mo)=0
      ma = max(-la, -lb-mo) # bounds on ma and mb=-(mo+ma)
      if (la == lb): # avoid duplicate work
        if (mo == 0):
          max_ma = 0
        elif (mo > 0):
          max_ma = -1
        else: # mo < 0
          max_ma = min(la, lb-mo)
          ma = 1 # avoid duplicate pairs
      else:
        max_ma = min(la, lb-mo) # bounds on ma and mb=-(mo+ma)

      best_ma = ma
      best_koab = 0.
      absbest_koab = 0.
      while ma <= max_ma:
        mb = -mo-ma
        ### exact computation
        koab = compute_kabc(lo, mo, la, ma, lb, mb, P=Po, k_hat=k_hat)
        if abs(koab) > absbest_koab:
          best_ma = ma
          best_koab = koab
          absbest_koab = abs(best_koab)
        ### assymptotic limit?
#        koab = gm.compute_kabc_assymptotic(lo, mo, la, ma, lb, -mo-ma, P=Po, k_hat=k_hat)
        ma += 1

      ma = best_ma
      mb = -mo-ma
      koab = best_koab

      # compute N dependence once for each la, lb pair
      na = min_na
      while na <= max_na:
        ka_fact = na*(1+1./la)**0.5

        min_nb = math.ceil(max([1, (ka_fact - ko_fact)/lb_1, bandwidth_min_nb])) # requirement that ko >= | ka - kb| and min_w <= wb <= max_w
        if (la == lb) and (ma == mb):
          min_nb = max(na, min_nb) # avoid duplicate pairs
        max_nb = math.floor(min( [(ka_fact + ko_fact)/lb_1, bandwidth_max_nb] ))

        if min_nb > max_nb: # avoid stupid work
          na += 1
          continue

        func_lalbna = func_lalb.subs(_na, na)
        min_func_lalbna = min_func_lalb.subs(_na, na)

        wa = alpha*la/na
        ya = cwo3a2 * ka_fact**2

        ans = sympy.solve(min_func_lalbna, nb) # numerically solve for nb that are local minima
        print ans
        for _nb in sorted(ans):
          if sympify(_nb).is_real:
            if (min_nb <= _nb) and (_nb <= max_nb):
              floor_nb = math.floor(_nb)
              E_floor  = func_lalbna.subs(nb, floor_nb)
              ceil_nb = floor_nb + 1
              E_ceil  = func_lalbna.subs(nb, ceil_nb)
              if E_floor < E_ceil:
                this_E = E_floor/koab**2
                this_nb = int(floor_nb)
              else:
                this_E = E_ceil/koab**2
                this_nb = int(ceil_nb)

              if this_E <= Emax:
                if filename:
                  print >>f, this_E, "\t", (na, la, ma, wa, ya), "\t", (this_nb, lb, mb, alpha*lb/this_nb, cwo3a2*(this_nb*lb_1)**2), "\t", koab
                  f.flush()
                else:
                  print this_E, "\t", (na, la, ma, wa, ya), "\t", (this_nb, lb, mb, alpha*lb/this_nb, cwo3a2*(this_nb*lb_1)**2), "\t", koab

        na += 1
      lb += 2
    la += 1

  if filename:
    f.close()
    return filename
  else:
    return True

###################################################################################################
def compute_min_pairs_heuristic(parent_mode, O, min_l=1, max_l=100, min_w=0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, filename=False, Emax=np.infty):
  """
  compute 3mode heuristic for daughter pairs. Daughter pairs are generated on the fly, with 
    min_l <= la <= max_l
    min_w <= wa <= max_w

  and all other parameters are limited by the following constraints
    |lo - la| <= lb <= lo + la
    lo + la + lb = even
    mo + ma + mb = 0
    ko >= |ka - kb|

  WARNING! mode parameters are computed internally without reference to compute_w() or compute_y()

  """
  no, lo, mo = parent_mode.get_nlm()
  Po = 2*np.pi/abs(parent_mode.w)
  cwo3a2 = c*wo**3*alpha**-2

  ko_fact = no*(1+1./lo)**0.5

  absO=abs(O)

  if filename:
    f=open(filename, "w")
    print >>f, absO
    print >>f, parent_mode.to_str_nlmwy()
  else:
    print absO
    print parent_mode.to_str_nlmwy()

  # we solve numerically for nb | na, la, lb
  eta = absO/alpha
  sigma = c*wo**3/alpha**3

  la = max(1, min_l) # instantiate la
  while la <= max_l:
    La = (1+1./la)
    la_1 = La**0.5
    min_na = max(1, int(math.ceil(alpha*la/max_w))) # lower bound on na
    max_na = math.floor(alpha*la/min_w) # upper bound on na
    # demand that lb >= la (computational reasons)
    # demand that (la+lb+lo) is even
    if lo%2 == 0:
      lb = la
    else:
      lb = la + 1
    max_lb = la + lo # enforce triangle inequality

    while lb <= max_lb:
      Lb = (1+1./lb)
      lb_1 = (1+1./lb)**0.5
      bandwidth_min_nb = alpha*lb/max_w
      bandwidth_max_nb = alpha*lb/min_w

      # maximize koab over ma, mb

      # require that (ma+mb+mo)=0
      ma = max(-la, -lb-mo) # bounds on ma and mb=-(mo+ma)
      if (la == lb): # avoid duplicate work
        if (mo == 0):
          max_ma = 0
        elif (mo > 0):
          max_ma = -1
        else: # mo < 0
          max_ma = min(la, lb-mo)
          ma = 1 # avoid duplicate pairs
      else:
        max_ma = min(la, lb-mo) # bounds on ma and mb=-(mo+ma)

      best_ma = ma
      best_koab = 0.
      absbest_koab = 0.
      while ma <= max_ma:
        mb = -mo-ma
        ### exact computation
        koab = compute_kabc(lo, mo, la, ma, lb, mb, P=Po, k_hat=k_hat)
        if abs(koab) > absbest_koab:
          best_ma = ma
          best_koab = koab
          absbest_koab = abs(best_koab)
        ### assymptotic limit?
#        koab = gm.compute_kabc_assymptotic(lo, mo, la, ma, lb, -mo-ma, P=Po, k_hat=k_hat)
        ma += 1

      ma = best_ma
      mb = -mo-ma
      koab = best_koab

      # compute N dependence once for each la, lb pair
      na = min_na
      while na <= max_na:
        ka_fact = na*la_1

        min_nb = math.ceil(max([1, (ka_fact - ko_fact)/lb_1, bandwidth_min_nb])) # requirement that ko >= | ka - kb| and min_w <= wb <= max_w
        if (la == lb) and (ma == mb):
          min_nb = max(na, min_nb) # avoid duplicate pairs
        max_nb = math.floor(min( [(ka_fact + ko_fact)/lb_1, bandwidth_max_nb] ))

        if min_nb > max_nb: # avoid stupid work
          na += 1
          continue

        wa = gm.compute_w(na, la, alpha)
        ya = cwo3a2 * ka_fact**2 ### faster than delegation?

        ### solve for nb(na,la,lb)
        nb_roots = np.roots( [2*sigma**2*Lb**2 , 0 , 2*sigma**2*Lb*La*na**2 , 0, 0, lb*(eta-1.0*la/na), -lb**2] )
        nb_roots = nb_roots[(nb_roots >= min_nb)*(nb_roots <= max_nb)] # only keep positive definite roots

        for nb in nb_roots: ### iterate over all allowed roots and pick the best neighbors
          nb_f = int(math.floor(nb.real))
          wb_f = gm.compute_w(nb_f, lb, alpha)
          yb_f = cwo3a2 * nb_f**2 * Lb # faster than delegating?

          nb_c = int(math.ceil(nb.real))
          wb_c = gm.compute_w(nb_c, lb, alpha)
          yb_c = cwo3a2 * nb_c**2 * Lb

          h_f = ms.compute_heuristic(-absO, wa, wb_f, ya, yb_f)
          h_c = ms.compute_heuristic(-absO, wa, wb_c, ya, yb_c)
          if h_f <= h_c:
            h = h_f
            nb, wb, yb = nb_f, wb_f, yb_f
          else:
            h = h_c
            nb, wb, yb = nb_c, wb_c, yb_c

          if (h < Emax): # even if Emax=np.infty, we handle this correctly
            if filename:
              print >>f, h, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab
              f.flush()
            else:
              print h, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab

        na += 1
      lb += 2 # sum must remain even
    la += 1

  if filename:
    f.close()
    return filename
  else:
    return True

##################################################
def flawed_compute_min_pairs_heuristic(parent_mode, O, min_l=1, max_l=100, min_w=0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, filename=False, Emax=np.infty):
  """
  compute 3mode heuristic for daughter pairs. Daughter pairs are generated on the fly, with 
    min_l <= la <= max_l
    min_w <= wa <= max_w

  and all other parameters are limited by the following constraints
    |lo - la| <= lb <= lo + la
    lo + la + lb = even
    mo + ma + mb = 0
    ko >= |ka - kb|

  WARNING! mode parameters are computed internally without reference to compute_w() or compute_y()

  """
  no, lo, mo = parent_mode.get_nlm()
  Po = 2*np.pi/abs(parent_mode.w)
  cwo3a2 = c*wo**3*alpha**-2

  ko_fact = no*(1+1./lo)**0.5

  absO=abs(O)

  if filename:
    f=open(filename, "w")
    print >>f, absO
    print >>f, parent_mode.to_str_nlmwy()
  else:
    print absO
    print parent_mode.to_str_nlmwy()

  # we solve numerically for nb | na, la, lb
  eta = absO/alpha
  sigma = c*wo**3/alpha**3

  la = max(1, min_l) # instantiate la
  while la <= max_l:
    la_1 = (1+1./la)
    min_na = max(1, int(math.ceil(alpha*la/max_w))) # lower bound on na
    max_na = math.floor(alpha*la/min_w) # upper bound on na
    # demand that lb >= la (computational reasons)
    # demand that (la+lb+lo) is even
    if lo%2 == 0:
      lb = la
    else:
      lb = la + 1
    max_lb = la + lo # enforce triangle inequality

    while lb <= max_lb:
      lb_1 = (1+1./lb)
      min_nb = alpha*lb/max_w
      max_nb = alpha*lb/min_w

      # maximize koab over ma, mb

      # require that (ma+mb+mo)=0
      ma = max(-la, -lb-mo) # bounds on ma and mb=-(mo+ma)
      if (la == lb): # avoid duplicate work
        if (mo == 0):
          max_ma = 0
        elif (mo > 0):
          max_ma = -1
        else: # mo < 0
          max_ma = min(la, lb-mo)
          ma = 1 # avoid duplicate pairs
      else:
        max_ma = min(la, lb-mo) # bounds on ma and mb=-(mo+ma)

      best_ma = ma
      best_koab = 0.
      absbest_koab = 0.
      while ma <= max_ma:
        mb = -mo-ma
        ### exact computation
        koab = compute_kabc(lo, mo, la, ma, lb, mb, P=Po, k_hat=k_hat)
        if abs(koab) > absbest_koab:
          best_ma = ma
          best_koab = koab
          absbest_koab = abs(best_koab)
        ### assymptotic limit?
#        koab = gm.compute_kabc_assymptotic(lo, mo, la, ma, lb, -mo-ma, P=Po, k_hat=k_hat)
        ma += 1

      ma = best_ma
      mb = -mo-ma
      koab = best_koab

      ### compute best na, nb allowed
      na_nb = ((la*lb_1)/(lb*la_1))**(1.0/3) ### ratio of na/nb
      nb_roots = np.roots( [2*sigma**2 * (na_nb**2 *la_1 + lb_1) * lb_1, 0 , 0 , 0, 0, -eta*lb, la*lb*na_nb + lb**2] ) ### 6th order polynomial for nb
      nb_roots = nb_roots[nb_roots.imag == 0.0] # keep only real roots
      print nb_roots
      print nb_roots*na_nb
      nb_roots = nb_roots[(nb_roots >= min_nb)*(nb_roots <= max_nb)] # only keep positive definite roots
      na_roots = nb_roots*na_nb
      print nb_roots
      print na_roots
      logical = (na_roots >= min_na)*(na_roots <= max_na)
      na_roots = na_roots[logical]
      nb_roots = nb_roots[logical]

      ### for each na,nb pair, select the best neighbor with na, nb integers
      for na, nb in zip(na_roots, nb_roots):
        if abs(na*la_1**0.5 - nb*lb_1**0.5) > ko_fact: # skip because coupling goes to zero
          continue 

        ### iterate over all possible combos
        h = np.infty
        for _na in [int(math.floor(na.real)), int(math.ceil(na.real))]:
          wa = gm.compute_w(_na, la, alpha)
          ya = gm.compute_y(_na, la, c, wo, alpha)
          for _nb in [int(math.floor(nb.real)), int(math.ceil(nb.real))]:
            wb = gm.compute_w(_nb, lb, alpha)
            yb = gm.compute_y(_nb, lb, c, wo, alpha)

            _h = ms.compute_heuristic(-absO, wa, wb, ya, yb)
            if _h < h:
              h = _h
              this_nwy = (_na, wa, ya, _nb, wb, yb)

        if (h < Emax): # even if Emax=np.infty, we handle this correctly
          na, wa, ya, nb, wb, yb = this_nwy
          if filename:
            print >>f, h, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab
            f.flush()
          else:
            print h, "\t", (na, la, ma, wa, ya), "\t", (nb, lb, mb, wb, yb), "\t", koab

      lb += 2 # sum must remain even
    la += 1

  if filename:
    f.close()
    return filename
  else:
    return True

###################################################################################################
def depredated_compute_min_pairs_heuristic(parent_mode, O, min_l=1, max_l=100, min_w=0, max_w=1.0, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, filename=False, Emax=np.infty):
  """
  compute 3mode heuristic for daughter pairs. Daughter pairs are generated on the fly, with 
    min_l <= la <= max_l
    min_w <= wa <= max_w

  and all other parameters are limited by the following constraints
    |lo - la| <= lb <= lo + la
    lo + la + lb = even
    mo + ma + mb = 0
    ko >= |ka - kb|

  WARNING! mode parameters are computed internally without reference to compute_w() or compute_y()

  """
  no, lo, mo = parent_mode.get_nlm()
  Po = 2*np.pi/abs(parent_mode.w)
  cwo3a2 = c*wo**3*alpha**-2

  ko_fact = no*(1+1./lo)**0.5

  absO=abs(O)

  if filename:
    f=open(filename, "w")
    print >>f, absO
    print >>f, parent_mode.to_str_nlmwy()
  else:
    print absO
    print parent_mode.to_str_nlmwy()

  # we solve numerically for nb | na, la, lb
  nb, _na, _la, _lb = sympy.symbols("nb na la lb")
  La = sympy.sqrt(_la*(_la+1))
  Lb = sympy.sqrt(_lb*(_lb+1))
  eta = absO/alpha
  sigma = c*wo**3/alpha**3

  func = sigma**2 * ( (_na*La/_la)**2 + (nb*Lb/_lb)**2 )**2 + (eta - _la/_na - _lb/nb)**2
  min_func = 2*sigma**2 *nb**4 *_na * Lb**2 * ( (_na*La/_la)**2 + (nb*Lb/_lb)**2 ) + (eta*_na*nb - _la*nb - _lb*_na)*_lb**3

  la = max(1, min_l) # instantiate la
  while la <= max_l:
    func_la = func.subs(_la, la)
    min_func_la = min_func.subs(_la, la)

    min_na = max(1, int(math.ceil(alpha*la/max_w))) # lower bound on na
    max_na = math.floor(alpha*la/min_w) # upper bound on na
    # demand that lb >= la (computational reasons)
    # demand that (la+lb+lo) is even
    if lo%2 == 0:
      lb = la
    else:
      lb = la + 1
    max_lb = la + lo # enforce triangle inequality

    while lb <= max_lb:
      func_lalb = func_la.subs(_lb, lb)
      min_func_lalb = min_func_la.subs(_lb, lb)

      lb_1 = (1+1./lb)**0.5
      bandwidth_min_nb = alpha*lb/max_w
      bandwidth_max_nb = alpha*lb/min_w

      # maximize koab over ma, mb

      # require that (ma+mb+mo)=0
      ma = max(-la, -lb-mo) # bounds on ma and mb=-(mo+ma)
      if (la == lb): # avoid duplicate work
        if (mo == 0):
          max_ma = 0
        elif (mo > 0):
          max_ma = -1
        else: # mo < 0
          max_ma = min(la, lb-mo)
          ma = 1 # avoid duplicate pairs
      else:
        max_ma = min(la, lb-mo) # bounds on ma and mb=-(mo+ma)

      best_ma = ma
      best_koab = 0.
      absbest_koab = 0.
      while ma <= max_ma:
        mb = -mo-ma
        ### exact computation
        koab = compute_kabc(lo, mo, la, ma, lb, mb, P=Po, k_hat=k_hat)
        if abs(koab) > absbest_koab:
          best_ma = ma
          best_koab = koab
          absbest_koab = abs(best_koab)
        ### assymptotic limit?
#        koab = gm.compute_kabc_assymptotic(lo, mo, la, ma, lb, -mo-ma, P=Po, k_hat=k_hat)
        ma += 1

      ma = best_ma
      mb = -mo-ma
      koab = best_koab

      # compute N dependence once for each la, lb pair
      na = min_na
      while na <= max_na:
        ka_fact = na*(1+1./la)**0.5

        min_nb = math.ceil(max([1, (ka_fact - ko_fact)/lb_1, bandwidth_min_nb])) # requirement that ko >= | ka - kb| and min_w <= wb <= max_w
        if (la == lb) and (ma == mb):
          min_nb = max(na, min_nb) # avoid duplicate pairs
        max_nb = math.floor(min( [(ka_fact + ko_fact)/lb_1, bandwidth_max_nb] ))

        if min_nb > max_nb: # avoid stupid work
          na += 1
          continue

        func_lalbna = func_lalb.subs(_na, na)
        min_func_lalbna = min_func_lalb.subs(_na, na)

        wa = alpha*la/na
        ya = cwo3a2 * ka_fact**2

        ans = sympy.solve(min_func_lalbna, nb) # numerically solve for nb that are local minima
        for _nb in sorted(ans):
          if sympify(_nb).is_real:
            if (min_nb <= _nb) and (_nb <= max_nb):
              floor_nb = math.floor(_nb)
              heuristic_floor  = func_lalbna.subs(nb, floor_nb)
              ceil_nb = floor_nb + 1
              heuristic_ceil  = func_lalbna.subs(nb, ceil_nb)
              if heuristic_floor < heuristic_ceil:
                this_heuristic = heuristic_floor
                this_nb = int(floor_nb)
              else:
                this_heuristic = heuristic_ceil
                this_nb = int(ceil_nb)

              if this_heuristic <= Emax:
                if filename:
                  print >>f, this_heuristic, "\t", (na, la, ma, wa, ya), "\t", (this_nb, lb, mb, alpha*lb/this_nb, cwo3a2*(this_nb*lb_1)**2), "\t", koab
                  f.flush()
                else:
                  print this_heuristic, "\t", (na, la, ma, wa, ya), "\t", (this_nb, lb, mb, alpha*lb/this_nb, cwo3a2*(this_nb*lb_1)**2), "\t", koab

        na += 1
      lb += 2
    la += 1

  if filename:
    f.close()
    return filename
  else:
    return True

####################################################################################################
#
#
#                 high-level methods to define the scope of mode list building, etc.
#
#
####################################################################################################
def compute_pairs(metric, parent_mode, O, catalog_dir, min_l=1, max_l=100, min_w=0, max_w=1e-4, alpha=4e-3, c=2e-11, wo=1e-5, k_hat=5e4, Emax=np.infty, maxp=1, verbose=True, max_num_pairs=-1):
  """
  defines the scope of jobs needed to build a mode list. Then submits those jobs, waits, and returns the resulting mode list.

  metric : the metric used to rank modes. This identifies the catalog with the function used to generate it [Ethr, min_Ethr, heuristic, min_heuristic]
   ==> func : the function with which we build the list [compute_paris_Ethr, compute_pairs_heuristic, compute_min_pairs_Ethr, compute_min_pairs_heuristic]
  no, lo, mo : mode numbers of parent mode
  O : driving frequency
  catalog_dir : path to existing catalog files directory, used to determine the minimum amount of work required.
  min_l : min_l for daughter mode A
  max_l : max_l for daughter mode A
  min_w : min frequency for daughter mode A
    because of the peculiarities derived from comparing min_w to the value from a filename, we only allow for certain (rounded) precision here. The rounding is done before any jobs are launched. It allow for precision up to 1e-11 Hz
  max_w : max frequency for daughter mode A
    because of the peculiarities derived from comparing max_w to the value from a filename, we only allow for certain (rounded) precision here. The rounding is done before any jobs are launched. It allow for precision up to 1e-11 Hz
  alpha, c, wo : parameters of the primary star used to compute g-mode parameters
  k_hat : prefactor for coupling constant (k_oab)
  Emax : the maximum value of the ranking statistic recorded
  
  We split the current job into sub-jobs that have (max_l=min_l). For each sub-job, we then define the frequency ranges that need to be computed. These jobs are all submitted in parallel in batches of size "maxp"

  returns the list of all useful filenames (including 

  """
  catalog_filenames = sorted(glob.glob(catalog_dir+"/*ctg"))
  if verbose:
#    print "searching for existing .ctg files in %s\nfound:" % catalog_dir
#    for catalog_filename in catalog_filenames:
#      print "\t%s" % catalog_filename
    print "searching for existing .ctg files in %s\nfound %d files:" % (catalog_dir, len(catalog_filenames))

  no, lo, mo = parent_mode.get_nlm()
  absO = abs(O)

  # make sure min_w and max_w are compatible with string formats
  min_w = float( "%fe-5" % (min_w*1e5) )
  max_w = float( "%fe-5" % (max_w*1e5) )

  ### determine what new work is required (for each l)
  good_filenames = []
  new_filenames = []
  new_jobs = []
  for l in range(min_l, max_l+1):
    bandwidths = [(min_w, max_w)]
    for filename in catalog_filenames: # check all existing files
      if check_filename(filename, metric, no, lo, mo, absO, min_l=l, max_l=l, min_w=min_w, max_w=max_w, verbose=False):
        good_filenames.append( filename ) # add this to list of useful files
        _metric, _no, _lo, _mo, _absO, _min_l, _max_l, _min_w, _max_w = parse_filename(filename)
        new_bandwidths = [] # new list of bandwidths
        for w_min, w_max in bandwidths: # determine if new work is still required
          if (_min_w <= w_min) and (w_min < _max_w): # check left endpoint
            w_min = _max_w
          if (_min_w < w_max) and (w_max <= _max_w): # check right endpoint
            w_max = _min_w
          if w_max > w_min:
            new_bandwidths.append( (w_min, w_max) ) # if this is a non-trivial range, we add it to the new bandwidths
        bandwidths = new_bandwidths[:] # update
    ### define new jobs by the argument strings  
    for w_min, w_max in bandwidths:
      new_filename = catalog_dir+"/"+generate_filename(metric, no, lo, mo, absO, l, l, w_min, w_max)
      new_filenames.append( new_filename )
      new_jobs.append( (parent_mode, O, l, l, w_min, w_max, alpha, c, wo, k_hat, new_filename, Emax) )

  if verbose:
    print "existing (useful) catalog files:"
    for catalog_filename in good_filenames:
      print "\t%s" % catalog_filename
    print "generating new catalog files:"
    for catalog_filename in new_filenames:
      print "\t%s" % catalog_filename
    import time
    tos = []


  ### submit new work
  procs = []
  # sort jobs by their expected length (1+Lmax-Lmin)*(wmax-wmin)
  new_jobs.sort(key=lambda job: (1+job[3]-job[2])*(job[5]-job[4]) )
  for args in new_jobs:
    if verbose: 
      print "LAUNCHING\n\tmetric=%s\n\targs=" % (metric), args[0].get_nlmwy(), args[1:]
      tos.append(time.time())

    p = mp.Process(target=metricD[metric], args=args)
    p.start()
    procs.append((p, args))

    while len(procs) >= maxp: # wait
      for ind, (p, _) in enumerate(procs):
        if not p.is_alive():
          break
      else:
        continue
      p, p_args = procs.pop(ind)
      if verbose:
        print "RECEIVING from\n\tmetric=%s\n\targs=" % (metric), p_args[0].get_nlmwy(), p_args[1:]
        print "done: %f seconds" % (time.time()-tos.pop(ind))

      if p.exitcode < 0:
        sys.exit("\njob Failed! with status "+p.exitcode() )

  while len(procs):
    for ind, (p, _) in enumerate(procs):
      if not p.is_alive():
        break
    else:
      continue
    p, p_args = procs.pop(ind)
    if verbose:
      print "RECEIVING from\n\tmetric=%s\n\targs=" % (metric), p_args[0].get_nlmwy(), p_args[1:]
      print "done: %f seconds" % (time.time()-tos.pop(ind))

    if p.exitcode < 0:
      sys.exit("\njob Failed! with status "+p.exitcode() )

#    if len(procs) >= maxp: # wait
#      p, p_args = procs.pop(0)
#      if verbose: print "WAITING for\n\tmetric=%s\n\targs=" % (metric), p_args[0].get_nlmwy(), p_args[1:]
#
#      p.join()
#      if p.exitcode < 0:
#        sys.exit("\njob Failed!\n\tmetric=%s\n\targs=" % (metric), p_args[0].get_nlmwy(), p_args[1:], "\nexited with status "+p.exitcode() )
#
#      if verbose: 
#        print "\tdone: %f seconds" % (time.time()-tos.pop(0))
#
#  for p, p_args in procs: # wait for the rest of the jobs
#    if verbose: print "WAITING for\n\tmetric=%s\n\targs=" % (metric), p_args[0].get_nlmwy(), p_args[1:]
#
#    p.join()
#    if p.exitcode < 0:
#      sys.exit("\njob Failed!\n\tmetric=%s\n\targs=" % (metric), p_args[0].get_nlmwy(), p_args[1:], "\nexited with status "+p.exitcode() )
#
#    if verbose: print "\tdone: %f seconds" % (time.time()-tos.pop(0))

  if (max_num_pairs >= 0) and ("min" not in metric): # max_num_pairs is a non-negative number ==> we attempt to downsample the new_filename
    from clean_catalogs import clean
    for new_filename in new_filenames:
      clean(new_filename, new_filename, max_num_pairs, alpha, c, wo, k_hat, coupling_type="ggg", verbose=verbose)

  return new_filenames+good_filenames

####################################################################################################
#
#
#            look-up dictionaries for the entire module
#                put down here so everything is defined
#
####################################################################################################
### define a dictionary that relates metrics to the appropriate functions
global metricD
metricD = {"Ethr":compute_pairs_Ethr, "min_Ethr":compute_min_pairs_Ethr, "heuristic":compute_pairs_heuristic, "min_heuristic":compute_min_pairs_heuristic}

