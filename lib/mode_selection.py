usage="""written to provide common functions for all mode selection and coupling lists"""

import os, shutil
import numpy as np
import nmode_utils as nm_u
import copy

####################################################################################################
#
#
#                            utility functions
#
#
####################################################################################################
def compute_Ethr(O, w1, w2, y1, y2, k):
  """ 
  compute the linear threshold energy for two modes. 
  use this function for self-coupled modes, with w1=w2 and y1=y2
  """
  if k == 0:
    return np.Inf
  else:
    return (y1*y2)/(4.0*k**2*w1*w2) * ( 1 + (O + w1 + w2)**2/float(y1 + y2)**2 ) # Eqn (A14) from Nevin's paper

##################################################
def compute_heuristic(O, w1, w2, y1, y2):
  """
  a heuristic for the 3mode selection criterion
  """
  return (y1+y2)**2 + (O+w1+w2)**2

##################################################
def compute_collE( sorted_Ethrs ):
  """
  computes the collective instability threshold for this set of Ethrs.
    Ethrs is the list of 3mode Ethr for all coupling in which this mode participates as a daughter sorted in order of increasing Ethr
  """
  if not len(sorted_Ethrs):
    return np.infty
#  return sorted_Ethrs[-1]
#  return 1.0*sorted_Ethrs[-1]/len(sorted_Ethrs)**2
  return np.min( np.array(sorted_Ethrs) / np.arange(1,len(sorted_Ethrs)+1)**2 )


##################################################
def compute_Elin(O, wo, yo, U):
  return (wo*U)**2 / ((O-wo)**2 + yo**2)

##################################################
def compute_wo(Mprim, Rprim):
  """ wo = (G*(Mprim*Msun)/(Rprim*Rusn)**3)**0.5 """
  return (nm_u.G*(Mprim*nm_u.Msun)/(Rprim*nm_u.Rsun)**3)**0.5

##################################################
def compute_linearized_growth_rate(w1, w2, y1, y2, O, k, Ao):
  """
  y+ = y1+y2
  y- = y1-y2
  D  = O+w1+w2

  X = (y+)**2 - (D)**2 + 16*wa*wb*k**2*Ao**2
  Y = 2*(y-)*(D))

  Re{s} = -(y+) + ( X**2 + Y**2 )**(1/4) * cos( (1/2)* atan( Y/X ) )
        = -(y+) + ( 2*X**2 + Y**2 + 2*X*( X**2 + Y**2 )**(1/2) )**(1/4)

  This is the expected growth rate of unstable 3mode pairs (from linearization around the linear equilib)
  """
  X = (y1+y2)**2 - (O+w1+w2)**2 +16*w1*w2*(k*Ao)**2
  Y = 2*(y1-y2)*(O+w1+w2)

  return -(y1+y2) + ( 2*X**2 + Y**2 + 2*X*( X**2 + y**2 )**(1./2) )**(1./4)

####################################################################################################
#
#
#         basic functions to read mode lists to and from ASCii files
#
#
####################################################################################################
#########################
def insert_in_place(item, mylist):
  """
  puts the item in place in the list. 
  item must have the form: (rank, stuff)
  mylist must have the form: [item1, item2, ...]
    with itemj = (rankj, stuffj)

  puts things in order of increasing rank
  """
  irank, stuff = item
  for ind, (rank, mystuff) in enumerate(mylist):
    if irank < rank:
      mylist.insert(ind, [irank, stuff] )
      break
  else:
    mylist.append( [irank, stuff] )

  return mylist

#########################
def insert_sortedList_in_place(sorted_list, mylist, sorted=True):
  """
  puts the items in sorted list in place in mylist.
  sorted_list must have the form [(rank, mystuff), ...]
  mylist must have the form [(rank, mystuff), ...]

  puts things in order of increasing rank
  if sorted == False, we first sort the "sorted list"
  """
  if not sorted:
    sorted_list.sort(key=lambda l: l[0])

  ml_ind = 0 # index of mylist
  while len(sorted_list):
    sl_rank, stuff = sorted_list.pop(0) # get the first element
    while ml_ind < len(mylist):
      rank, _ = mylist[ml_ind]
      if sl_rank < rank: # insert at this index
        break
      else:
        ml_ind += 1 
    mylist.insert(ml_ind, [sl_rank, stuff] )
    ml_ind += 1 # increment because we added something to mylist

  return mylist

#########################
def check_mode(mode, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False, verbose=False):
  """
  check a single mode and ensure that it obeys the selection rules defined as input arguments.

  WARNING! we check abs(w), not w
  """
  n, l, m, w, y = mode.get_nlmwy()
  absw = abs(w)
  if min_n and (n < min_n):
    if verbose: print "min_n fails"
    return False
  if max_n and (n > max_n):
    if verbose: print "max_n fails"
    return False
  if min_l and (l < min_l):
    if verbose: print "min_l fails"
    return False
  if max_l and (l > max_l):
    if verbose: print "max_l fails"
    return False
  if min_w and (absw < min_w):
    if verbose: print "min_w fails"
    return False
  if max_w and (absw > max_w):
    if verbose: print "max_w fails"
    return False

  return True

####################################################################################################
#
#
#                   coupling_list class
#
#
####################################################################################################
class coupling_list:
  """
  this class handles coupling lists, reading and writing them to/from disk, etc.

  contains basic functionality, like sorting and combining different lists

  data contains:
    filenames : list of filenames associated with this list
    parent_mode : an instance of networks.mode()
    O : forcing freq against which the daughters couple
    num_pairs : the maximum number of pairings read from each file in filenames
      we should not trust this list if we request more couplings than num_pairs
  """

  def __init__(self, parent_mode=None):
    self.filenames = []
    self.couplings = []
    self.O = None
    self.parent_mode = parent_mode
    self.num_pairs = -1 # an unphysical number so we don't get confused.

  ###
  def to_unique_couplings(self):
    """
    takes the list of couplings to a unique list, in which each coupling appears only once
    returns None

    if couplings are repeated, it keeps the one that appears first in the list (lowest metric_value)
    """
    print "to_unique_couplings"
    new_couplings = []
    pairs = []

    for coupling in self.couplings:
      nlms1 = coupling.dmode1.get_nlms()
      nlms2 = coupling.dmode2.get_nlms()
      if nlms1 > nlms2:
        this_pair = (nlms1, nlms2)
      else:
        this_pair = (nlms2, nlms1)
      for ind, pair in enumerate(pairs):
        if pair == this_pair:
          break # pair already exists
        elif pair > this_pair: # pair does not exist, but would have by now
          pairs.insert(ind, this_pair) # add it to pairs
          new_couplings.append( coupling )
          break
      else:
        pairs.append( this_pair ) # reached end of list without finding or adding coupling, so we append
        new_couplings.append( coupling )

    self.couplings = new_couplings
    return self

  ###
  def load_mode_lists(self, filenames, num_pairs=-1, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """
    loads data from sorted files into a sorted mode list. Looks for the first num_pairs modes in each list and loads them. 

    the internal list is kept sorted.

    a subset of modes can be specified with the optional arguments
    """

    # check for basic compatibility
    if self.num_pairs == -1:
      self.num_pairs = num_pairs
    elif self.num_pairs < num_pairs:
      sys.exit("already loaded lists with num_pairs=%d. Cannot load more lists with num_pairs=%d>%d" %(self.num_pairs, num_pairs, self.num_pairs))

    # set up iteration
    if isinstance(filenames, str):
      filenames = [filenames]
    couplings = [ (c.metric_value, c) for c in self.couplings]

    for filename in filenames:

      self.filenames.append(filename)

      f = open(filename, 'r')
      absO = float(f.readline().strip())
      Mo = networks.mode().from_str_nlmwy(f.readline().strip())

      # check for forcing freq compatibility
      if self.O == None:
        self.O = absO
      elif self.O != absO:
        sys.exit("incompatible forcing frequencies!")
      # check for parent mode compatibility
      if self.parent_mode == None:
        self.parent_mode = Mo
      elif (self.parent_mode.nlm()[:-1] != Mo.get_nlm()[:-1]) or ((self.parent_mode.m) != abs(Mo.m)) or (self.parent_mode.get_wyU()[:-1] != Mo.get_wyU()[:-1]):
        sys.exit("incompatible parent modes!\n\texisting: %s\n\tnew: %s" % (str(self.parent_mode.get_nlmwy()), str(Mo.get_nlmwy())) )

      # read in modes
      ind == 0
      _couplings = []
      for line in f:
        if (num_pair != -1) and (ind > num_pair): # only load num_pair couplings
          break

        if line[0] != "#":
          coupling = coupling_list_element().from_string( line.strip() ).setup( self.parent_mode )
          if check_mode(coupling.dmode1, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w) and check_mode(coupling.dmode2, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w):
            _couplings = insert_in_place( (coupling.metric_value, coupling), _couplings )

      f.close()
      couplings = ms.insert_sortedList_in_place( _couplings, couplings, sorted=True )

    self.couplings = [c for metric_value, c in couplings]

    return self

  ###
  def load_unsorted_mode_lists(self, filenames, num_pairs=-1, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """
    loads data from unsorted files into a sorted mode list. Looks for the first num_pairs modes in each list and loads them. 

    the internal list is kept sorted.

    a subset of modes can be specified with the optional arguments
    """

    # check for basic compatibility
    if self.num_pairs == -1:
      self.num_pairs = num_pairs
    elif self.num_pairs < num_pairs:
      sys.exit("already loaded lists with num_pairs=%d. Cannot load more lists with num_pairs=%d>%d" %(self.num_pairs, num_pairs, self.num_pairs))

    # set up iteration
    if isinstance(filenames, str):
      filenames = [filenames]
    couplings = [ (c.metric_value, c) for c in self.couplings]

    for filename in filenames:

      self.filenames.append(filename)

      f = open(filename, 'r')
      absO = float(f.readline().strip())
      Mo = networks.mode().from_str_nlmwy(f.readline().strip())

      # check for forcing freq compatibility
      if self.O == None:
        self.O = absO
      elif self.O != absO:
        sys.exit("incompatible forcing frequencies!")
      # check for parent mode compatibility
      if self.parent_mode == None:
        self.parent_mode = Mo
      elif (self.parent_mode.nlm()[:-1] != Mo.get_nlm()[:-1]) or ((self.parent_mode.m) != abs(Mo.m)) or (self.parent_mode.get_wyU()[:-1] != Mo.get_wyU()[:-1]):
        sys.exit("incompatible parent modes!\n\texisting: %s\n\tnew: %s" % (str(self.parent_mode.get_nlmwy()), str(Mo.get_nlmwy())))

      # read in modes
      _couplings = []
      for line in f:
        if line[0] != "#":
          coupling = coupling_list_element().from_string( line.strip() ).setup( self.parent_mode )
          if check_mode(coupling.dmode1, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w) and check_mode(coupling.dmode2, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w):
            if num_pairs == -1: # load everything 
              _couplings = insert_in_place( (coupling.metric_value, coupling), couplings )
            else:
              if len(_couplings) < num_pairs: # add the mode
                _couplings = insert_in_place( (coupling.metric_value, coupling), _couplings)
              else:
                if coupling.metric_value < _couplings[-1][0]:
                  _couplings.pop()
                  _couplings = insert_in_place( (coupling_metric_value, coupling), _couplings)

      f.close()
#      for c in _couplings:
#        couplings = insert_in_place( c, couplings )
      couplings = ms.insert_sortedList_in_place( _couplings, couplings, sorted=False )

    self.couplings = [c for metric_value, c in couplings]

    return self

  ###
  def combine_with(self, clist):
    """
    combines the two coupling_lists
    checks to make sure they are compatible (parent_mode and absO match exactly)

    if they are compatible, returns a NEW coupling_list object
      things like parent_mode, O are taken from self

    WARNING: if you combine partially redundant lists, you may get two identical entries.
    """
    if (self.parent_mode.get_nlmwy() != clist.parent_mode.get_nlmwy()) or (self.O != clist.O): # we're good
      sys.exit("incompatible mode lists")

    combined_coupling_list = coupling_list()
    combined_coupling_list.filenames = list(set(self.filenames[:] + clist.filenames[:]))
    combined_coupling_list.O = self.O
    combined_coupling_list.parent_mode = copy.deepcopy( self.parent_mode )
    if (self.num_pairs == -1) and (clist.num_pairs == -1):
      combined_coupling_list.num_pairs = -1
    elif (self.num_pairs == -1):
      combined_coupling_list.num_pairs = clist.num_pairs
    elif (clist.num_pairs == -1):
      combined_coupling_list.num_pairs = self.num_pairs
    else:
      combined_coupling_list.num_pairs = min(self.num_pairs, clist.num_pairs)

    ### iterate over lists and add in sorted order
    couplings = [ (c.metric_value, copy.deepcopy(c)) for c in self.couplings ]
    for coupling in clist.couplings:
      coupling = copy.deepcopy( coupling ) 
      couplings = insert_in_place( (coupling.metric_value, coupling), couplings ) 

    combined_coupling_list.couplings = [c for metric_value, c in couplings]

    return combined_coupling_list

  ###
  def write_mode_list(self, filename):
    """
    dumps the expected coupling_list structure to ASCii file
    """
    f = open(filename, "w")
    print >>f, abs(self.O)
    print >>f, self.parent_mode.to_str_nlmwy()
    for coupling in self.couplings:
      print >>f, coupling.to_str()
    f.close()
    return filename

  ###
  def within_bounds(self, min_n=False, max_n=False, min_l=False, max_l=False, min_w=False, max_w=False):
    """
    returns a new mode list containing modes that only live within a certain parameter space
    """
    new_coupling_list = coupling_list()
    new_coupling_list.filenames = copy.deepcopy(self.filenames)
    new_coupling_list.O = copy.deepcopy(couplings.O)
    new_coupling_list.parent_mode = copy.deepcopy( self.parent_mode )
    new_coupling_list.num_pairs = copy.deepcopy( self.num_pairs )

    new_coupling_list.couplings = []

    for coupling in self.couplings:
      if check_mode(coupling.dmode1, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w) and check_mode(coupling.dmode2, min_n=min_n, max_n=max_n, min_l=min_l, max_l=max_l, min_w=min_w, max_w=max_w):
        new_coupling_list.couplings.append( copy.deepcopy(coupling) )

    return new_coupling_list

  ###
  def to_triples(self, num_pairs):
    """
    converts the mode list to a triples format compatible with networks.network.add_couplings()
    """
    if (self.num_pairs != -1) and (num_pairs > self.num_pairs):
      sys.exit("cannot report more pairs than the maximum read from any given file")

    triples = []
    for coupling in couplings[:num_pairs]:
      triples.append( (parent_mode, coupling.dmode1, coupling.dmode2, coupling.k) )  
    return triples

##################################################
class coupling_list_element:
  """
  a class representing the coupling_list elements
  data contains:
    metric_value : used to sort
    dmode1 : first daughter mode. must be an instance of networks.mode()
    dmode2 : second daughter mode. must be an instance of networks.mode()
    k : coupling coefficient

  coupling_list_element does not know about the parent mode. That's stored by the coupling_list as a whole.
  """

  def __init__(self):
    self.metric_value = np.infty
    self.dmode1 = None
    self.dmode2 = None
    self.k = None

  ###
  def to_str(self):
    """
    converts to standard string form for mode list ASCii files
    """
    x, expx = nm_u.float_to_scientific(self.metric_value)
    k, expk = nm_u.float_to_scientific(self.k)
    return "%.14fe%d\t%s\t%s\t%.14fe%d" % (x, expx, self.dmode1.to_str_nlmwy(tuple=True), self.dmode2.to_str_nlmwy(tuple=True), k, expk) 

  ###
  def from_str(self, string):
    """
    parse dpair from catalog ASCII files
    """
    metric_value, dmode1str, dmode2str, k = string.split()
    self.metric_value = float(metric_value)
    self.k = float(self.k)

    n, l, m, w, y = dmode1str.strip("(").strip(")").split(",")
    self.dmode1 = networks.mode(int(n), int(l), int(m), float(w), float(y))

    n, l, m, w, y = dmode2str.strip("(").strip(")").split(",")
    self.dmode2 = networks.mode(int(n), int(l), int(m), float(w), float(y))

    line = string.strip().split("(")
    line = [item.strip().strip(")") for item in line[:-1] + line[-1].split(")")]

    return self

  ###
  def set_sign(self, parent_sign):
    """
    daughter modes must have opposite sign of the parent
    """
    self.dmode1.w = (-parent_sign)*abs(self.dmode1.w)
    self.dmode2.w = (-parent_sign)*abs(self.dmode2.w)
    return self

  ###
  def set_m(self, parent_m):
    """
    require that parent_m + dmode1.m + dmode2.m = 0
    """
    if (self.dmode1.m + self.dmode2.m) == parent_m:
      self.dmode1.m = -1*self.dmode1.m
      self.dmode2.m = -1*self.dmode2.m
      return self
    if (self.dmode1.m + self.dmode2.m) == -parent_m:
      return self
    else:
      sys.exit("incompatible daughter modes and parent mode!")

  ### 
  def setup(self, parent_mode):
    """
    sets up a mode so that it is compatible with parent_mode.
    also leaves the default mode forcing as U=[]
    """
    self.set_sign( parent_mode.w/abs(parent_mode.w) )
    self.set_m( parent_mode.m )
    return self
