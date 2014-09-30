#!python_alias
usage="""an executable to grow networks by adding daughters. "b" is for "building" """

import sys, glob, copy
import nmode_utils as nm_u
import nmode_state as nm_s
import prune_network as pn
import numpy as np
import mode_selection as ms
import ggg

import ConfigParser

from optparse import OptionParser

####################################################################################################
#
#
# helper functions
#
#
####################################################################################################
def new_parent_selection(config, system, q=False, verbose=False):
  """
  takes the config file and selects new parents based on the criteria therein
  returns [(O,mode), ...]
  """
#  selection = dict( config.items("parents") )
  selection = dict( [(key.lower(),value) for key, value in config.items("parents")] )
  network = system.network
  freqs = system.compute_3mode_freqs()

  if verbose: print "selection : ", selection

  ### modeNos
  if selection.has_key("modeNos".lower()):
    if verbose: print "selecting parents based on modeNo"
    modeNos = [int(l) for l in selection["modeNos".lower()].strip().split()]
  else:
    modeNos = range(len(network))

  ### genNos
  if selection.has_key("genNos".lower()):
    _modeNos = []
    if verbose: print "selecting parents based on genNo"
    gens,_ = network.gens()
    for genNo in [int(l) for l in selection['genNos'.lower()].split()]:
      for modeNo in gens[genNo]:
        if modeNo in modeNos:
          _modeNos.append( modeNo )
    modeNos = _modeNos

  ### l
  has_min = selection.has_key("min_l".lower())
  has_max = selection.has_key("max_l".lower())
  if has_min or has_max:
    if verbose: print "selecting parents based on l"
    _modeNos = []
    for modeNo in modeNos:
      l = network.modes[modeNo].l
      if (not has_min) or (has_min and (l >= int(selection["min_l".lower()]))):
        if (not has_max) or (has_max and (l <= int(selection["max_l".lower()]))):
          _modeNos.append( modeNo )
    modeNos = _modeNos

  ### frac_Oorb
  has_min = selection.has_key("min_frac_Oorb".lower())
  has_max = selection.has_key("max_frac_Oorb".lower())
  if has_min or has_max:
    if verbose: print "selecting parents based on frac_Oorb"
    _modeNos = []
    for modeNo in modeNos:
      absO_Oorb = abs(freqs[modeNo]/system.Oorb)
      if (not has_min) or (has_min and (absO_Oorb >= float(selection["min_frac_Oorb".lower()]))):
        if (not has_max) or (has_max and (absO_Oorb <= float(selection["max_frac_Oorb".lower()]))):
          _modeNos.append( modeNo )
    modeNos = _modeNos

  ### Amean
  if selection.has_key("min_Amean".lower()):
    if verbose: print "selecting parents based on min_Amean"
    keep, _ = pn.downselect_Amean(q, float(selection["min_Amean".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
  if selection.has_key("max_Amean"):
    if verbose: print "selecting parents based on max_Amean"
    _, keep = pn.downselect_Amean(q, float(selection["max_Amean".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
    
  ### Alast
  if selection.has_key("min_Alast".lower()):
    if verbose: print "selecting parents based on min_Alast"
    keep, _ = pn.downselect_Alast(q, float(selection["min_Alast".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
  if selection.has_key("max_Alast"):
    if verbose: print "selecting parents based on max_Alast"
    _, keep = pn.downselect_Alast(q, float(selection["max_Alast".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]

  ### Afit
  if selection.has_key("min_Afit".lower()):
    if verbose: print "selecting parents based on min_Afit"
    keep, _ = pn.downselect_Afit(q, float(selection["min_Afit".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
  if selection.has_key("max_Afit"):
    if verbose: print "selecting parents based on max_Afit"
    _, keep = pn.downselect_Afit(q, float(selection["max_Afit".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]

  ### num_k
  if selection.has_key("min_num_k".lower()):
    if verbose: print "selecting parents based on min_num_k"
    keep, _ = pn.downselect_num_k(q, float(selection["min_num_k".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
  if selection.has_key("max_num_k"):
    if verbose: print "selecting parents based on max_num_k"
    _, keep = pn.downselect_num_k(q, float(selection["max_num_k".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]

  ### min_heuristic
  if selection.has_key("min_heuristic".lower()):
    if verbose: print "selecting parents based on min_heuristic"
    keep, _ = pn.downselect_min_heuristic(q, float(selection["min_heuristic".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
  if selection.has_key("max_heuristic"):
    if verbose: print "selecting parents based on max_heuristic"
    _, keep = pn.downselect_min_heuristic(q, float(selection["max_heuristic".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]

  ### min_Ethr
  if selection.has_key("min_Ethr".lower()):
    if verbose: print "selecting parents based on min_Ethr"
    keep, _ = pn.downselect_min_Ethr(q, float(selection["min_Ethr".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
  if selection.has_key("max_Ethr"):
    if verbose: print "selecting parents based on max_Ethr"
    _, keep = pn.downselect_min_Ethr(q, float(selection["max_Ethr".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]

  ### collE
  if selection.has_key("min_collE".lower()):
    if verbose: print "selecting parents based on min_collE"
    keep, _ = pn.downselect_collE(q, float(selection["min_collE".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]
  if selection.has_key("max_collE"):
    if verbose: print "selecting parents based on max_collE"
    _, keep = pn.downselect_collE(q, float(selection["max_collE".lower()]), network, mode_nums=modeNos)
    modeNos = [network.modeNoD[nlms] for nlms in keep]

  return [(freqs[modeNo], network.modes[modeNo]) for modeNo in sorted(modeNos)]

####################################################################################################
#
#
#                        Parse input options
#
#
####################################################################################################

parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-F", "--outfilename", default=False, type="string", help="an integration output file used to compute Eo for collective daughter modes")

parser.add_option("-l", "--logfilename", default=False, type="string", help="the log file containing a network you wish to manipulate.")
parser.add_option("-L", "--new-logfilename", default=False, type="string", help="the log file into which the new network will be written. This should be a base on which the number of daughter modes will be added.")

parser.add_option("-c", "--config", default=False, type="string", help="the config file for adding new modes")

### network parameters
#parser.add_option("", "--alpha", default=4e-3, type="float", help="alpha for gmode generation")
#parser.add_option("", "--c", default=2e-11, type="float", help="c for gmode generation")
#parser.add_option("", "--k-hat", default=5e4, type="float", help="k_hat for 3mode coupling coefficients")

### add to an existing network based on integration output
#parser.add_option("", "--freqfilename", default=False, type="string", help="the file containing freq-domain data")
#parser.add_option("", "--min-f_forb", default=False, type="float", help="minimum frequency read from freqfilename")
#parser.add_option("", "--max-f_forb", default=False, type="float", help="maximum frequency read form freqfilename")
#parser.add_option("", "--freqfile-downsample", default=False, type="float", help="only read in one out of this many entries from freqfilename")

#parser.add_option("", "--daughter-min-frac-Oorb", default=-2.0, type="float", help="the minimum fraction of the orbital frequency above which we define \"new parents\"")
#parser.add_option("", "--daughter-max-frac-Oorb", default= 2.0, type="float", help="the maximum fraction fo the orbital frequency below which we define \"new parents\"")

#parser.add_option("", "--daughter-selection", default=False, type="string", help="the daughter selection method used: [Ethr, heuristic, min_Ethr, min_heuristic, collective]")
#parser.add_option("", "--num-pairs", default="1", type="string", help="the number of couplings to add to a network. For multiple values, supply a space-delimited string")
#parser.add_option("", "--daughter-min-l", default=1, type="int", help="the minimum angular momentum quantum number for daughter modes.")
#parser.add_option("", "--daughter-max-l", default=4, type="int", help="the maximum angular momentum quantum number for daughter modes.")
#parser.add_option("", "--daughter-min-frac-w", default=0.0, type="float", help="the minimum fraction of the parent's oscillation frequency above which we look for daughter modes.")
#parser.add_option("", "--daughter-max-frac-w", default=1.0, type="float", help="the maximum fraction of the parent's oscillation frequency below which we look for daughter modes.")
parser.add_option("", "--Emax", default=np.infty, type="float", help="the maximum Ethr that will be recorded in any new coupling list ASCii files")
parser.add_option("", "--maxp", default=1, type="int", help="the maximum number of coupling list generation jobs that will be submitted in parallel")

#parser.add_option("", "--to-unique-couplings", default=False, action="store_true", help="makes coupling lists explicitly check and remove redundant couplings")

parser.add_option("", "--max-num-pairs", default=50000, type="int", help="the maximum number of couplings that will be written into new ctg files. DEFAULT=50000")
parser.add_option("", "--catalog-dir", default="./", type="string", help="path to the directory in which catalog files (mode lists) are stored")

#parser.add_option("", "--parent-forcing", default=False, action="store_true", help="if True, we add forcing terms to the 'new parents'")
#parser.add_option("", "--daughter-forcing", default=False, action="store_true", help="if True, we add forcing terms to the 'new daughters'")
#parser.add_option("", "--intercouple", default=False, action="store_true", help="if True, we compute and include all possible couplings in this network.")

opts, args = parser.parse_args()

if not opts.config:
  opts.config = raw_input("config = ")

if not opts.new_logfilename:
  opts.new_logfilename = raw_input("new logfilename = ")

if not opts.logfilename:
  opts.logfilename = raw_input("old logfilename = ")

### read in config
if opts.verbose: print "reading config from %s" % opts.config
config = ConfigParser.ConfigParser()
config.read(opts.config)

### requre outfilename with certain parent selection criteria
selection = dict( config.items("parents") )
require_outfile = selection.has_key("min_Amean") or selection.has_key("max_Amean") or selection.has_key("min_Alast") or selection.has_key("max_Alast") or selection.has_key("min_Afit") or selection.has_key("max_Afit")
if require_outfile:
  if not opts.outfilename:
    opts.outfilename = raw_input("outfilename = ")

### get the number of daughters requested into a convenient form
num_pairs = sorted([int(l) for l in config.get("children","num_pairs").strip().split()])

if "collective" in config.get("children","selection"):
  if config.get("children","selection") == "collective":
    if not opts.outfilename:
      opts.outfilename = raw_input("outfilename = ")
    min_num_pairs = config.getint("children","min_num_pairs")

  if len(num_pairs) > 1:
    raise ValueError, "collective instabilities currently only support a single num-pair argument."
  num_pairs = num_pairs[:1]

####################################################################################################
#
#
#                              grow the network based on integration data
#                                        (daughter selection)
#
####################################################################################################

if opts.verbose: print "reading in network parameters from %s" % opts.logfilename
system = nm_u.load_log(opts.logfilename)
Mprim = system.Mprim
Mcomp = system.Mcomp
Rprim = system.Rprim
Porb = system.Porb
Oorb = 2*np.pi/Porb
eccentricity = system.eccentricity

wo = ms.compute_wo(Mprim, Rprim)

### general parameters
alpha = config.getfloat("general", "alpha")
c = config.getfloat("general", "c")
k_hat = config.getfloat("general", "k_hat")

if opts.outfilename and require_outfile:
  if opts.verbose: print "reading integration data from ", opts.outfilename
  t_P, q, N_m = nm_u.load_out_f(opts.outfilename)
else:
  q = False

### new parent selection
modes = new_parent_selection(config, system, q=q, verbose=opts.verbose)

### make new system objects
new_systems = [copy.deepcopy(system) for npairs in num_pairs]

### new child selection
children = dict( config.items("children") )
daughter_selection = children["selection"]

### default values
if children.has_key("min_l"):
  daughter_min_l = int(children["min_l"])
else:
  daughter_min_l = 1
  if opts.verbose: print "using DEFAULT daughter_min_l=%d"%daughter_min_l

if children.has_key("max_l"):
  daughter_max_l = int(children["max_l"])
else:
  daughter_max_l = 100
  if opts.verbose: print "using DEFAULT daughter_max_l=%d"%daughter_max_l

if children.has_key("min_frac_w"):
  daughter_min_frac_w = float(children["min_frac_w"])
else:
  daughter_min_frac_w = 0.1
  if opts.verbose: print "using DEFAULT daughter_min_frac_w=%.3f"%daughter_min_frac_w

if children.has_key("max_frac_w"):
  daughter_max_frac_w = float(children["max_frac_w"])
else:
  daughter_max_frac_w = 0.9
  if opts.verbose: print "using DEFAULT daughter_max_frac_w=%.3f"%daughter_max_frac_w

parent_forcing = config.getboolean("inclusion","parent_forcing")
daughter_forcing = config.getboolean("inclusion","daughter_forcing")
intercouple = config.getboolean("inclusion","intercouple")
to_unique_couplings = config.getboolean("inclusion","to_unique_couplings")

### compute the threshold energies only once.
if daughter_selection == "collective":
  if opts.verbose: print "reading in integration data from %s and computing energies" % opts.outfilename
  E = [np.mean(a)**2 for a in nm_s.compute_A(nm_u.load_out(opts.outfilename)[1], Eo=1.0)] # growth rate depends on amplitude, so that's what we should average (not energy)
  
new_systems = [copy.deepcopy(system) for npairs in num_pairs]

# iterate over new parents
n_modes = len(modes)
for n_m, (O, mode) in enumerate(modes):
  n, l, m, w, y = mode.get_nlmwy()
  if opts.verbose: print "%d / %d :finding daughter modes\n\tfor mode: %s\n\toscillating at %fe-5 rad/%s\n\tusing %s" % (n_m+1, n_modes, mode.to_str_nlmwy(), O*1e5, nm_u.units['time'], daughter_selection)

  absO = abs(O)
  signO = O/absO
  min_w = daughter_min_frac_w*absO
  max_w = daughter_max_frac_w*absO

  if daughter_selection == "collective":  ### collective modes
    modeNo = system.network.modeNoD[mode.get_nlms()]
    ### identify mode number in current network
    my_triples = ggg.multiple_collective_instabilities(mode, O, E[modeNo], maxp=opts.maxp, Nmin=min_num_pairs, Nmax=num_pairs[0], alpha=alpha, c=c, wo=wo, k_hat=k_hat, verbose=opts.verbose, min_l=daughter_min_l, max_l=daughter_max_l, min_absw=min_w, max_absw=max_w, catalogdir=opts.catalog_dir)
    if opts.verbose: 
      print "found %d triples" % (len(my_triples))
    ### add triples to the network
    new_systems[0].network.add_couplings(my_triples, opts.verbose)

  elif daughter_selection == "fast_collective": ### collective sets a la Weinberg 2012
    my_triples = ggg.fast_collective_instability(mode, O, min_l=daughter_min_l, max_l=daughter_max_l, alpha=alpha, c=c, wo=wo, k_hat=k_hat, Nmax=num_pairs[0])
    if opts.verbose: print "found %d triples" % len(my_triples)
    new_systems[0].network.add_couplings(my_triples, opts.verbose)

  else: # 3mode selection criteria
    # parameters are set, we now call ggg.compute_pairs()
    useful_filenames = ggg.compute_pairs(daughter_selection, mode, O, opts.catalog_dir, min_l=daughter_min_l, max_l=daughter_max_l, min_w=min_w, max_w=max_w, alpha=alpha, c=c, wo=wo, k_hat=k_hat, Emax=opts.Emax, maxp=opts.maxp, verbose=opts.verbose, max_num_pairs=opts.max_num_pairs)

    # read in data from the useful_filenames
    if opts.verbose: 
      print "reading in %d couplings from:" % max(num_pairs)
      for filename in useful_filenames:
        print "\t%s" % filename
  
    if "min_" == daughter_selection[0:4]:
      my_coupling_list = ggg.ggg_minima_coupling_list(alpha, c, wo, k_hat, parent_mode=mode).load_unsorted_mode_lists(daughter_selection, useful_filenames, min_n=False, max_n=False, min_l=daughter_min_l, max_l=daughter_max_l, min_w=min_w, max_w=max_w)
      if to_unique_couplings:
        if opts.verbose: print "checking for unique couplings"
        my_coupling_list = my_coupling_list.to_unique_couplings()

      for ind, npairs in enumerate(num_pairs):
        if opts.verbose: print "grabbing %d triples" % npairs
        my_triples = my_coupling_list.to_triples(npairs, min_n=False, max_n=False, min_l=daughter_min_l, max_l=daughter_max_l, min_w=min_w, max_w=max_w, parent_forcing=parent_forcing, daughter_forcing=daughter_forcing, Mprim=Mprim, Mcomp=Mcomp, Porb=Porb, eccentricity=eccentricity)
        ### add triples to the network
        new_systems[ind].network.add_couplings(my_triples, opts.verbose)

    else:
      my_coupling_list = ggg.ggg_coupling_list(alpha, c, wo, k_hat, parent_mode=mode).load_unsorted_mode_lists(daughter_selection, useful_filenames, num_pairs=max(num_pairs), min_n=False, max_n=False, min_l=daughter_min_l, max_l=daughter_max_l, min_w=min_w, max_w=max_w)
      if to_unique_couplings:
        if opts.verbose: print "checking for unique couplings"
        my_coupling_list = my_coupling_list.to_unique_couplings()

      for ind, npairs in enumerate(num_pairs):
        if opts.verbose: print "grabbing %d triples" % npairs
        my_triples = my_coupling_list.to_triples(npairs, parent_forcing=parent_forcing, daughter_forcing=daughter_forcing, Mprim=Mprim, Mcomp=Mcomp, Porb=Porb, eccentricity=eccentricity)
        ### add triples to the network
        new_systems[ind].network.add_couplings(my_triples, opts.verbose)

# write new network to disk
for ind, npairs in enumerate(num_pairs):
  new_logfilename = opts.new_logfilename+"_%d.log" % npairs

  if intercouple:
    if opts.verbose: print "computing all possible couplings for %s" % new_logfilename
    new_systems[ind].network = ggg.intercouple_network(new_systems[ind].network, k_hat=k_hat, verbose=opts.verbose)

  if opts.verbose: print "writing newtork to %s" % new_logfilename
  nm_u.write_log(new_logfilename, new_systems[ind])


