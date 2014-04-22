#!python_alias
usage="""an executable to grow networks by adding daughters. "b" is for "building" """

import sys, glob, copy
import nmode_utils as nm_u
import nmode_state as nm_s
import prune_network as pn
import numpy as np
import mode_selection as ms
import ggg

from optparse import *

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

### network parameters
parser.add_option("", "--alpha", default=4e-3, type="float", help="alpha for gmode generation")
parser.add_option("", "--c", default=2e-11, type="float", help="c for gmode generation")
parser.add_option("", "--k-hat", default=5e4, type="float", help="k_hat for 3mode coupling coefficients")

### add to an existing network based on integration output
#parser.add_option("", "--freqfilename", default=False, type="string", help="the file containing freq-domain data")
#parser.add_option("", "--min-f_forb", default=False, type="float", help="minimum frequency read from freqfilename")
#parser.add_option("", "--max-f_forb", default=False, type="float", help="maximum frequency read form freqfilename")
#parser.add_option("", "--freqfile-downsample", default=False, type="float", help="only read in one out of this many entries from freqfilename")

parser.add_option("", "--daughter-min-frac-Oorb", default=-2.0, type="float", help="the minimum fraction of the orbital frequency above which we define \"new parents\"")
parser.add_option("", "--daughter-max-frac-Oorb", default= 2.0, type="float", help="the maximum fraction fo the orbital frequency below which we define \"new parents\"")

parser.add_option("", "--daughter-selection", default=False, type="string", help="the daughter selection method used: [Ethr, heuristic, min_Ethr, min_heuristic, collective]")
parser.add_option("", "--num-pairs", default="1", type="string", help="the number of couplings to add to a network. For multiple values, supply a space-delimited string")
parser.add_option("", "--daughter-min-l", default=1, type="int", help="the minimum angular momentum quantum number for daughter modes.")
parser.add_option("", "--daughter-max-l", default=4, type="int", help="the maximum angular momentum quantum number for daughter modes.")
parser.add_option("", "--daughter-min-frac-w", default=0.0, type="float", help="the minimum fraction of the parent's oscillation frequency above which we look for daughter modes.")
parser.add_option("", "--daughter-max-frac-w", default=1.0, type="float", help="the maximum fraction of the parent's oscillation frequency below which we look for daughter modes.")
parser.add_option("", "--Emax", default=np.infty, type="float", help="the maximum Ethr that will be recorded in any new coupling list ASCii files")
parser.add_option("", "--maxp", default=1, type="int", help="the maximum number of coupling list generation jobs that will be submitted in parallel")

parser.add_option("", "--max-num-pairs", default=50000, type="int", help="the maximum number of couplings that will be written into new ctg files. DEFAULT=50000")
parser.add_option("", "--catalog-dir", default="./", type="string", help="path to the directory in which catalog files (mode lists) are stored")

parser.add_option("", "--parent-forcing", default=False, action="store_true", help="if True, we add forcing terms to the 'new parents'")
parser.add_option("", "--daughter-forcing", default=False, action="store_true", help="if True, we add forcing terms to the 'new daughters'")

parser.add_option("", "--intercouple", default=False, action="store_true", help="if True, we compute and include all possible couplings in this network.")

opts, args = parser.parse_args()

if not opts.new_logfilename:
  opts.new_logfilename = raw_input("new logfilename = ")

### daughter selection
if not opts.logfilename:
  opts.logfilename = raw_input("old logfilename = ")
#if not opts.freqfilename:
#  opts.freqfilename = raw_input("freq-domain filename = ")
if not opts.daughter_selection:
  opts.daughter_selection = raw_input("daughter_selection = ")

num_pairs = sorted([int(l) for l in opts.num_pairs.split()])

if opts.daughter_selection == "collective":
  if not opts.outfilename:
    opts.outfilename = raw_input("outfilename = ")
  if len(num_pairs) > 1:
    raise ValueError, "collective instabilities currently only support a single num-pair arguement."
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

#if opts.verbose: print "reading frequency data from %s" % opts.freqfilename
#freq, fq, N_m = nm_u.load_out(opts.freqfilename, tmin=opts.min_f_forb, tmax=opts.max_f_forb, downsample=opts.freqfile_downsample)

if opts.verbose: print "identifying modes within bandwidth"
min_w = opts.daughter_min_frac_Oorb*Oorb
max_w = opts.daughter_max_frac_Oorb*Oorb

modes = pn.within_bandwidth_analytic(min_w, max_w, system)
#modes = pn.within_bandwidth_maxPSD(freq, fq, min_w, max_w, system.network, Oorb=Oorb) # frequencies in "freq" are stored as f/forb
#modes = pn.within_bandwidth_lorentzian(freq, fq, min_w, max_w, system.network, Oorb=Oorb)

if opts.daughter_selection == "collective":
  if opts.verbose: print "reading in integration data from %s and computing energies" % opts.outfilename
#  E = [np.max(a)**2 for a in nm_s.compute_A(nm_u.load_out(opts.outfilename)[1], Eo=1.0)]
#  E = [np.min(a)**2 for a in nm_s.compute_A(nm_u.load_out(opts.outfilename)[1], Eo=1.0)]
  E = [np.mean(a)**2 for a in nm_s.compute_A(nm_u.load_out(opts.outfilename)[1], Eo=1.0)] # growth rate depends on amplitude, so that's what we should average (not energy)
  
new_systems = [copy.deepcopy(system) for npairs in num_pairs]

# iterate over new parents
n_modes = len(modes)
for n_m, (O, mode) in enumerate(modes):
  n, l, m, w, y = mode.get_nlmwy()
  if opts.verbose: print "%d / %d :finding daughter modes\n\tfor mode: %s\n\toscillating at %fe-5 rad/%s\n\tusing %s" % (n_m+1, n_modes, mode.to_str_nlmwy(), O*1e5, nm_u.units['time'], opts.daughter_selection)

  absO = abs(O)
  signO = O/absO
  min_w = opts.daughter_min_frac_w*absO
  max_w = opts.daughter_max_frac_w*absO

  if opts.daughter_selection == "collective":
    modeNo = system.network.modeNoD[mode.get_nlms()]
    ### identify mode number in current network
    my_triples = ggg.multiple_collective_instabilities(mode, O, E[modeNo], maxp=opts.maxp, Nmin=0, Nmax=num_pairs[0], alpha=opts.alpha, c=opts.c, wo=wo, k_hat=opts.k_hat, verbose=opts.verbose, min_l=opts.daughter_min_l, max_l=opts.daughter_max_l, min_absw=min_w, max_absw=max_w)
    if opts.verbose: 
      print "found %d triples" % (len(my_triples))
    ### add triples to the network
    new_systems[0].network.add_couplings(my_triples, opts.verbose)

  else: # 3mode selection criteria
    # parameters are set, we now call ggg.compute_pairs()
    useful_filenames = ggg.compute_pairs(opts.daughter_selection, mode, O, opts.catalog_dir, min_l=opts.daughter_min_l, max_l=opts.daughter_max_l, min_w=min_w, max_w=max_w, alpha=opts.alpha, c=opts.c, wo=wo, k_hat=opts.k_hat, Emax=opts.Emax, maxp=opts.maxp, verbose=opts.verbose, max_num_pairs=opts.max_num_pairs)

    # read in data from the useful_filenames
    if opts.verbose: 
      print "reading in %d couplings from:" % max(num_pairs)
      for filename in useful_filenames:
        print "\t%s" % filename
  
    if "min_" == opts.daughter_selection[0:4]:
      my_coupling_list = ggg.ggg_minima_coupling_list(opts.alpha, opts.c, wo, opts.k_hat, parent_mode=mode).load_unsorted_mode_lists(opts.daughter_selection, useful_filenames, min_n=False, max_n=False, min_l=opts.daughter_min_l, max_l=opts.daughter_max_l, min_w=min_w, max_w=max_w).to_unique_couplings()

      for ind, npairs in enumerate(num_pairs):
        if opts.verbose: print "grabbing %d triples" % npairs
        my_triples = my_coupling_list.to_triples(npairs, min_n=False, max_n=False, min_l=opts.daughter_min_l, max_l=opts.daughter_max_l, min_w=min_w, max_w=max_w, parent_forcing=opts.parent_forcing, daughter_forcing=opts.daughter_forcing, Mprim=Mprim, Mcomp=Mcomp, Porb=Porb, eccentricity=eccentricity)
        ### add triples to the network
        new_systems[ind].network.add_couplings(my_triples, opts.verbose)


    else:
      my_coupling_list = ggg.ggg_coupling_list(opts.alpha, opts.c, wo, opts.k_hat, parent_mode=mode).load_unsorted_mode_lists(opts.daughter_selection, useful_filenames, num_pairs=max(num_pairs), min_n=False, max_n=False, min_l=opts.daughter_min_l, max_l=opts.daughter_max_l, min_w=min_w, max_w=max_w).to_unique_couplings()

      for ind, npairs in enumerate(num_pairs):
        if opts.verbose: print "grabbing %d triples" % npairs
        my_triples = my_coupling_list.to_triples(npairs, parent_forcing=opts.parent_forcing, daughter_forcing=opts.daughter_forcing, Mprim=Mprim, Mcomp=Mcomp, Porb=Porb, eccentricity=eccentricity)
        ### add triples to the network
        new_systems[ind].network.add_couplings(my_triples, opts.verbose)

# write new network to disk
for ind, npairs in enumerate(num_pairs):
  new_logfilename = opts.new_logfilename+"_%d.log" % npairs

  if opts.intercouple:
    if opts.verbose: print "computing all possible couplings for %s" % new_logfilename
    new_systems[ind].network = ggg.intercouple_network(new_systems[ind].network, k_hat=opts.k_hat, verbose=opts.verbose)

  if opts.verbose: print "writing newtork to %s" % new_logfilename
  nm_u.write_log(new_logfilename, new_systems[ind])


