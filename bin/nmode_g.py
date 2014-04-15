#!/usr/local/bin/python

usage="""an executable to generate new networks (parent modes). "g" is for "generate" """

import sys, glob
import nmode_utils as nm_u
import prune_network as pn
import numpy as np
import mode_selection as ms

import networks, ggg

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

parser.add_option("-L", "--new-logfilename", default=False, type="string", help="the log file into which the new network will be written.")

### network parameters
parser.add_option("", "--Mprim", default=1.0, type="float", help="Mprim in units of Msun")
parser.add_option("", "--Mcomp", default=1.0, type="float", help="Mcomp in units of Mjup")
parser.add_option("", "--Rprim", default=1.0, type="float", help="Rprim in units of Rsun")
parser.add_option("", "--Porb", default=3*86400, type="float", help="the orbital period in seconds")
parser.add_option("", "--eccentricity", default=0.0, type="float", help="the orbital eccentricity")

parser.add_option("", "--alpha", default=4e-3, type="float", help="alpha for gmode generation")
parser.add_option("", "--c", default=2e-11, type="float", help="c for gmode generation")
parser.add_option("", "--k-hat", default=5e4, type="float", help="k_hat for 3mode coupling coefficients")

### generate a new network (find parents)
parser.add_option("", "--parent-selection", default=False, type="string", help="the algorithm used to identify parent modes: Elin or detuning")
parser.add_option("", "--num-parents", default=2, type="int", help="the minimum number of parents to be added to the network")
parser.add_option("", "--parent-lm", default=False, type="string", help="the values of l and m allowed for parent modes. This should have the format: \"l1:m11,m12 l2:m21,m22\"")
parser.add_option("", "--parent-min-frac-Oorb", default=0.0, type="float", help="minimum fraction of orbital frequency for parent selection")
parser.add_option("", "--parent-max-frac-Oorb", default=4.0, type="float", help="maximum fraction of orbital frequency for parent selection")

opts, args = parser.parse_args()

if not opts.new_logfilename:
  opts.new_logfilename = raw_input("new logfilename = ")

### ensure that all required fields are present for parent selection
if not opts.parent_lm:
  opts.parent_lm = raw_input("parent lm = ")
# set up bounds for parent selection
bounds = {}
for item in opts.parent_lm.split():
  _l, _ms = item.split(":")
  bounds[int(_l)] = [ int(l) for l in _ms.split(",") ]

### check eccentricity value
if (opts.eccentricity < 0) or (opts.eccentricity > 1):
  sys.exit("bad value for eccentricity. Please supply a value between 0 and 1")

##### WE CURRENTLY ONLY SUPPORT eccentricity=0 SYSTEMS ##########
if opts.eccentricity != 0:
  sys.exit("WARNING!: we currently only support eccentricity == 0 sytems")

####################################################################################################
#
#
#                                      generate a new network
#                                        (parent selection)
#
####################################################################################################
if opts.verbose: print "generating a new network (parent selection)"

wo = ms.compute_wo(opts.Mprim, opts.Rprim)

Oorb = 2*np.pi/opts.Porb
system = networks.system(opts.Mprim, opts.Mcomp, opts.Rprim, opts.Porb, opts.eccentricity)
network = system.network

min_w = opts.parent_min_frac_Oorb*Oorb
max_w = opts.parent_max_frac_Oorb*Oorb

if opts.verbose: print "computing parent modes using %s" % opts.parent_selection
if opts.parent_selection == "Elin":
  parents = ggg.compute_parents_Elin(Oorb, bounds, N=opts.num_parents, min_w=min_w, max_w=max_w, alpha=opts.alpha, c=opts.c, wo=wo, Mprim=opts.Mprim, Mcomp=opts.Mcomp, Porb=opts.Porb, eccentricity=opts.eccentricity)
elif opts.parent_selection == "detuning":
  parents = ggg.compute_parents_detuning(Oorb, bounds, N=opts.num_parents, min_w=min_w, max_w=max_w, alpha=opts.alpha, c=opts.c, wo=wo, forcing=True, Mprim=opts.Mprim, Mcomp=opts.Mcomp, Porb=opts.Porb, eccentricity=opts.eccentricity)
else:
  sys.exit("unknown parent selection algorithm: %s\nplease supply a parent selection algorithm from the following:\n\tEthr\n\tdetuning" % (opts.parent_selection) )

if opts.verbose: print "attempting to add %d parents to the network" % len(parents)
network = network.add_modes(parents, verbose=opts.verbose)

if opts.verbose: print "writing network to %s" % opts.new_logfilename
nm_u.write_log(opts.new_logfilename, system)
  

