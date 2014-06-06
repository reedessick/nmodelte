#!python_alias
usage="""an executable to remove modes from a network based on integration data. "r" is for "remove" """

import sys, glob
import nmode_utils as nm_u
import prune_network as pn
import numpy as np
import mode_selection as ms

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

parser.add_option("-l", "--logfilename", default=False, type="string", help="the log file containing a network you wish to manipulate.")
parser.add_option("-L", "--new-logfilename", default=False, type="string", help="the log file into which the new network will be written.")

### remove modes from the network based on amplitude cuts, etc.
parser.add_option("", "--filename", default=False, type="string", help="the file containing time-domain data")
parser.add_option("", "--min-t_Porb", default=False, type="float", help="minimum t/Porb read from filename")
parser.add_option("", "--max-t_Porb", default=False, type="float", help="maximum t/Porb read from filename")
parser.add_option("", "--file-downsample", default=False, type="float", help="only read in one out of this many entries from filename")

parser.add_option("", "--downselection-methods", default=False, type="string", help="downselection methods used: Alast or Amean or Afit")
parser.add_option("", "--Athrs", default=False, type="string", help="amplitude thresholds for downselection methods. These must be supplied in the same order as downselection-methods")

opts, args = parser.parse_args()

if not opts.new_logfilename:
  opts.new_logfilename = raw_input("new logfilename = ")


### downselection
downselection_methods = opts.downselection_methods.split()
if not opts.logfilename:
  opts.logfilename = raw_input("old logfilename = ")

if not opts.filename:
  opts.filename = raw_input("time-domain filename = ")

Athrs = [float(l) for l in opts.Athrs.split()]
if len(Athrs) != len(downselection_methods):
  print "len(Athrs) conflicts with len(downselection_methods). Please supply one Athr for each of the following:"
  Athrs = []
  for method in downselection_methods:
    Athrs.append( float( raw_input("Athr for %s = " % method) ) )

####################################################################################################
#
#
#                             network downselection based on integration data
#
#
####################################################################################################

if opts.verbose: print "downselecting a subset of modes based on integration data"

if opts.verbose: print "reading in network parameters from %s" % opts.logfilename
system = nm_u.load_log(opts.logfilename)
network = system.network
Oorb = system.Oorb
Porb = 2*np.pi/Oorb

if opts.verbose: print "reading time-domain data from %s" % opts.filename
t_P, q, N_m = nm_u.load_out(opts.filename, tmin=opts.min_t_Porb, tmax=opts.max_t_Porb, downsample=opts.file_downsample)

remove = []
for Athr, method in zip(Athrs, downselection_methods):
  if opts.verbose:  _athr, _athrexp = nm_u.float_to_scientific(Athr)
  if method == "Alast":
    if opts.verbose: print "downselecting modes using %s and Athr=%fe%d" % (method, _athr, _athrexp)
    remove += pn.downselect_Alast(q, Athr, network)[1]
  elif method == "Amean":
    if opts.verbose: print "downselecting modes using %s and Athr=%fe%d" % (method, _athr, _athrexp)
    remove += pn.downselect_Amean(q, Athr, network)[1]
  elif method == "Afit":
    if opts.verbose: print "downselecting modes using %s and Athr=%fe%d" % (method, _athr, _athrexp)
    remove += pn.downselect_Afit(q, Athr, network, t_P=t_P)[1]
  else:
    raise ValueError, "unknown downselection method: %s\nplease supply at last one downselection method from the following:\n\tAlst\n\tAmean\n\tAfit"

remove = list(set(remove)) # generate a unique set of modes to remove

if opts.verbose: print "attempting to remove %d/%d modes from the network" % (len(remove), N_m)
network = network.remove_modes(remove, verbose=opts.verbose)

if opts.verbose: print "writing network to %s" % opts.new_logfilename
nm_u.write_log(opts.new_logfilename, system)


