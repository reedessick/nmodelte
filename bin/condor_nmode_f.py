#!python_alias
usage="""a wrapper for nmode_f for condor submission. This will help us spare the nfs from constant output. "f" is for "fast integration" """

import subprocess

from optparse import *

####################################################################################################
#
#
#                        Parse input options
#
#
####################################################################################################

parser=OptionParser(usage=usage)

parser.add_option("", "--tmpfilename", default=False, type="string", help="the temporary file into which we write to spare the nfs from our persistent output")

parser.add_option("", "--outfilename", default=False, type="string", help="file into which integration data is written. If not supplied, data is printed to the terminal.")

### the rest of these are simply passed to nmode_f.py

# options about integration
parser.add_option("-P", "--num-periods", dest="N_P", default=False, type="float", help="the number of periods for simulation")
parser.add_option("-f", "--function", dest="function", default="", type="string", help="the name of the derivative function used in integration. Names come from network_flow.py")
parser.add_option("-n", "--num-proc", dest="num_proc", default=1, type="int", help="the number of processors available for parallel computation.")
parser.add_option("", "--steps-per-P", default=1, type="float", help="report the system's state this many times per period.")

parser.add_option("", "--rel-err", default=1e-6, type="float", help="relative error for integrator")
parser.add_option("", "--abs-err", default=0.0, type="float", help="absolute error for integrator")

parser.add_option("", "--equilib-IC", default=False, type="string", help="either \"lin\" for \"3mode\" and will compute the associated equilibrium amplitudes for the system and use them as initial conditions. WILL BE IGNORED IF --init-conds IS SUPPLIED.")
parser.add_option("", "--init-conds", default=False, type="string", help="set initial conditions for the integration. Supply a string separated by spaces. WILL BE IGNORED IF --onward IS SUPPLIED.")
parser.add_option("", "--init-time", default="none", type="string", help="set initial time for the integration in units of orbital periods.")
parser.add_option("", "--default-value", default=1e-12, type="float", help="the initial value assigned to real&imag parts of a mode when no other information is supplied.")

parser.add_option("", "--onward", default=False, type="string", help="continue an integration where it left off. supply the filename of the output file and use -L for the log-file.")
parser.add_option("", "--onward-logfilename", default=False, type="string", help="use the last line of the file specified by --onward and this log file to map old network's final conditions into the new network's initial conditions. If there are modes in the new network that were not present before, they start from a default value.")

parser.add_option("-l", "--logfilename", dest="logfilename", default=False, type="string", help="the system parameters will be taken from this log file.")

parser.add_option("-t", "--time-integration", default=False, action="store_true", help="time the duration of the integration loop and report the result. This does not include overhead like setting up initial conditions or integration objects")

parser.add_option("", "--flush", default=False, action="store_true", help="flush outputfile after each print statment")

opts, args = parser.parse_args()

if not opts.outfilename:
	opts.outfilename = raw_input("outfilenmae = ")

if not opts.tmpfilename:
	opts.tmpfilename = opts.outfilename

####################################################################################################
### build command for nmode_f.py

cmd = "nmode_f.py --outfilename %s -l %s -f %s --num-proc %d -P %f --steps-per-P %f --rel-err %s --abs-err %s "%(opts.tmpfilename, opts.logfilename, opts.function, opts.num_proc, opts.num_periods, opts.steps_per_P, str(opts.rel_err), str(opts.abs_err) )

cmd = cmd.split()

if opts.equilib_IC:
	cmd += ["--equilib-IC", opts.equilib_IC]
if opts.init_conds:
	cmd += ["--init-conds", opts.init_conds]
if opts.default_value:
	cmd += ["--default-value", str(opts.default_value)]
if opts.onward:
	cmd += ["--onward", opts.onward]
if opts.onward_logfilename:
	cmd += ["--onward-logfilename", opts.onward]
if opts.time-integration:
	cmd += ["--time-integration"]
if opts.flush:
	cmd += ["--flush"]

### launch command through subprocess
p = subprocess.Popen(cmd) ### instantiate process
exitcode = p.wait() ### wait for procesis to finish
if exitcode:
	import sys
	sys.exit(exitcode)

### move tmpfilename to outfilename
if opts.tmpfilename != opts.outfilename:
	if opts.onward and (opts.outfilename==opts.onward): 
 		### cat opts.tmpfilename >> opts.outfilename
		outfile_obj = open(opts.outfilename, "a") ### we want to append
		subprocess.call(["cat", opts.tmpfilename], stdout=outfile_obj)
		outfile_obj.close()
		### rm opts.tmpfilename
		subprocess.call(["rm", opts.tmpfilename])
	else:
		### mv opts.tmpfilename opts.outfilename
		subprocess.call(["mv", opts.tmpfilename, opts.outfilename])

    



