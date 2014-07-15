#!python_alias
usage="""an executable to generate integration data. "f" is for "fast integration" """

import gc

import sys, pygsl
import pygsl.odeiv as odeiv

import network_flow as nf
import nmode_utils as nm_u

import numpy as np

from optparse import *

####################################################################################################
#
#
#                        Parse input options
#
#
####################################################################################################

parser=OptionParser(usage=usage)

parser.add_option("", "--outfilename", default=False, type="string", help="file into which integration data is written. If not supplied, data is printed to the terminal.")

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

opts, args = parser.parse_args()

if opts.outfilename:
  outfile = open(opts.outfilename, "a")
else:
  outfile = sys.stdout

if not opts.logfilename:
  opts.logfilename = raw_input("logfilename = ")

if opts.onward_logfilename and (not opts.onward):
  opts.onward = raw_input("onward filename = ")

if opts.init_time != "none":
  opts.init_time = float(opts.init_time) # stupid fix to get control sequence correct later on

system = nm_u.load_log(opts.logfilename)
Porb = system.Porb
Oorb = 2*np.pi/Porb

N_m = len(system.network) # number of modes
#dimension = np.array(2*(N_m), dtype="i") # dimensionality of the problem
dimension = 2*(N_m) # dimensionality of the problem

if opts.equilib_IC:
  if "dxdt" in opts.function:
    tcurrent = "x"
  elif "dqdt" in opts.function:
    tcurrent = "q"
  else:
    raise ValueError, "could not determine tcurrent from --function=%s"%opts.function

####################################################################################################
#
#
#                          set up initial and termination conditions
#
#
####################################################################################################
### set termination conditions
t1 = opts.N_P*Porb # final time (period times the number of periods)
h = 0.05*Porb # initial step size
sample_step = float(Porb)/opts.steps_per_P # phase between consecutive reports of the intergrator
sample_time = sample_step

### set up initial conditions
if opts.onward:
  f=open(opts.onward, 'r')
  for line in f:
    if line[0] != '#':
      l = line
    else:
      pass
  f.close()
  l = [float(_l) for _l in l.strip().split()]
  t = l[0]*Porb # convert from t/P --> t
  q = l[1:]
  sample_time = t + sample_step

  if opts.onward_logfilename:
    old_system = nm_u.load_log(opts.onward_logfilename)
    map12, map21 = nm_u.modeNo_map([mode.get_nlms() for mode in system.network.modes], [mode.get_nlms() for mode in old_system.network.modes])

    old_q = q[:]
    q = []
    for modeNo in range(N_m):
      old_modeNo = map12[modeNo]
      if old_modeNo == "none":
        q += [ opts.default_value, opts.default_value ]
      else:
        q += [ old_q[2*old_modeNo], old_q[2*old_modeNo+1] ]

  if opts.init_time != "none": # allow the user to restart at t_P=0 if we've changed the network
    t = opts.init_time*Porb
    sample_time = t + sample_step

else:
  if opts.init_time != "none":
    t = opts.init_time*Porb
  else:
    t = 0.0

  if opts.init_conds:
    q = [float(l) for l in opts.init_conds.strip().split()]
    if len(q) == 2:
      print >>outfile, ("# taking %s to be initial conditions for all modes" % str(q))
      q_o = q[:]
      for m in range(1,dimension/2):
        q += q_o
    elif len(q) != dimension:
      raise ValueError, "#IC's supplied do not match number of modes."
  else:
    if opts.equilib_IC: # compute equilibrium state
      if opts.equilib_IC == "lin":
        q = system.compute_lin_eq(t=t, default=opts.default_value, tcurrent="x")
      elif opts.equilib_IC == "3mode":
        q = system.compute_3mode_eq(t=t, default=opts.default_value, tcurrent="x")
      else:
        raise ValueError, "unknown --equilib-IC=%s"%opts.equilib_IC
    else:
      q = np.empty((dimension,)); q[:] = opts.default_value

use_phase = opts.function in ["dxdt_no_NLT_withPHI"]

if use_phase: # phase gets tacked onto the end of the mode amplitude vector
  if t%Porb == 0:
    q = np.concatenate((q, np.array([0])))
  elif t%Porb == 0.5*Porb:
    q = np.concatenate((q, np.array([np.pi])))
  else:
    sys.exit("initial time must be an integer or half-integer multiple of Porb, otherwise we don't know how to initialize phase")

# cast as an array
t = np.array(t, dtype="f")
q = np.array(q, dtype="f")

####################################################################################################
#
#
#                         Set up itegrator objects
#
#
####################################################################################################

if opts.function == "dxdt_no_NLT_mpi":
  if opts.num_proc == 1: # default to single core functions for speed
    step = odeiv.step_rkf45(dimension, nf.dxdt_no_NLT, args=system)
  else:
    from mpi4py import MPI

    ### set up communicators
    comm = MPI.COMM_SELF.Spawn(sys.executable, args=[__file__.strip("nmode_f.py")+"child-dxdt_no_NLT_mpi.py"], maxprocs=opts.num_proc)
#    print >>outfile, "set up communicators with child processes : ", __file__.strip("nmode_f.py")+"child-dxdt_no_NLT_mpi.py"

    ### broadcast basic system parameters, do this only once
    comm.Bcast([np.array(dimension, dtype="i"), MPI.INT], root=MPI.ROOT) # np.array objects (fast)
#    print >>outfile, "dimension"
    comm.Bcast([np.array(Oorb, dtype="f"), MPI.FLOAT], root=MPI.ROOT)
#    print >>outfile, "Oorb"
    comm.Bcast([np.array(system.network.nlm, dtype="i"), MPI.INT], root=MPI.ROOT)
#    print >>outfile, "nlm"
    comm.bcast(system.network.wyU, root=MPI.ROOT) # python lists (slow), but they have weird shapes...
#    print >>outfile, "wyU"
    comm.bcast(system.network.K, root=MPI.ROOT) # this has different variable types...
#    print >>outfile, "K"

    ### iterate over children and set up persistent communicators
    dxmdt = np.zeros(dimension, dtype="f")
    snds = []
    rcvs = []
    for child in xrange(opts.num_proc):
      snds.append( (comm.Send_init([t, MPI.FLOAT], dest=child), comm.Send_init([q, MPI.FLOAT], dest=child) ) )
      rcvs.append( comm.Recv_init([dxmdt, MPI.FLOAT], source=child) )#, tag=child) )

#    print >>outfile, "succesfully set up persistent communicators"

    step = odeiv.step_rkf45(dimension, nf.dxdt_no_NLT_mpi, args=(dimension, snds, rcvs, dxmdt))

elif opts.function == "dxdt_no_NLT_mp":
  if opts.num_proc == 1:
    step = odeiv.step_rkf45(dimension, nf.dxdt_no_NLT, args=(dimension, system))
  else:
    num_k_par = int(np.ceil(1.0*sum([len(K) for K in system.network.K])/opts.num_proc))
    Msets = []
    n_k = 0
    Mset = []
    for modeNo, K in enumerate(system.network.K):
      Mset.append(modeNo)
      n_k += len(K)
      if n_k >= num_k_par:
        n_k = 0
        Msets.append(Mset)
        Mset = []
    if len(Mset):
      Msets.append(Mset)
    conns = []
    procs = []
    for Mset in Msets:
      con1, con2 = nf.mp.Pipe()
      proc = nf.mp.Process(target=nf.__dxmdt_no_NLT_mp, args=(Mset, system, con2))
      proc.start()
      conns.append(con1)
      procs.append(proc)
    step = odeiv.step_rkf45(dimension, nf.dxdt_no_NLT_mp, args=(dimension, Msets, conns))

elif opts.function == "dxdt_no_NLT_p":
  if opts.num_proc == 1: # default to single core functions for speed
    step = odeiv.step_rkf45(dimension, nf.dxdt_no_NLT, args=(dimension, system))
  else:
    num_k_par = int(np.ceil(1.0*sum([len(K) for K in system.network.K])/opts.num_proc))
    Msets = []
    n_k = 0
    Mset = []
    for modeNo, K in enumerate(system.network.K):
      Mset.append(modeNo)
      n_k += len(K)
      if n_k >= num_k_par:
        n_k = 0
        Msets.append(Mset)
        Mset = []
    if len(Mset):
      Msets.append(Mset)
    step = odeiv.step_rkf45(dimension, nf.dxdt_no_NLT_p, args=(dimension, system, Msets))

elif   opts.function == "dqdt_no_NLT":
  step = odeiv.step_rkf45(dimension, nf.dqdt_no_NLT, args=(dimension, system))

elif opts.function == "dxdt_no_NLT":
  step = odeiv.step_rkf45(dimension, nf.dxdt_no_NLT, args=(dimension, system))

elif opts.function == "dxdt_no_NLT_withPHI":
  step = odeiv.step_rkf45(dimension+1, nf.dxdt_no_NLT_withPHI, args=(dimension, system))

else:
  sys.exit("unkown derivative function: "+ opts.function)

control = odeiv.control_y_new(step, opts.abs_err, opts.rel_err)
if use_phase:
  evolve = odeiv.evolve(step, control, dimension+1)
else:
  evolve = odeiv.evolve(step, control, dimension)

####################################################################################################
#
#
#                                run integration
#
#
####################################################################################################
if opts.time_integration:
  import time
  to=time.time()

# integrate the system and report when necessary 
if use_phase:
  if (not opts.onward) or (opts.init_time != "none"):
    print >>outfile, nm_u.report_func(t, q[:-1], Porb)
  while t < t1:
    while t < sample_time:
      t, h, q = evolve.apply(t, sample_time, h, q) # evolve network and update variables
    sample_time += sample_step
    print >>outfile, nm_u.report_func(t, q[:-1], Porb) # phase is tacked on the end

else:
  if (not opts.onward) or (opts.init_time != "none"):
    print >>outfile, nm_u.report_func(t, q, Porb)
  while t < t1:
    while t < sample_time:
      t, h, q = evolve.apply(t, sample_time, h, q) # evolve network and update variables
    sample_time += sample_step
    print >>outfile, nm_u.report_func(t, q, Porb)
    
if opts.time_integration:
  print >>outfile, "#total integration time = ", time.time()-to, " sec"

####################################################################################################
#
#
#                     clean up any open processes after integration has finished
#
#
####################################################################################################


if opts.function == "dxdt_no_NLT_mp":
  for proc in procs:
    proc.terminate()
    
