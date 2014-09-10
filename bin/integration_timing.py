#!python_alias
usage="""generates integration timing data from *out files"""

import nmode_utils as nmu
import nmode_plotting as nmp

import numpy as np

from optparse import OptionParser

#=================================================

parser=OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("-F", "--outfilename", type="string")

parser.add_option("", "--tmin", default=0.0, type="float")
parser.add_option("", "--tmax", default=np.infty, type="float")
parser.add_option("", "--downsample", default=1, type="float")

parser.add_option("", "--hist", default=False, action="store_true")
parser.add_option("", "--n-per-bin", default=10, type="float")
parser.add_option("", "--log-bin", default=False, action="store_true")

parser.add_option("", "--dt-time", default=False, action="store_true")

parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("-g", "--grid", default=False, action="store_true")

opts, args = parser.parse_args()

if opts.tag:
	opts.tag = "."+opts.tag

if not (opts.hist or opts.dt_time):
  raise StandardError, "nothing to do..."

#=================================================
### read in outfiles, get dt
if opts.verbose: print "reading data from ", opts.outfilename
t_P, q, N_m, dt = nmu.load_out(opts.outfilename, tmin=opts.tmin, tmax=opts.tmax, downsample=opts.downsample, timing=True)

if not len(dt):
	raise StandardError, "could not find any integration timing information in ", opts.outfilename

#=================================================
### generate plots
### histogram dt
if opts.hist:
	if opts.verbose: print "\thist"
	fig = nmp.plt.figure()
	ax = nmp.plt.subplot(1,1,1)

	num_bin = len(dt)/opts.n_per_bin+1

	if opts.log_bin:
		bins = np.logspace(np.log10(np.min(dt)/1.1), np.log10(np.max(dt)*1.1), num_bin+1)
	else:
		bins = np.linspace(np.min(dt)/1.1, np.max(dt)*1.1, num_bin+1)

	n, b, p = ax.hist(dt, bins=bins, histtype="step")

	if opts.log_bin:
		ax.set_xscale("log")

	ax.set_xlabel("duration [sec]")
	ax.set_ylabel("count")

	ax.set_ylim(ymin=min(n)/1.1, ymax=max(n)*1.1)

	ax.grid(opts.grid, which="both")

	figname = opts.outfilename+".dt-hist%s.png"%opts.tag
	if opts.verbose: print figname
	fig.savefig(figname)
	nmp.plt.close(fig)

### plot dt vs time
if opts.dt_time:
	if opts.verbose: print "\tdt_time"
	
	### check lengths of dt and t_P
	npts = len(t_P)
	ndt = len(dt)
	if ndt != npts:
		if ndt == npts-1:
			t_P=t_P[1:]
		else:
			raise ValueError, "len(t_P) = %d != %d = len(dt)"%(npts, ndt)

	fig = nmp.plt.figure()
	ax = nmp.plt.subplot(1,1,1)

	ax.plot(t_P, dt, label="duration")

	ax.set_xlabel("$t/P_{\mathrm{orb}}$")
	ax.set_ylabel("duration [sec]")

	ax.grid(opts.grid, which="both")

	figname = opts.outfilename+".dt.lin_time%s.png"%opts.tag
	if opts.verbose: print figname
	fig.savefig(figname)
	nmp.plt.close(fig)






