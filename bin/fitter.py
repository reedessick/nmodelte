#!python_alias
usage = """ an executable that will produce fits for derived data products """

import fitting

import nmode_utils as nmu
from nmode_plotting import plt

from optparse import OptionParser

#=================================================
parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")

parser.add_option("", "--sweeps_powerlaw", default=[], type="string", action="append", help="stepkl filenames that will be used for sweeps_powerlaw fitting")
parser.add_option("", "--Porb_window", default=1000, type="float", help="the clustering window for harmonic averages during sweeps fitting")

parser.add_option("-o", "--output-dir", default="./", type="string")
parser.add_option("-t", "--tag", default="", type="string")

opts, args = parser.parse_args()

if opts.tag:
	opts.tag = "_%s"%opts.tag

#=================================================
### sweeps_powerlaw fits!
if opts.verbose:
	print "generating sweeps_powerlaw fits"

(logp, logedot), (m, b) = fitting.sweeps_powerlaw(opts.sweeps_powerlaw, Porb_window=opts.Porb_window)

if opts.verbose:
	print "\tlog(Edot) = %f * log(P) + %f"%(m, b)
	print "plotting data..."

fig = plt.figure()
ax = plt.subplot(1,1,1)

### plot data from ste file directly
mEdot = []
sEdot = []
Porb = []
for filename in opts.sweeps_powerlaw:
	sdata, mdata = nmu.load_ste(filename)

	Eorb = abs(sdata["system"]["Eorb"])
	porb = sdata["system"]["Porb"]
	Porb.append( porb )
	mEdot.append( sdata["stats"]["mean{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/porb )
	sEdot.append( sdata["stats"]["stdv{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/porb )

ax.plot(Porb, mEdot, marker="*", markerfacecolor="none", markeredgecolor="b", linestyle="none")

### plot fitting function
p = np.exp(logp)
ax.plot(p, np.exp(logedot), marker="o", markerfacecolor="none", markeredgecolor="r", linestyle="none")
ax.plot(p, m*p + b, marker="none", linestyle="-", color="g")

### decoration, etc
ax.set_xlabel("$P_\mathrm{orb}$ [sec]")
ax.set_ylabel("$\left< \partial_t E_\mathrm{orb} \\right>$ [energy units???/sec]")
ax.grid(True, which="both")

fig.figtext(0.9, 0.9, "$\log \left< \partial_t E_\mathrm{orb} \\right> = %f * \log P_\mathrm{orb} + %f"%(m, b), ha="right", va="top")

### save, etc
figname = "%s/sweeps_powerlaw%s.png"%(opts.output_dir, opts.tag)
if opts.verbose:
	print "saveing ", figname
fig.savefig(figname)
plt.close(fig)

#================================================

### other fits?
### fits using Mcomp as well?
### do this more by hand once we have Porb sweeps fits?



