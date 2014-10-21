#!python_alias
usage = """ an executable that will produce fits for derived data products """

import numpy as np
import pickle

import fitting
import nmode_utils as nmu
import mode_selection as ms
import gmodes as gm
from nmode_plotting import plt

from optparse import OptionParser

#=================================================
colors = ["b", "r", "g", "m", "c", "y", "k"]
markers = ["o", "s", "^", ">", "<"]

#=================================================
parser = OptionParser(usage=usage)

parser.add_option("-v", "--verbose", default=False, action="store_true")
parser.add_option("-V", "--vverbose", default=False, action="store_true")

parser.add_option("-c", "--clusters", default=[], type="string", action="append", help="filename of pickled lists of clusters")

parser.add_option("-i", "--inset-clusters", default=[], type="string", action="append", help="clusters that will be plotted as a freq sweep in an inset")

parser.add_option("", "--Edot", default=False, action="store_true")
parser.add_option("", "--E", default=False, action="store_true")

parser.add_option("","--unit-system", default="SI", type="string", help="the system of units used in the plot. Currently support either \"SI\" or \"CGS\"")

parser.add_option("-o", "--output-dir", default="./", type="string")
parser.add_option("-t", "--tag", default="", type="string")

parser.add_option("-g", "--grid", default=False, action="store_true")

parser.add_option("", "--raw-data", default=False, action="store_true", help="plot the raw data as well as the fits")
parser.add_option("", "--legend", default=False, action="store_true")

parser.add_option("", "--time-unit", default=False, type="string")

parser.add_option("", "--Mprim", default=1.0, type="float")
parser.add_option("", "--Rprim", default=1.0, type="float")
parser.add_option("", "--n-lin", default=10, type="int", help="number of modes to include in computation of linear dampling rates")
parser.add_option("", "--n-pts", default=1001, type="int", help="the number of points used to construct the linear dampoing rates curve")
parser.add_option("", "--alpha", default=4e-3, type="float")
parser.add_option("", "--c", default=2e-11, type="float")

opts, args = parser.parse_args()

nmu.set_units(system=opts.unit_system) ### set our system of units
energy_unit = nmu.units["energy"]
if opts.time_unit:
	time_unit = opts.time_unit
else:
	time_unit = nmu.units["time"]

if opts.tag:
	opts.tag = "_%s"%opts.tag

#=================================================
if opts.Edot:
	### sweeps_powerlaw fits!
	if opts.verbose:
		print "generating P vs Edot"

	fig = plt.figure()
	ax = plt.subplot(1,1,1)

	### generate a fit for each list of clusters
	for ind, filename in enumerate(opts.clusters):
		color = colors[ind%len(colors)]
		marker = markers[ind%len(markers)]

		label = filename.split("clusters")[-1].strip(".pkl")
		if label[0] == "_":
			label = label[1:]
		label = label.replace("_","\_")
	
		### read in list of clusters
		if opts.verbose: print filename
		file_obj = open(filename, "r")
		clusters = pickle.load( file_obj )
		file_obj.close()

		### load data for each clusters
		for cluster in clusters:
			cluster.load(verbose=opts.vverbose)

		### compute power law fits
		(logp, logedot), (m, b) = fitting.sweeps_Edot_powerlaw(clusters, unit_system=opts.unit_system)

		if opts.verbose:
			print "\tlog(Edot) = %f * log(P) + %f"%(m, b)

		### plot data from ste file directly
		if opts.raw_data:
			mEdot = []
			Porb = []
			for cluster in clusters:
				mEdot += cluster.get_Edot(unit_system=opts.unit_system)[0]
				Porb += cluster.get_Porb(unit_system=opts.unit_system)

			ax.plot(Porb, mEdot, marker="*", markerfacecolor="none", markeredgecolor=color, markersize=4, linestyle="none")

		### plot fitting function
		p = np.exp(logp)
		ax.plot(p, np.exp(logedot), marker=marker, markerfacecolor="none", markeredgecolor=color, markersize=6, linestyle="none")
		ax.plot(p, np.exp(m*logp + b), marker="none", linestyle="-", color=color, label=label)

	### inset clusters (for zoom of a particular region)
	if opts.inset_clusters:
		iax = fig.add_axes([0.55, 0.55, 0.40, 0.40])
		### compute linear expectation
	        wo = ms.compute_wo(opts.Mprim, opts.Rprim)

		if opts.verbose: print "determining frequency range"
		minO = np.infty
		maxO = -np.infty
		Mcomp = False
		for filename in opts.inset_clusters:
			file_obj = open(filename, "r")
                        clusters = pickle.load( file_obj )
                        file_obj.close()

                        for cluster in clusters:
                                cluster.load(verbose=opts.vverbose)
				for sdata, _ in cluster.data:
					O = 2*np.pi / sdata["system"]["Porb"]
					if minO > O:
						minO = O
					if maxO < O:
						maxO = O
					if Mcomp and Mcomp != sdata["system"]["Mcomp/Mjup"]:
						raise ValueError, "Mcomp mismatch"
					else:
						Mcomp = sdata["system"]["Mcomp/Mjup"]
			del clusters

		Edot_lin = []
		resonances = set()
		O = np.linspace(minO, maxO, opts.n_pts)
	        l = 2 ### dominant harmonic of the tide
	        m = 2
	        ecc = 0.0 ### take zero eccentricity limit
	        for o in O:
	                n_res = int(round(opts.alpha*l/(m*o), 0))
	                edot_lin = 0.0
	                for n in range(max(1, n_res-opts.n_lin), n_res+opts.n_lin):
	                        mode = gm.gmode(n=n, l=l, m=m, alpha=opts.alpha, c=opts.c, wo=wo)
	                        mode.compute_forcing(2*np.pi/o, ecc, opts.Mprim, Mcomp, opts.Rprim)
	                        w, y, U = mode.get_wyU()
	                        resonances.add( w/m )
	                        for u, k in U:
	                                edot_lin += 2*y*(w*u)**2 / ( (abs(w)-m*k*o)**2 + y**2 ) ### add linear result
	                Edot_lin.append( edot_lin )

		nmu.set_units(opts.unit_system)
		Eo = nmu.G*(opts.Mprim*nmu.Msun)**2 / (opts.Rprim*nmu.Rsun)

		Edot_lin = np.array(Edot_lin)*Eo

	        ### plot linear result
	        #ax.semilogy(O, Edot_lin, color='k', alpha=0.5)
	        iax.plot((2*np.pi/O), Edot_lin, color='k', alpha=0.5)

		for ind, filename in enumerate(opts.inset_clusters):
			color = colors[ind%len(colors)]
			marker = markers[ind%len(markers)]

			if opts.verbose: print filename
			file_obj = open(filename, "r")
			clusters = pickle.load( file_obj )
			file_obj.close()

			for cluster in clusters:
				cluster.load(verbose=opts.vverbose)

			mEdot = []
                        Porb = []
                        for cluster in clusters:
                                mEdot += cluster.get_Edot(unit_system=opts.unit_system)[0]
                                Porb += cluster.get_Porb(unit_system=opts.unit_system)

                        iax.plot(Porb, mEdot, marker=marker, markerfacecolor=color, markeredgecolor=color, markersize=4, linestyle="none")

		iax.set_xlabel("$P_\mathrm{orb}$ [%s]"%time_unit)
		iax.set_ylabel("$\left< \partial_t E_\mathrm{orb} \\right>$ [%s/%s]"%(energy_unit, time_unit))
		iax.grid(opts.grid, which="both")

		iax.set_yscale('log')

	### decoration, etc
	ax.set_xlabel("$P_\mathrm{orb}$ [%s]"%time_unit)
	ax.set_ylabel("$\left< \partial_t E_\mathrm{orb} \\right>$ [%s/%s]"%(energy_unit, time_unit))
	ax.grid(opts.grid, which="both")

	ax.set_xscale('log')
	ax.set_yscale('log')

	if opts.legend:
		ax.legend(loc="lower left")

	### save, etc
	figname = "%s/sweeps_Edot-powerlaw%s.png"%(opts.output_dir, opts.tag)
	if opts.verbose:
		print "saving ", figname
	fig.savefig(figname)
	plt.close(fig)

#================================================
if opts.E:
        ### sweeps_powerlaw fits!
        if opts.verbose:
                print "generating P vs E"

        fig = plt.figure()
        ax = plt.subplot(1,1,1)

        ### generate a fit for each list of clusters
        for ind, filename in enumerate(opts.clusters):
                color = colors[ind%len(colors)]
                marker = markers[ind%len(markers)]

                label = filename.split("clusters")[-1].strip(".pkl")
                if label[0] == "_":
                        label = label[1:]
                label = label.replace("_","\_")

                ### read in list of clusters
                if opts.verbose: print filename
                file_obj = open(filename, "r")
                clusters = pickle.load( file_obj )
                file_obj.close()

                ### load data for each clusters
                for cluster in clusters:
                        cluster.load(verbose=opts.vverbose)

                ### compute power law fits
                (logp, loge), (m, b) = fitting.sweeps_E_powerlaw(clusters, unit_system=opts.unit_system)

                if opts.verbose:
                        print "\tlog(E) = %f * log(P) + %f"%(m, b)

                ### plot data from ste file directly
                if opts.raw_data:
                        mE = []
                        Porb = []
                        for cluster in clusters:
                                mE += cluster.get_E(unit_system=opts.unit_system)[0]
                                Porb += cluster.get_Porb(unit_system=opts.unit_system)

                        ax.plot(Porb, mE, marker="*", markerfacecolor="none", markeredgecolor=color, markersize=4, linestyle="none")

                ### plot fitting function
                p = np.exp(logp)
                ax.plot(p, np.exp(loge), marker=marker, markerfacecolor="none", markeredgecolor=color, markersize=6, linestyle="none")
                ax.plot(p, np.exp(m*logp + b), marker="none", linestyle="-", color=color, label=label)

        ### inset clusters (for zoom of a particular region)
        if opts.inset_clusters:
                iax = fig.add_axes([0.55, 0.55, 0.40, 0.40])
                ### compute linear expectation
                wo = ms.compute_wo(opts.Mprim, opts.Rprim)

                if opts.verbose: print "determining frequency range"
                minO = np.infty
                maxO = -np.infty
                Mcomp = False
                for filename in opts.inset_clusters:
                        file_obj = open(filename, "r")
                        clusters = pickle.load( file_obj )
                        file_obj.close()

                        for cluster in clusters:
                                cluster.load(verbose=opts.vverbose)
                                for sdata, _ in cluster.data:
                                        O = 2*np.pi / sdata["system"]["Porb"]
                                        if minO > O:
                                                minO = O
                                        if maxO < O:
                                                maxO = O
                                        if Mcomp and Mcomp != sdata["system"]["Mcomp/Mjup"]:
                                                raise ValueError, "Mcomp mismatch"
                                        else:
                                                Mcomp = sdata["system"]["Mcomp/Mjup"]
                        del clusters

                E_lin = []
                resonances = set()
                O = np.linspace(minO, maxO, opts.n_pts)
                l = 2 ### dominant harmonic of the tide
                m = 2
                ecc = 0.0 ### take zero eccentricity limit
                for o in O:
                        n_res = int(round(opts.alpha*l/(m*o), 0))
                        e_lin = 0.0
                        for n in range(max(1, n_res-opts.n_lin), n_res+opts.n_lin):
                                mode = gm.gmode(n=n, l=l, m=m, alpha=opts.alpha, c=opts.c, wo=wo)
                                mode.compute_forcing(2*np.pi/o, ecc, opts.Mprim, Mcomp, opts.Rprim)
                                w, y, U = mode.get_wyU()
                                resonances.add( w/m )
                                for u, k in U:
                                        e_lin += (w*u)**2 / ( (abs(w)-m*k*o)**2 + y**2 ) ### add linear result
                        E_lin.append( e_lin )

                nmu.set_units(opts.unit_system)
                Eo = nmu.G*(opts.Mprim*nmu.Msun)**2 / (opts.Rprim*nmu.Rsun)

                E_lin = np.array(E_lin)*Eo

                ### plot linear result
                #ax.semilogy(O, Edot_lin, color='k', alpha=0.5)
                iax.plot((2*np.pi/O), E_lin, color='k', alpha=0.5)

                for ind, filename in enumerate(opts.inset_clusters):
                        color = colors[ind%len(colors)]
                        marker = markers[ind%len(markers)]

                        if opts.verbose: print filename
                        file_obj = open(filename, "r")
                        clusters = pickle.load( file_obj )
                        file_obj.close()

                        for cluster in clusters:
                                cluster.load(verbose=opts.vverbose)

                        mE = []
                        Porb = []
                        for cluster in clusters:
                                mE += cluster.get_E(unit_system=opts.unit_system)[0]
                                Porb += cluster.get_Porb(unit_system=opts.unit_system)

                        iax.plot(Porb, mE, marker=marker, markerfacecolor=color, markeredgecolor=color, markersize=4, linestyle="none")

                iax.set_xlabel("$P_\mathrm{orb}$ [%s]"%time_unit)
                iax.set_ylabel("$\left< E_{\\ast} \\right>$ [%s]"%(energy_unit))
                iax.grid(opts.grid, which="both")

                iax.set_yscale('log')

        ### decoration, etc
        ax.set_xlabel("$P_\mathrm{orb}$ [%s]"%time_unit)
        ax.set_ylabel("$\left< E_{\\ast} \\right>$ [%s]"%(energy_unit))
        ax.grid(True, which="both")

        ax.set_xscale('log')
        ax.set_yscale('log')

        if opts.legend:
                ax.legend(loc="lower left")

        ### save, etc
        figname = "%s/sweeps_E-powerlaw%s.png"%(opts.output_dir, opts.tag)
        if opts.verbose:
                print "saving ", figname
        fig.savefig(figname)
        plt.close(fig)


#================================================
### other fits?
### fits using Mcomp as well?
### do this more by hand once we have Porb sweeps fits?



