usage="""written to attempt to decompose networks meaningfully"""

import nmode_utils as nmu
import prune_network as pn
import nmode_plotting as nmp
import nmode_diagnostic as nmd
import nmode_state as nms

import numpy as np
import pickle

#==================================================
### define utilities
#==================================================
def plot(t_P, q, system, gens, mode_nums=None, verbose=False, tc="x", xmin=False, xmax=False, ymin=False, ymax=False):
	""" a wrapper for plotting routines we'll use a lot """

	n_l_m = system.network.nlm

	ans = []

	### amp-time
	if verbose: print "amp_time"
	fig, ax = nmp.amp_plot(t_P, q, n_l_m=n_l_m, mode_nums=mode_nums)
	ax.set_ylabel(r'$|'+tc+'_i|$')
	ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')
	if ymin:
		ax.set_ylim(ymin=ymin)
	if ymax:
		ax.set_ylim(ymax=ymax)
	if xmin:
		ax.set_xlim(xmin=xmin)
	if xmax:
		ax.set_xlim(xmax=xmax)
	ans.append( (fig, ax) )

	### multi-gen-amp-time
	if verbose: print "multi-gen_amp_time"
	fig, axs = nmp.multi_gen_amp_plot(t_P, q, gens, n_l_m=n_l_m, mode_nums=mode_nums)
	for ax in axs:
		ax.set_ylabel(r'$|'+tc+'|$')
	ax.set_xlabel(r'$t/P_{\mathrm{orb}}$')
	if ymin:
		ax.set_ylim(ymin=ymin)
	if ymax:
		ax.set_ylim(ymax=ymax)
	if xmin:
		ax.set_xlim(xmin=xmin)
	if xmax:
		ax.set_xlim(xmax=xmax)
	ans.append( (fig, axs) )

	### coupling diagrams	
	colormap = nmd.plt.get_cmap("jet")

	### nl-placement disp
	if verbose: print "nl-placmenet disp"
	myE = -np.array([np.mean(_) for _ in nms.viscous_disp(q, network, Eo=1.0)[-1]])
	mode_colors = [colormap(_) for _ in myE]
	mode_order = myE.argsort()
	fig = nmd.plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.750, 0.8])
	cbax = fig.add_axes([0.875, 0.1, 0.025, 0.8])
	fig, ax = nmd.coupling_tree_nl(system, tree_type="placement", genNos=[], verbose=verbose, mode_colors=mode_colors, mode_order=mode_order, fig_ax=(fig,ax), mode_nums=mode_nums)
#	colorbar = nmd.matplotlib.colorbar.ColorbarBase(cbax, cmap=colormap, orientation='vertical', norm=nmd.matplotlib.colors.LogNorm())
	colorbar = nmd.matplotlib.colorbar.ColorbarBase(cbax, cmap=colormap, orientation='vertical')
	colorbar.ax.set_ylabel(r"$\gamma_i A_i^2/\mathrm{max}\{\gamma_i A_i^2\}$")
	ans.append( (fig, ax, colorbar) )

	### nl-placement E
	if verbose: print "nl-placement E"
	mE = np.array([np.mean(_) for _ in nms.compute_E(q, Eo=1.0)])
	mode_colors = [colormap(_) for _ in mE/max(mE)]
	mode_order = mE.argsort()

	fig = nmd.plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.750, 0.8])
	cbax = fig.add_axes([0.875, 0.1, 0.025, 0.8])
	fig, ax = nmd.coupling_tree_nl(system, tree_type="placement", genNos=[], verbose=verbose, mode_colors=mode_colors, mode_order=mode_order, fig_ax=(fig,ax), mode_nums=mode_nums)
#	colorbar = nmd.matplotlib.colorbar.ColorbarBase(cbax, cmap=colormap, orientation='vertical', norm=nmd.matplotlib.colors.LogNorm())
	colorbar = nmd.matplotlib.colorbar.ColorbarBase(cbax, cmap=colormap, orientation='vertical')
	colorbar.ax.set_ylabel(r"$A_i^2/\mathrm{max}\{A_i^2\}$")
	ans.append( (fig, ax, colorbar) )

	### nl-siblings gens
	for genNo in xrange(1, len(gens)):
		if verbose: print "nl-siblings gen%d"%genNo
		fig, ax = nmd.coupling_tree_nl(system, tree_type="siblings", genNos=[genNo], verbose=verbose, mode_nums=mode_nums)
		ans.append( (fig, ax) )

	return ans

def savefig(ans, outfilename, logfilename, tag="", no_legend=False, grid=False, verbose=False):
	""" computes fignames and saves figures within ans """

	if tag:
		tag = "."+tag

	### amp_time
	fig, ax = ans.pop(0)
	if not no_legend:
		ax.legend()
	ax.grid(grid, which="both")
	figname = "%s.Lamp.lin_time%s.png"%(outfilename, tag)
	if verbose: print figname
	fig.savefig(figname)

	ax.set_yscale('log')
	figname = "%s.amp.lin_time%s.png"%(outfilename, tag)
	if verbose: print figname
	fig.savefig(figname)
	nmp.plt.close(fig)

	### multi-gen_amp_time
	fig, axs = ans.pop(0)
	for ax in axs:
		if not no_legend:
			ax.legend()
		ax.grid(grid, which="both")
	figname = "%s.multi-gen_Lamp.lin_time%s.png"%(outfilename, tag)
	if verbose: print figname
	fig.savefig(figname)

	for ax in axs:
		ax.set_yscale('log')
	figname = "%s.multi-gen_amp.lin_time%s.png"%(outfilename, tag)
        if verbose: print figname
        fig.savefig(figname)
	nmp.plt.close(fig)

	### nl_placement-disp
	fig, ax, colorbar = ans.pop(0)
	ax.grid(grid, which="both")
	figname = "%s.nl_placement-disp%s.png"%(logfilename, tag)
	if verbose: print figname
	fig.savefig(figname)
	nmd.plt.close(fig)

	### nl_placmeent-disp
	fig, ax, colorbar = ans.pop(0)
	ax.grid(grid, which="both")
	figname = "%s.nl_placement-E%s.png"%(logfilename, tag)
	if verbose: print figname
	fig.savefig(figname)
	nmd.plt.close(fig)

	### nl_siblings gens
	genNo = 1
	while len(ans):
		fig, ax = ans.pop(0)
		ax.grid(grid, which="both")
		figname = "%s.nl_siblings-gen%d%s.png"%(logfilename, genNo, tag)
		if verbose: print figname
		fig.savefig(figname)
		nmd.plt.close(fig)
		genNo += 1 

#=================================================

if __name__ == "__main__":

	from optparse import OptionParser

	parser = OptionParser(usage=usage)

	parser.add_option("-v", "--verbose", default=False, action="store_true")
	parser.add_option("", "--vverbose", default=False, action="store_true")

	parser.add_option("-F", "--outfilename", default=None, type="string")
	parser.add_option("", "--tcurrent", dest="tc", default=None, type="string")

	parser.add_option("-l", "--logfilename", default=None, type="string")

	parser.add_option("-p", "--pklfilename", default=None, type="string", help="filename into which summary information will be pickled")

	parser.add_option("-P", "--plot", default=False, action="store_true")

	parser.add_option("", "--time-ymin", default=False, type="float", help="the lower limit on the y-axis in time domain figures.")
	parser.add_option("", "--time-ymax", default=False, type="float", help="the upper limit on the y-axis in time domain figures.")
	parser.add_option("", "--time-xmin", default=False, type="float", help="the lower limit on the x-axis in time domain figures.")
	parser.add_option("", "--time-xmax", default=False, type="float", help="the upper limit on the x-axis in time domain figures.")

	parser.add_option("-d", "--downsample", default=False, type="int", help="plot only 1 point out of this many")

	parser.add_option("", "--no-legend", default=False, action="store_true")
	parser.add_option("-g", "--grid", default=False, action='store_true')

	parser.add_option("-t", "--tag", default="", type="string")

	parser.add_option("-k", "--num-k", default=[], action="append", type="float")
	parser.add_option("-e", "--Ethr", default=[], action="append", type="float")
	parser.add_option("-c", "--collE", default=[], action="append", type="float")

	opts, args = parser.parse_args()
	
	if not opts.outfilename:
		opts.outfilename = raw_input("outfilename = ")
	if not opts.tc:
		opts.tc = raw_input("tcurrent = ")
	if not opts.logfilename:
		opts.logfilename = raw_input("logfilename = ")

	if opts.tag:
		opts.tag = "."+opts.tag

#	if not opts.num_k:
#		opts.num_k = np.arange(1, 25, 1)
#	if not opts.Ethr:
#		opts.Ethr = np.logspace(-21, -15, 100)
#	if not opts.collE:
#		opts.collE = np.logspace(-21, -15, 100)

	#=================================================
	### load data
	#=================================================

	if opts.verbose: print "loading data from", opts.outfilename
	t_P, q, N_m = nmu.load_out(opts.outfilename, tmin=opts.time_xmin, tmax=opts.time_xmax, downsample=opts.downsample)

	if opts.verbose: print "loading system from", opts.logfilename
	system = nmu.load_log(opts.logfilename)
	network = system.network

	if opts.verbose: print "decomposing network into generations"
	gens, coups = network.gens()

	#=================================================
	### define subsets of modes
	#=================================================

	### pull apart by number of couplings
	if opts.verbose: print "separating modes by num_k:"
	num_k_modeNos = {}
	for num_k in opts.num_k:
		if opts.verbose: print "\tnum_k:", num_k
		num_k_modeNos[num_k] = [[network.modeNoD[nlms] for nlms in _] for _ in pn.downselect_num_k(num_k, network, mode_nums=None)]

	### pull apart by Ethr
	if opts.verbose: print "separating modes by Ethr:"
	ethr_modeNos = {}
	for ethr in opts.Ethr:
		if opts.verbose: print "\tEthr:", ethr
		ethr_modeNos[ethr] = [[network.modeNoD[nlms] for nlms in _] for _ in pn.downselect_Ethr(ethr, system, freqs=None, mode_nums=None)]

	### pull apart by collE
	if opts.verbose: print "separating modes by collective Ethr:"
	collE_modeNos = {}
	for collE in opts.collE:
		if opts.verbose: print "\tcollE:", collE
		collE_modeNos[collE] = [[network.modeNoD[nlms] for nlms in _] for _ in pn.downselect_collE(collE, system, freqs=None, mode_nums=None)]


	#=================================================
	### plot data!
	#=================================================

	### compute stuff for statistics
	myE = -np.array([np.mean(_) for _ in nms.viscous_disp(q, network, Eo=1.0)[-1]])
	myE /= sum(myE)
	mE = np.array([np.mean(_) for _ in nms.compute_E(q, Eo=1.0)])
	mE /= sum(mE)

	pkl_obj = {'myE':myE, 'mE':mE}

	#=== total data set
	if opts.plot:
		if opts.verbose: print "total data"
		tag = ""+opts.tag
		savefig( plot(t_P, q, system, gens, mode_nums=None, verbose=opts.vverbose, tc=opts.tc), opts.outfilename, opts.logfilename, tag=tag, no_legend=opts.no_legend, grid=opts.grid, verbose=opts.vverbose)

	#=== num_k
	if opts.verbose: print "num_k"
	pkl_obj['num_k'] = {}
	for num_k in sorted(num_k_modeNos.keys()):
		g, l = num_k_modeNos[num_k]
                pkl_obj['num_k'][num_k] = [(g, sum(mE[g]), sum(myE[g])), (l, sum(mE[l]), sum(myE[l]))]
		if opts.verbose: print "\tnum_k=",num_k,"\n\tN>=",len(g),"\n\t\tmE=",sum(mE[g]),"\n\t\tmyE=",sum(myE[g]),"\n\tN<=",len(l),"\n\t\tmE=",sum(mE[l]),"\n\t\tmyE=",sum(myE[l])

		if opts.plot:
			tag = "num_k-g-%d"%num_k+opts.tag
			if opts.vverbose: print tag
			savefig( plot(t_P, q, system, gens, mode_nums=g, verbose=opts.vverbose, tc=opts.tc, xmin=opts.time_xmin, xmax=opts.time_xmax, ymin=opts.time_ymin, ymax=opts.time_ymax), opts.outfilename, opts.logfilename, tag=tag, no_legend=opts.no_legend, grid=opts.grid, verbose=opts.vverbose)

			tag = "num_k-l-%d"%num_k+opts.tag
			if opts.vverbose: print tag
			savefig( plot(t_P, q, system, gens, mode_nums=l, verbose=opts.vverbose, tc=opts.tc, xmin=opts.time_xmin, xmax=opts.time_xmax, ymin=opts.time_ymin, ymax=opts.time_ymax), opts.outfilename, opts.logfilename, tag=tag, no_legend=opts.no_legend, grid=opts.grid, verbose=opts.vverbose)

	#=== Ethr
	if opts.verbose: print "Ethr"
	pkl_obj['Ethr'] = {}
	for ethr in sorted(ethr_modeNos.keys()):
		g, l = ethr_modeNos[ethr]
		pkl_obj['Ethr'][ethr] = [(g, sum(mE[g]), sum(myE[g])), (l, sum(mE[l]), sum(myE[l]))]
		if opts.verbose: print "\tEthr=",ethr,"\n\tN>=",len(g),"\n\t\tmE=",sum(mE[g]),"\n\t\tmyE=",sum(myE[g]),"\n\tN<=",len(l),"\n\t\tmE=",sum(mE[l]),"\n\t\tmyE=",sum(myE[l])

		if opts.plot:
			tag = "Ethr-g-%.3fe%d"%nmu.float_to_scientific(ethr)+opts.tag
			if opts.vverbose: print tag
			savefig( plot(t_P, q, system, gens, mode_nums=g, verbose=opts.vverbose, tc=opts.tc, xmin=opts.time_xmin, xmax=opts.time_xmax, ymin=opts.time_ymin, ymax=opts.time_ymax), opts.outfilename, opts.logfilename, tag=tag, no_legend=opts.no_legend, grid=opts.grid, verbose=opts.vverbose)

			tag = "Ethr-l-%.3fe%d"%nmu.float_to_scientific(ethr)+opts.tag
			if opts.vverbose: print tag
			savefig( plot(t_P, q, system, gens, mode_nums=l, verbose=opts.vverbose, tc=opts.tc, xmin=opts.time_xmin, xmax=opts.time_xmax, ymin=opts.time_ymin, ymax=opts.time_ymax), opts.outfilename, opts.logfilename, tag=tag, no_legend=opts.no_legend, grid=opts.grid, verbose=opts.vverbose)

	#=== collE
	if opts.verbose: print "collE"
	pkl_obj['collE'] = {}
	for collE in sorted(collE_modeNos.keys()):
		g, l = collE_modeNos[collE]
		pkl_obj['collE'][collE] = [(g, sum(mE[g]), sum(myE[g])), (l, sum(mE[l]), sum(myE[l]))]
		if opts.verbose: print "\tcollE=",collE,"\n\tN>=",len(g),"\n\t\tmE=",sum(mE[g]),"\n\t\tmyE=",sum(myE[g]),"\n\tN<=",len(l),"\n\t\tmE=",sum(mE[l]),"\n\t\tmyE=",sum(myE[l])

		if opts.plot:
			tag = "collE-g-%.3fe%d"%nmu.float_to_scientific(collE)+opts.tag
			if opts.vverbose: print tag
			savefig( plot(t_P, q, system, gens, mode_nums=g, verbose=opts.vverbose, tc=opts.tc, xmin=opts.time_xmin, xmax=opts.time_xmax, ymin=opts.time_ymin, ymax=opts.time_ymax), opts.outfilename, opts.logfilename, tag=tag, no_legend=opts.no_legend, grid=opts.grid, verbose=opts.vverbose)

			tag = "collE-l-%.3fe%d"%nmu.float_to_scientific(collE)+opts.tag
			if opts.vverbose: print tag
			savefig( plot(t_P, q, system, gens, mode_nums=l, verbose=opts.vverbose, tc=opts.tc, xmin=opts.time_xmin, xmax=opts.time_xmax, ymin=opts.time_ymin, ymax=opts.time_ymax), opts.outfilename, opts.logfilename, tag=tag, no_legend=opts.no_legend, grid=opts.grid, verbose=opts.vverbose)

	#================================================
	### pkl data
	#================================================
	if opts.pklfilename:
		if opts.verbose: print "writing data to", opts.pklfilename
		f=open(opts.pklfilename, "w")
		pickle.dump(pkl_obj, f)
		f.close()

	#================================================
	### scatters
	#================================================
	if opts.plot:
		ax_pos = [0.10, 0.15, 0.8, 0.8]
		diag = np.linspace(0,1,101)
		for key in ['Ethr', 'collE', 'num_k']:
			fig_g = nmp.plt.figure()
			ax_g_me = fig_g.add_axes(ax_pos)
			ax_g_mye = ax_g_me.twinx()

			fig_l = nmp.plt.figure()
			ax_l_me = fig_l.add_axes(ax_pos)
			ax_l_mye = ax_l_me.twinx()

			for _k, ((g, g_me, g_mye), (l, l_me, l_mye)) in pkl_obj[key].items():
				N_g = 1.0*len(g)/N_m
				ax_g_me.plot(N_g, g_me, marker="o", markerfacecolor='none', markeredgecolor='b', markersize=8)				
				ax_g_mye.plot(N_g, g_mye, marker="s", markerfacecolor='none', markeredgecolor='g', markersize=8)				

				N_l = 1.0*len(l)/N_m				
				ax_l_me.plot(N_l, l_me, marker="o", markerfacecolor='none', markeredgecolor='b', markersize=8)				
				ax_l_mye.plot(N_l, l_mye, marker="s", markerfacecolor='none', markeredgecolor='g', markersize=8)				
				
			ax_g_me.set_xlabel('fraction of modes')
			ax_g_me.set_ylabel('fraction of Energy', color="b")
			ax_g_mye.set_ylabel('fraction of dissapation', color="g")
			ax_g_mye.yaxis.tick_right()

			ax_g_me.grid(True, which="both")

			ax_g_me.set_ylim(ymin=0, ymax=1)
			ax_g_mye.set_ylim(ymin=0, ymax=1)
	
			ax_g_me.plot(diag, diag, 'k--')

			ax_l_me.set_xlabel('fraction of modes')
			ax_l_me.set_ylabel('fraction of Energy', color="b")
			ax_l_mye.set_ylabel('fraction of dissapation', color="g")
			ax_l_mye.yaxis.tick_right()

			ax_l_me.grid(True, which="both")

			ax_l_me.set_ylim(ymin=0, ymax=1)
			ax_l_mye.set_ylim(ymin=0, ymax=1)

			ax_l_me.plot(diag, diag, 'k--')

			fig_g.text(0.10,0.95, "large: "+key, ha='left', va='top')
			fig_g.text(0.90,0.15, "small: "+key, ha='right', va='bottom')

			figname_g = key+"_g.png"
			if opts.verbose: print figname_g
			fig_g.savefig(figname_g)
			nmp.plt.close(fig_g)

			fig_l.text(0.10,0.95, "small: "+key, ha='left', va='top')
			fig_l.text(0.90,0.15, "large: "+key, ha='right', va='bottom')

			figname_l = key+"_l.png"
			if opts.verbose: print figname_l
			fig_l.savefig(figname_l)
			nmp.plt.close(fig_l)









