usage="""utilities for fitting derived data products"""

import numpy as np
import pickle
import nmode_utils as nmu


#=================================================
# sweeps fitting
#=================================================
def sweeps_powerlaw(stepkl_filenames, Porb_window=1000):
	"""
	a function that loads state data from stepkl_filenames and attempts to fit the data with a standard functional
	data are first clustered with Porb_window 
		we compute harmonic averages over the clustered data to compute the average dissipation rate (time-weighted)
	We then fit the resulting data points with a function
		laurant series?
		power law?

	WARNING: sorts data by Porb and assumes there is a unique system for each Porb. Will overwrite data if fed more than one system for any Porb.
	ALSO ASSUMES all system data is identical EXCEPT for Porb so the fit makes sense.
	"""

	### load data from files
	data = {}
	for filename in ste_pkl_filenames:
		sdata, mdata = nmu.load_ste(filename)
		data[sdata["system"]["Porb"]] = (sdata, mdata)

	### cluster data by Porb_window
	Porbs = sorted(data.keys())
	clusters = []
	cluster = [Porbs[0]]
	for Porb in Porbs[1:]:
		if Porb-old_Porb < Porb_window: ### we're in the same cluster
			cluster.append( data[Porb] )
		else:
			clusters.append( cluster )
			cluster = [data[Porb]]

	### compute harmonic average (via delegation)
	p = np.array([np.mean([sdata["system"]["Porb"] for sdata, _ in cluster] for cluster in clusters])
	edot = np.array( [harmonic_average(cluster)) for cluster in clusters] )

	### fit to function!
	### simple power law!
	logp = np.log(p)
	logedot = np.log(edot)

	### linear fit in log-space
	yx = np.sum( logp*logedot )
	y = np.sum( logedot )
	xx = np.sum( logedot*logedot )
	x = np.sum( logedot )
	n = np.len( logedot )

	det = xx*n - x*x
	m = (n*yx - x*y)/det
	b = (-x*yx + xx*y)/det

	return (logp, logedot), (m, b)

###
def harmonic_average(cluster):
	"""
	computes the harmonic average over the data within cluster
        we compute <Edot> = int dt Edot / int dt
        but we take the integral over Porb, so changing the integration variable to Porb yields

        <Edot> = int dP ( E/P * Edot**-1) * Edot / int dP ( E/P * Edot**-1 )

	ASSUMES spacing is even in Porb-space (within each cluster!), so the integral's measure drops out
	"""
	num = 0.0
	den = 0.0
	for sdata, mdata in cluster:
		Porb = sdata["system"]["Porb"]
		Eorb = sdata["system"]["Eorb"]
		Edot = sdata["stats"]["mean{|sum{Edot}|*(Porb/|Eorb|)}"]
		
		num += Eorb/Porb
		den += 1.0/Edot

	return num/den

