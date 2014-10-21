usage="""utilities for fitting derived data products"""

import numpy as np
import pickle
import nmode_utils as nmu

#=================================================
# cluster object
#=================================================
class Cluster(object):
	def __init__(self, filenames):
		self.filenames = filenames
		self.data = None

	def load(self, verbose=False):
		self.data = []
		for filename in self.filenames:
			if verbose:
				print filename
			self.data.append( nmu.load_ste(filename) )

	def get_Edot(self, unit_system="SI"):
		nmu.set_units(system=unit_system)
		energy_unit = nmu.units["energy"]
		time_unit = nmu.units["time"]
		
		mEdot = []
		sEdot = []
		for sdata, _ in self.data:
			nmu.set_units(sdata["unit_system"])

			Eorb = nmu.convert_energy(abs(sdata["system"]["Eorb"]), nmu.units["energy"], energy_unit)
			porb = nmu.convert_time(sdata["system"]["Porb"], nmu.units["time"], time_unit)

			mEdot.append( sdata["stats"]["mean{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/porb )		
			sEdot.append( sdata["stats"]["stdv{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/porb )
	
		return mEdot, sEdot

	def get_Porb(self, unit_system="SI"):
		nmu.set_units(system=unit_system)
		energy_unit = nmu.units["energy"]
		time_unit = nmu.units["time"]
		
		Porb = []
		for sdata, _ in self.data:
			nmu.set_units(sdata["unit_system"])
			Porb.append( nmu.convert_time(sdata["system"]["Porb"], nmu.units["time"], time_unit) )

		return Porb

	def get_E(self, unit_system="SI"):
                nmu.set_units(system=unit_system)
                energy_unit = nmu.units["energy"]

                mE = []
		sE = []
                for sdata, _ in self.data:
                        nmu.set_units(sdata["unit_system"])

			Eorb = nmu.convert_energy(abs(sdata["system"]["Eorb"]), nmu.units["energy"], energy_unit)
			mE.append( sdata["stats"]["mean{sum{E}/|Eorb|}"] * Eorb )
			sE.append( sdata["stats"]["stdv{sum{E}/|Eorb|}"] * Eorb )

                return mE, sE

#=================================================
### harmonic average
def harmonic_average(cluster, key, unit_system="SI"):
        """
        computes the harmonic average over the data within cluster
        we compute <Edot> = int dt Edot / int dt
        but we take the integral over Porb, so changing the integration variable to Porb yields

        <Edot> = int dP ( E/P * Edot**-1) * key / int dP ( E/P * Edot**-1 )

        ASSUMES spacing is even in Porb-space (within each cluster!), so the integral's measure drops out
        """
        nmu.set_units(unit_system)
        time_unit = nmu.units["time"]
        energy_unit = nmu.units["energy"]

        num = 0.0
        den = 0.0
        p = 0.0
        for sdata, mdata in cluster.data:
                nmu.set_units(sdata["unit_system"])

		### pull out conversion to Porb measure
                Porb = nmu.convert_time(sdata["system"]["Porb"], nmu.units["time"], time_unit)
                Eorb = nmu.convert_energy(sdata["system"]["Eorb"], nmu.units["energy"], energy_unit)
		Edot = sdata["stats"]["mean{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/Porb

		### pull out statistic we're averaging
		if key == "Edot":
	                val = Edot
		elif key == "E":
			val = sdata["stats"]["mean{sum{E}/|Eorb|}"] * Eorb
		else:
			raise KeyError, "key=%s not understood"%key

		### compute contribution to the average
		integrand = Eorb/(Porb*Edot)
                num += integrand * val
                den += integrand 

                p += integrand * Porb

        return num/den, p/den

#=================================================
# sweeps fitting
#=================================================
def sweeps_Edot_powerlaw(clusters, unit_system="SI"):
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
	### compute harmonic average (via delegation)
	edot, p = np.transpose(np.array( [harmonic_average(cluster, "Edot", unit_system=unit_system) for cluster in clusters] ))

	### fit to function!
	### simple power law!
	logp = np.log(p)
	logedot = np.log(np.abs(edot))

	### linear fit in log-space
	yx = np.sum( logp*logedot )
	y = np.sum( logedot )
	xx = np.sum( logp*logp )
	x = np.sum( logp )
	n = len( logp )

	det = xx*n - x*x
	m = (n*yx - x*y)/det
	b = (-x*yx + xx*y)/det

	return (logp, logedot), (m, b)

def sweeps_E_powerlaw(clusters, unit_system="SI"):
	e, p = np.transpose(np.array( [harmonic_average(cluster, "E", unit_system=unit_system) for cluster in clusters] ))

        ### fit to function!
        ### simple power law!
        logp = np.log(p)
        loge = np.log(np.abs(e))

        ### linear fit in log-space
        yx = np.sum( logp*loge )
        y = np.sum( loge )
        xx = np.sum( logp*logp )
        x = np.sum( logp )
        n = len( logp )

        det = xx*n - x*x
        m = (n*yx - x*y)/det
        b = (-x*yx + xx*y)/det

        return (logp, loge), (m, b)

