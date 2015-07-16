usage="""utilities for fitting derived data products"""

import numpy as np
import pickle
import nmode_utils as nmu

import pygsl.multifit_nlin as pygsl_mfN # fitting functions
import pygsl.errno as pygsl_errno

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

	def get_Mcomp(self, unit_system="SI"):
		nmu.set_units(system=unit_system)
		mass_unit = nmu.units["mass"]

		Mcomp = []
		for sdata, _ in self.data:
			Mcomp.append( sdata["system"]["Mcomp/Mjup"]*nmu.Mjup )

		return Mcomp

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

	def get_Eorb(self, unit_system="SI"):
		nmu.set_units(system=unit_system)
		energy_unit = nmu.units["energy"]

		Eorb = []
		for sdata, _ in self.data:
			nmu.set_units(sdata["unit_system"])
			Eorb.append( nmu.convert_energy(sdata["system"]["Eorb"], nmu.units["energy"], energy_unit) )

		return Eorb

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

	def get_Mprim(self, unit_system="SI"):
		nmu.set_units(system=unit_system)
		mass_unit = nmu.units["mass"]

		Mprime = []
		for sdata, _ in self.data:
			nmu.set_units(sdata["unit_system"])
			Mprime.append( nmu.convert_mass(sdata["system"]["Mprim/Msun"]*nmu.Msun, nmu.units["mass"], mass_unit) )
		return Mprime

	def get_Rprim(self, unit_system="SI"):
		nmu.set_units(system=unit_system)
		length_unit = nmu.units["length"]

		Rprime = []
		for sdata, _ in self.data:
			nmu.set_units(sdata["unit_system"])
			Rprime.append( nmu.convert_length(sdata["system"]["Rprim/Rsun"]*nmu.Rsun, nmu.units["length"], length_unit) )
		return Rprime

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

	P = []
	E = []
	mEdot = []
	sEdot = []
	x = []
	sx = []

	for sdata, mdata in cluster.data:
		nmu.set_units(sdata["unit_system"])

		### convertion to Porb measure
		Porb = nmu.convert_time(sdata["system"]["Porb"], nmu.units["time"], time_unit)
		Eorb = nmu.convert_energy(sdata["system"]["Eorb"], nmu.units["energy"], energy_unit)

		P.append( Porb )
		E.append( Eorb )

		mEdot.append( sdata["stats"]["mean{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/Porb )
		sEdot.append( sdata["stats"]["stdv{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/Porb )

		if key == "Edot":
			x.append( sdata["stats"]["mean{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/Porb )
			sx.append( sdata["stats"]["stdv{|sum{Edot}|*(Porb/|Eorb|)}"] * Eorb/Porb )
		elif key == "E":
			x.append( sdata["stats"]["mean{sum{E}/|Eorb|}"] * Eorb )
			sx.append( sdata["stats"]["stdv{sum{E}/|Eorb|}"] * Eorb )
		else:
			raise KeyError, "key=%s not understood"%key

	### cast to arrays
	P = np.array(P)
	E = np.array(E)
	mEdot = np.array(mEdot)
	sEdot = np.array(sEdot)
	x = np.array(x)
	sx = np.array(sx)

	### compute averages
	integrand = Eorb/(Porb*mEdot)

	num = np.sum( integrand * x )
	den = np.sum( integrand )

	p = np.sum( integrand * P )

	mx = num/den
	mp = p/den

	### compute uncertainties
	T = np.sum( integrand )
	integrate = E/(P*mEdot**2)
	a = np.sum( integrand**2 * sx**2 ) / T**2
	b = np.sum( (integrate * ( x*T - np.sum( x*integrand) ))**2 * sEdot**2 ) / T**4

	smx = a + b ### uncertainty in x
	smp = np.sum( (integrate * ( P*T - np.sum( P*integrand) ))**2 * sEdot**2 ) / T**4 ### uncertainty in p


	return (mx, smx**0.5) , (mp, smp**0.5)


	'''
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
	'''
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
	(edot, p), (sedot, sp) = np.transpose(np.array( [harmonic_average(cluster, "Edot", unit_system=unit_system) for cluster in clusters] ))
#	edot, p = np.transpose(np.array( [harmonic_average(cluster, "Edot", unit_system=unit_system) for cluster in clusters] ))

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

	return (logp, logedot), (sedot, sp), (m, b)

###
def sweeps_E_powerlaw(clusters, unit_system="SI"):
	""" fits stellar energy to a power law of orbital period """
	(e, p), (se, sp) = np.transpose(np.array( [harmonic_average(cluster, "E", unit_system=unit_system) for cluster in clusters] ))
#	e, p = np.transpose(np.array( [harmonic_average(cluster, "E", unit_system=unit_system) for cluster in clusters] ))

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

        return (logp, loge), (se, sp), (m, b)

#========================
# helper methods for more general fitting
#========================
###
def chi2(y, yhat, reduced=False):
	"""
	computes the chi2 statistic for y, yhat

	if reduced, returns the reduced chi2 statistic
	"""

	chi2_stat = np.sum( (1.0 - y/yhat)**2 )

	if reduced:
		chi2_stat /= len(y)

	return chi2_stat

###
def steps(x, params):
	"""
	computes y(x) = sum C * x**a * (1 + e**(b*x))**(-1)

	and params = [ C, a, b, C1, a1, b1, ...]
	"""
	if not isinstance(x, np.ndarray):
		x = np.array( x )

	y = np.zeros_like(x, float)
	for n in xrange(len(params)/3):
		C, a, b = params[3*n:3*(n+1)]
		y += C * x**a / (1+np.exp(b*x))

	return y

def dsteps(x, s, params):
	"""
	computes the partial derivative matrix for steps
	"""
	if not isinstance(x, np.ndarray):
		x = np.array( x )

	dfit = []
	for n in xrange(len(params)/3):
		C, a, b = params[3*n:3*(n+1)]
		dfit.append( x**a / (1+np.exp(b*x)) ) ### dfit/dC
		dfit.append( C * x**a * np.log(x) / (1+np.exp(b*x)) ) ### dfit/da
		dfit.append( -C * x**(a+1) * np.exp(b*x) / (1+np.exp(b*x))**2 ) ### dfit/db

	return np.transpose( np.array(dfit) / s )

### 
def steps_fit(params, data):
	x, y, s = data
	return (y - steps(x, params)) / s

###
def steps_dfit(params, data):
	x, y, s = data
	return dsteps(x, s, params)

###
def steps_fitdfit(params, data):
	return steps_fit(params, data), steps_dfit(params, data)

###
def steps_fitter(x, y, sigma, n, max_iters=50, rtol=1e-8, verbose=False):
	"""
	fits data to "steps" with n terms
	"""
	if not isinstance(x, np.ndarray):
		x = np.array( x )
	if not isinstance(y, np.ndarray):
		y = np.array( y )
	if not isinstance(sigma, np.ndarray):
		sigma = np.array( sigma )


	n_pts = len(x) ### number of points
	p = 3*n ### number of fitting parameters

	### set up system
	if verbose:
		print "setting up objects"
	mysys = pygsl_mfN.gsl_multifit_function_fdf( steps_fit, steps_dfit, steps_fitdfit, np.array( [x, y, sigma] ), n_pts, p)
	solver = pygsl_mfN.lmsder(mysys, n_pts, p)

	### starting points
	if verbose: 
		print "setting up stating point"
	params = []
	for i in xrange(n):
		a = -8.2
		b = 2.6e-6
		C = y[0]/x[0]**a
#		C = np.mean(y)/np.mean(x)**a
		params += [ C, a, b]

	solver.set(params) ### set starting point

	if verbose:
		s = "# %5s"
		S = "  %5d"%0
		for i in xrange(n):
			s += " %9s %9s %9s"%("C", "a", "b")
			S += " %.7f %.7f %7f"%tuple(params[3*i:3*(i+1)])
		print s
		print S

	for iter in range(1,max_iters+1):
		status = solver.iterate() # move fit params
		_fit_params = solver.getx() # new guess
		dfit_params = solver.getdx() # change in fit params
		fits = np.array( solver.getf() ) # residuals at every data point

		J = solver.getJ() # jacobian of fit
		tdx = pygsl_mfN.gradient( J, fits ) # gradient at fit
		status = pygsl_mfN.test_delta(dfit_params, _fit_params, rtol, rtol) # just copied, not understood...

		fn = np.sum((fits)**2)**0.5 # sum square errors

		if status == pygsl_errno.GSL_SUCCESS:
			if verbose: 
				print "# Converged :"
			break
		if verbose: 
			S = "  %5d"%iter
			for i in xrange(n):
				S += " %.7f %.7f %7f"%tuple(_fit_params[3*i:3*(i+1)])
			S += " %.7f"%fn
			print S
	else:
		if verbose:
			print "WARNING! Number of Iterations exceeded in nmode_state.broken_PowLaw_fitter!"
			print "continuing with best guess after %d iterations" % max_iters

	# get error bars on fit params
	covar = pygsl_mfN.covar(solver.getJ(), 0.0) # covariance matrix
	red_chi2 = 1.0*sum( (steps(x, _fit_params) - y)**2 / y )/len(x)


	if verbose:
                s = "#"
                S = " "
		ss= " "
                for i in xrange(n):
                        s += " %9s %9s %9s"%("C", "a", "b")
                        S += " %.7f %.7f %7f"%tuple(_fit_params[3*i:3*(i+1)])
			ss+= " %.7f %.7f %7f"%tuple(covar[i][i] for i in xrange(3*i, 3*(i+1)))
		print s
		print S

	return _fit_params, [covar[i][i] for i in range(len(covar))], red_chi2
