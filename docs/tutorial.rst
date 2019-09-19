.. _tutorial:
  
Tutorial
============

This tutorial explains in detail the different features of ``catwoman``.

Initialising the model
----------------------

As shown in the :ref:`quickstart`, to start setting up the model, one has to initialise a variety of parameters:
::

	import catwoman
	import numpy as np
	import matplotlib as plt
	
	params = catwoman.TransitParams() 	#object to store transit parameters
	params.t0 = 0.				#time of inferior conjuction (in days)
	params.per = 1.				#orbital period (in days)
	params.rp = 0.1				#top semi-circle radius (in units of stellar radii)
	params.rp2 = 0.1005			#bottom semi-circle radius (in units of stellar radii)
	params.a = 15.				#semi-major axis (in units of stellar radii)
	params.inc = 90.			#orbital inclination (in degrees)
	params.ecc = 0.				#eccentricity
	params.w = 90.				#longitude of periastron (in degrees)
	params.u = [0.1, 0.3]                   #limb darkening coefficients [u1, u2]
	params.limb_dark = "quadratic"		#limb darkening model
	params.phi = 0.				#angle of rotation of top semi-circle

	t = np.linspace(-0.04,0.04,1000)	#array of times to calculate the light curves for
	mod = catwoman.TransitModel(params, t, max_err = 0.1)	#initialises the model

As in ``batman``, the initialisation step automatically calculates the array of separation of centres between the star and the planet and also pre-runs the light_curve function numerous times in order to find the approriate integration step size for a given ``max_err``. 
``catwoman`` does this for all the available limb darkenig laws ("quadratic", "logarithmic", "exponential", "nonlinear", "linear", "power2", "uniform" and "custom").
Note: The default for ``max_err`` is 1ppm and describes the allowed error (in ppm) between the smallest integration step size and the selected integration step size. The lower the specified ``max_err``, the smaller the step size and the longer this initialisation step will take to run.

Calculating light curves
-----------------------------  
				  


