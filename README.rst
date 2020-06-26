catwoman: A transit modelling Python package for asymmetric light curves
====================================================
.. image:: docs/cw.png

``catwoman`` is a Python package that models asymmetric transit lightcurves where planets are modelled as two semi-circles with different radii in any orientation, for any radially symmetric stellar limb darkening law. 

``catwoman`` uses the integration algorithm developed for the ``batman`` package (Kreidberg 2015), from which ``catwoman`` builds upon. 

For a detailed introduction and more information please visit https://catwoman.readthedocs.io/.

Installation
=============
You can install ``catwoman`` with pip (recommended):

::

	$ pip install catwoman


Quickstart
===========
Below is an example of how to create and plot a simple lightcurve using ``catwoman``. For a more detailed tutorial please see the documentation linked above.

The first step is to import ``catwoman`` and the packages needed for it to run and to plot the results:

::

        import catwoman
        import numpy as np
        import matplotlib.pyplot as plt

Next, following a similar procedure as to that in ``batman``, initialise a ``TransitParams`` object to store the input parameters of the transit:

::

        params  = catwoman.TransitParams()
        params.t0 = 0.                          #time of inferior conjuction (in days)
        params.per = 1.                         #orbital period (in days)
        params.rp = 0.1                         #top semi-circle radius (in units of stellar radii)
        params.rp2 = 0.1                        #bottom semi-circle radius (in units of stellar radii)
        params.a = 15.                          #semi-major axis (in units of stellar radii)
        params.inc = 90.                        #orbital inclination (in degrees)
        params.ecc = 0.                         #eccentricity
        params.w = 90.                          #longitude of periastron (in degrees)
        params.u = [0.1, 0.3]                   #limb darkening coefficients [u1, u2]
        params.limb_dark = "quadratic"          #limbs darkening model
        params.phi = 0.                         #angle of rotation of top semi-circle (in degrees)

Next make the time array to specify the times we want to calculate the model for:

::

        t = np.linspace(-0.05, 0.05, 1000)

Then, to initialise the model and calculate a light curve:

::

        model = catwoman.TransitModel(params,t)         #initalises model
        flux = model.light_curve(params)                #calculates light curve

To view the light curve:

::

        plt.plot(t, flux)
        plt.xlabel("Time from central transit/days")
        plt.ylabel("Relative flux")
        plt.show()


.. image:: docs/Simplesymmetric.png

To model an asymmetric planet, simply change ``params.rp`` and/or ``params.rp2`` and ``params.phi`` to change the orientation of the system.

