.. catwoman documentation master file, created by
   sphinx-quickstart on Tue Sep 17 13:30:53 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: cw.png

Welcome to catwoman's documentation!
====================================

When exoplanets pass in front of their stars, they imprint a transit signature on the stellar light curve which to date has been assumed to be symmetric in time, owing to the planet being modelled as a circular area occulting the stellar surface. However, this signature might be asymmetric due to different temperature/pressure and/or chemical compositions in the different terminator regions of the transiting planet.

``catwoman`` is a Python package that allows to model these asymmetric transit lightcurves, calculating light curves for any radially symmetric stellar limb darkening law, and where planets are modelled as two semi-circles, of different radii, using the integration algorithm developed in Kreidberg (2015) and implemented in the ``batman`` library, from which ``catwoman`` builds upon.

Please cite the paper Espinoza & Jones (in prep.) if you use catwoman in your research.   

If you find a bug or have any problems with catwoman, please `opening an issue <https://github.com/KathrynJones1/catwoman/issues>`_ on the project's GitHub and we will try to get back to you as soon as possible.  

Table of Contents
==================

.. toctree::
   :maxdepth: 2
	
   installation
   quickstart
   tutorial
   API
	

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
