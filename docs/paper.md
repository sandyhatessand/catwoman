---
title: 'Catwoman: A transit modelling Python package for asymmetric light curves'
tags:
  - Python
  - astronomy
  - exoplanets
  - transit 
authors:
  - name: Kathryn Jones
    orcid: 
    affiliation: 1
  - name: Néstor Espinoza
    orcid: 
    affiliation: 2
affiliations:
 - name: Max-Planck-Institut f\"ur Astronomie, K\"onigstuhl 17, 69117 Heidelberg, Germany
   index: 1
 - name: Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
   index: 2
date: 1 November 2019
bibliography: paper.bib

---

# Summary

When exoplanets pass in front of their stars from our point of view on Earth, they imprint a transit signature on the stellar light curve which, to date, has been assumed to be symmetric in time, owing to the planet being modelled as a circular area occulting the stellar surface (see, e.g., [@Mandel02]; [@Kreidberg15]; [@Luger19]). However this signature might be asymmetric due to several possible effects, one of which is the different temperature/pressure and/or chemical compositions the different terminator regions a transiting planet could have (see, e.g., [@Powell19]). Being able to model these asymmetric signatures directly from transit light curves could give us an unprecedented glimpse into planetary 3-dimensional structure, helping constrain models of atmospheric evolution, structure and composition.

``catwoman`` is a Python package that models these assymmetric transit light curves, calculating light curves for any radially symmetric stellar limb darkening law and where planets are modelled as two semi-circles, of different radii, using the integration algorithm developed in [@Kreidberg15] and implemented in the ``batman`` library, from which ``catwoman`` builds upon. It is fast and efficient and open source with full documentation available to view at (INSERT readthedocs?).
     
The light curves are modelled as follows: The decrease in flux, $\delta$, as a planet transits its star can be approximated by the sum 

$$\delta = \sum_{i=1}^{N} I\left(x_m\right)\Delta A(x_m,R_{p,1},R_{p,2},\varphi,d)$$,

splitting the semi-circles into iso-intensity bands centred on the star and for each intersectional segment (see Figure 1) you multiply its area, $\Delta A$, by the intensity of the star and then sum these strips to generate the full $\delta$ for a specific separation between the centre of the star and planet, $d$. The code then increments $d$ by a small pre-determined amount (based on the time array given by the user) and recalculates $\delta$.

![Figure 1](strips.png)

The width of the iso-intensity bands determines the truncation error of the model. The model is first initialised with parameters including a maximum truncation error either set by the user or taken as the pre-set value as 1ppm. As in ``batman``, ``catwoman`` first calculates many models, with varying widths and geometrically searches for a width that produces an error less than 1% away (and always less than) the specified level. The model then uses this width value to calculate the desired light curves. A lower specified error, and therefore thinner iso-intensity bands, produces more accurate light curves, however more steps are needed to calculate $\delta$ which takes more time.  

``catwoman`` also allows for $\varphi$, the angle of rotation of the semi-circles, to vary as a free parameter, which is something no other model has tried to implement, accounting for the possibility of spin-orbit misalignments of the planet. The two semi-circle radii, $R_{p,1}$ and $R_{p,2}$, and other orbital variables are also completely free parameters.

``catwoman`` was designed to be used by astronomical researchers. For a realistic light curve with 100 in-transit data points, ``catwoman`` takes around 340 seconds to produce 1 million quadratic-limb-darkened light curves on a single 1.3 GHz Intel Core i5 processor. It is used in Espinoza & Jones (in prep.).

# Acknowledgements


# References
