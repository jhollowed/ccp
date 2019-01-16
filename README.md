# ccp
ccp is a bad acronym that I came up with to mask the even worse title CPAC-convenience-plots. ccp provides a simple set of functions for generating convenient plots to aide in all your favorite cosmological activities.

`ccp.py` provides one class, `ccp.ccp()`, which is constructed as:

```
'''
Class to handle plotting of various convenience figures given 
the relevant parameters of some simulation of interest. Available 
plots are:

snapshot vs. redshift (and scale factor)
snapshot vs. comoving distance (and volume)
snapshot vs. lightcone storage size (and snapshot average)

Feel free to contribute and add more!

Params:
:param snapshots:  A list of snapshot integers. Defaults to None, in 
                   which case the snapshots will be read from AlphaQ
                   outout directories on the Mira filesystem (applies
                   to OuterRim and many other HACC sims)
:param step_range: The step range of interest. Expected a list in the
                   form [minStep, maxStep]. Defaults to [499, 247], 
                   which is a redshift range of [0, 3]. 
:param cosmo:      An Astropy cosmology object. Defaults to WMAP7
:param cosmo_name: Name to associate with the cosmo object as a string. 
                   Defaults to 'WMAP7'
:param z_init:     The initial redshift of the simulation of interest. 
                   Defaults to 200
:param tot_steps:  The total number of steps of the simulation of interest.
                   Defaults to 500
:param lc_dir:     Path to a lightcone top-level directory (expects snapshot-wise 
                   subdirectories). Defaults to the OuterRim full particle 
                   lightcone on the Mira filesystem
:param snap_dir:   Path to a snapshot top-level directory (expects snapshot-wise 
                   subdirectories). Defaults to the OuterRim full particle 
                   snapshots on the Mira filesystem
:param sim_name:   simulation name for figure text
:param rL:         box side length of the simulation of interest, in Mpc/h. 
                   Defaults to 3000 for Outer Rim
:param cm:         the colormap to sample from for all plots; can be any mpl
                   colormap object. Defaults to plt.cm.plasma
:return:           None
'''
```

## Usage

```
import ccp

# plot with default settings (otherwise specify any of the optional params listed above) 
outerRim_figs = ccp.ccp()

outerRim_figs.make_all()
```

