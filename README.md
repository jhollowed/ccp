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

Disclaimer: if you muck around with the too much, the plots may render with poor formatting (if the axes ranges are drastically changed from the default, for instance), but it should be easy enough to look toward the bottom section of each plotting function, under the `# formatting` comment heading, and clean it up.

```python
import ccp
import matplotlib.pyplot as plt

# plot with default settings (otherwise specify any of the optional params listed above) 
outerRim_figs = ccp.ccp()
outerRim_figs.make_all()

# or customize to a range of z=[0,1], with new color scheme
custom_figs = cpp.cpp(step_range=[247, 499], cm=plt.cm.viridis)
custom_figs.step_vs_z()
```

The result of the first call to `make_all()` is shown below, followed by the second customized call

![Outer Rim step vs redshift](sample_figs/or_step_vs_z.png?raw=true "Outer Rim step vs redshift")
![Outer Rim step vs comoving distance/volume](sample_figs/or_step_vs_comv.png?raw=true "Outer Rim step vs comoving distance/volume")
![Outer Rim step vs storage size](sample_figs/or_step_vs_mem.png?raw=true "Outer Rim step vs storage size")
![Customized Outer Rim step vs redshift](sample_figs/or_step_vs_z_alt.png?raw=true "Customized Outer Rim step vs redshift")
