import os
import pdb
import glob
import itertools
import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP7

# =========================================================================

class ccp:

    def __init__(self, snapshots = None, step_range = [121, 499],
                 cosmo = WMAP7, cosmo_name = 'WMAP7', 
                 z_init = 200, tot_steps = 500,
                 lc_dir = '/projects/DarkUniverse_esp/rangel/lc_data/outer_rim', 
                 snap_dir = '/projects/DarkUniverse_esp/heitmann/OuterRim/M000'
                            '/L4225/HACC000/output', 
                 sim_name = 'OuterRim', rL = 3000, cm = plt.cm.plasma):
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

        # get available snapshots from AlphaQ output
        if snapshots is None:
            snapshots = np.array([int(s.split('STEP')[-1]) for s in glob.glob('/projects/'
                                  'DarkUniverse_esp/heitmann/OuterRim/M000/L360/HACC001/'
                                  'analysis/Particles/STEP*')])
        range_mask = np.logical_and(snapshots >= step_range[0], 
                                    snapshots <= step_range[1])
        self.snaps = np.sort(snapshots[range_mask])
        self.a = np.linspace(1/(z_init+1), 1, tot_steps);
        self.z = 1/self.a-1

        self.cosmo = cosmo 
        self.cosmo_name = cosmo_name
        self.sim_name = sim_name
        self.lc_dir = lc_dir
        self.snap_dir = snap_dir
        self.rL = rL

        # mpl setup
        n = 8
        self.cm = cm(np.linspace(0.2, 0.8, n))
        rc('text', usetex=True)


    # =========================================================================


    def step_vs_z(self):
        '''
        Snapshot vs redshift on left y-axs, snapshot vs scale factor on 
        right y-axis. Generates a table listing all available snapshots
        to the right of the axes.
        '''

        # define 1x6 subplot grid layout. Leftmost column will be the
        # snapshot table. Right 5 columns will be the figure 
        tab = plt.subplot2grid((1,6), (0,0))
        ax1 = plt.subplot2grid((1,6), (0,1), colspan=5)
        ax2 = ax1.twinx()
        
        # make snapshot table
        # if len(snaps) is much larger(smaller) than 60, may need to increase(decrease) numCols
        tab.axis('tight')
        tab.axis('off')
        numCols = 3
        snap_strs = self.snaps.astype(str)[::-1]
        if(len(snap_strs)%numCols != 0): snap_strs = np.hstack([snap_strs, ['']])
        snap_strs = np.vstack(np.split(snap_strs, 3)).T
        snap_text = ('$\\mathrm{\\underline{all\\>snapshots}}$\n$' +
                    '$\n$'.join(['\\>\\>\\>'.join(s) for s in snap_strs]) + '$')
        tab.text(0,0.1,snap_text, transform=tab.transAxes, fontsize=12)
      
        # plot snap vs z and snap vs a curves
        l1 = ax1.plot(self.snaps, self.z[self.snaps], lw=2, 
                 label=r'$\mathrm{redshift}$', c=self.cm[0])
        l2 = ax2.plot(self.snaps, self.a[self.snaps], lw=2, 
                 label=r'$\mathrm{scale\>factor}$', c=self.cm[4])
       
        # formatting
        ax1.set_xticks(self.snaps[0::2])
        ax1.set_xticklabels([r'${}$'.format(s) for s in self.snaps[0::2].astype(str)], 
                            rotation= 270,fontsize=13)
        ax1.set_yticks(self.z[self.snaps][0::4])
        ax2.set_yticks(self.a[self.snaps][0::4])
        for tick in ax1.yaxis.get_major_ticks(): tick.label.set_fontsize(13)
        for tick in ax2.yaxis.get_major_ticks(): tick.label2.set_fontsize(13)
        
        orig_xlim = ax1.get_xlim()
        ticks = ax1.get_yticks()
        for i in range(len(ticks)):
            if(i%2==0): continue
            ax1.fill_between([orig_xlim[0], orig_xlim[1]], 
                             np.ones(2)*(ticks[i] - (ticks[i]-ticks[i-1])), 
                             np.ones(2)*ticks[i], 
                             color='k', alpha=0.075)
        ax1.set_xlim(orig_xlim)

        lns = l1+l2
        labs = [l.get_label() for l in lns]
        ax1.grid()
        ax2.grid(linestyle='--')
        ax1.set_xlabel(r'$\mathrm{snapshot}$', fontsize=16)
        ax1.set_ylabel(r'$z$', fontsize=16, color=self.cm[0])
        ax2.set_ylabel(r'$a$', fontsize=16, color=self.cm[4])
        ax2.legend(lns, labs, loc='upper center', fontsize=14)

        fig = plt.gcf()
        fig.set_size_inches(12.5, 6)
        plt.tight_layout()
        plt.show()
    

    # =========================================================================


    def step_vs_comv(self):
        '''
        Snapshot vs. comoving distance on the left y-axis, snapshot vs. 
        comoving volume on the right y-axis.
        '''
    
        # get comoving distance and volume
        d = self.cosmo.comoving_distance(self.z).value
        v = self.cosmo.comoving_volume(self.z).value
        d_box = self.rL / self.cosmo.h
        v_box = d_box ** 3
       
        # plot snap vs d and snap vs v curves
        ax1 = plt.subplot2grid((1,1), (0,0))
        ax2 = ax1.twinx()
        ld = ax1.plot(self.snaps, d[self.snaps], lw=2, 
                      label=r'$\mathrm{{comv.\>distance\>({}\>cosmology)}}$'.format(
                      self.cosmo_name), c=self.cm[2])
        ldb = ax1.plot(self.snaps, np.ones(len(self.snaps))*d_box, '--', c=self.cm[2], 
                       lw=1.5, label=r'$\mathrm{{{}\>box\>length}}$'.format(self.sim_name))
        lv = ax2.plot(self.snaps[0:-1], v[self.snaps][0:-1], lw=2, 
                      label=r'$\mathrm{{comv.\>volume\>({}\>cosmology)}}$'.format(
                      self.cosmo_name), c=self.cm[6])
        lvb = ax2.plot(self.snaps, np.ones(len(self.snaps))*v_box, '--', c=self.cm[6], 
                       lw=1.5, label=r'$\mathrm{{{}\>box\>volume}}$'.format(self.sim_name))
        ax2.set_yscale('log') 
        
        # formatting
        ax1.set_xticks(self.snaps[0::2])
        ax1.set_xticklabels([r'${}$'.format(s) for s in self.snaps[0::2].astype(str)], 
                            rotation= 270,fontsize=13)
        ax1.set_yticks(d[self.snaps][0::4])
        for tick in ax1.yaxis.get_major_ticks(): tick.label.set_fontsize(13)
        for tick in ax2.yaxis.get_major_ticks(): tick.label2.set_fontsize(13)
        
        orig_xlim = ax1.get_xlim()
        ticks = ax1.get_yticks()
        for i in range(len(ticks)):
            if(i%2==0): continue
            ax1.fill_between([orig_xlim[0], orig_xlim[1]], 
                             np.ones(2)*(ticks[i] - (ticks[i]-ticks[i-1])), 
                             np.ones(2)*ticks[i], 
                             color='k', alpha=0.075)
        ax1.set_xlim(orig_xlim)

        lns = lv+lvb+ld+ldb
        labs = [l.get_label() for l in lns]
        ax1.grid()
        ax1.set_xlabel(r'$\mathrm{snapshot}$', fontsize=16)
        ax1.set_ylabel(r'$\mathrm{comv.\>distance\>[Mpc]}$', fontsize=16, color=self.cm[2])
        ax2.set_ylabel(r'$\mathrm{comv.\>volume\>[Mpc}^3]$', fontsize=16, color=self.cm[6])
        ax2.legend(lns, labs, loc='lower left', fontsize=14)

        fig = plt.gcf()
        fig.set_size_inches(10.5, 6)
        fig.tight_layout()
        plt.show()
    
    
    # =========================================================================
    
    
    def step_vs_lcMem(self):
        '''
        Snapshot vs. lightcone storage size, including average snapshot
        storage size w/ 1-std band
        '''
    
        # get subdirectories
        # assume snapshot-wise subdirectories are found in lc_dir and snap_dir,
        # where each subdirectory has the same prefix (could be 'STEP487' or 'lc487'...)
        # other subdirectories present might mess this up if they share the same prefix,
        # or if there's lots of them
        lc_subdirs = [s.split('/')[-2] for s in glob.glob('{}/*/'.format(self.lc_dir))]
        lc_subdirs_split = np.array([["".join(x) for _, x in 
                          itertools.groupby(s, key=str.isdigit)] for s in lc_subdirs])
        mask = np.array([len(l) > 2 for l in lc_subdirs_split])
        lc_prefixes = [l[0] for l in lc_subdirs_split[~mask]]
        lc_uniq_prefixes, counts = np.unique(lc_prefixes, return_counts = True)
        lc_prfx = lc_uniq_prefixes[np.argmax(counts)]

        snap_subdirs = [s.split('/')[-2] for s in glob.glob('{}/*/'.format(self.snap_dir))]
        snap_subdirs_split = np.array([["".join(x) for _, x in 
                          itertools.groupby(s, key=str.isdigit)] for s in snap_subdirs])
        mask = np.array([len(l) > 2 for l in snap_subdirs_split])
        snap_prefixes = [l[0] for l in snap_subdirs_split[~mask]]
        snap_uniq_prefixes, counts = np.unique(snap_prefixes, return_counts = True)
        snap_prfx = snap_uniq_prefixes[np.argmax(counts)]
        
        # get correct subdirs, sort by step, get files, measure sizes.
        # it is assumed that lightcone files have the characters 'lc' in them somewhere...
        # anything else is ignored, including any .SubInput files
        lc_shells = np.array(glob.glob('{}/{}*'.format(self.lc_dir, lc_prfx)))
        lc_snaps = np.array([int(s.split(lc_prfx)[-1]) for s in lc_shells])
        lc_snap_mask = [s in self.snaps for s in lc_snaps]
        lc_shells = (lc_shells[lc_snap_mask])[np.argsort(lc_snaps[lc_snap_mask])]

        all_files = [np.array(glob.glob('{}/*'.format(dir_))) for dir_ in lc_shells]
        skip_files = [np.array([(('lc' in f.split('/')[-1]) and 
                                 ('SubInput' not in f.split('/')[-1])) 
                                 for f in files]) for files in all_files]
        
        lc_sizes = np.array([ sum([os.path.getsize(f) for f in all_files[i][skip_files[i]]])
                            for i in range(len(all_files)) ]) / 1e9 # in GB

        # it is assumed that snapshot files contain 'mpicosmo', and do not contain
        # 'full'-- this is very particular to Outer Rim, because full particle snapshots
        # have been partially deleted. If running for a different simulation, change 
        # skip_files, and the factor of 1/1e7 in snap_mean and snap_std, as needed
        snapshots = np.array(glob.glob('{}/{}*'.format(self.snap_dir, snap_prfx)))
        snap_nums = np.array([int(s.split(snap_prfx)[-1]) for s in snapshots])
        snap_mask = [s in self.snaps for s in snap_nums]
        snapshots = (snapshots[snap_mask])[np.argsort(snap_nums[snap_mask])]
        
        all_files = [np.array(glob.glob('{}/*'.format(dir_))) for dir_ in snapshots]
        skip_files = [np.array([(('mpicosmo' in f.split('/')[-1]) and 
                                 ('full' not in f.split('/')[-1])) 
                                 for f in files]) for files in all_files]
        snap_sizes = np.array([ sum([os.path.getsize(f) for f in all_files[i][skip_files[i]]])
                            for i in range(len(all_files)) ])
        snap_mean = np.mean(snap_sizes) / 1e7 # in GB (extra factor of 100 to account for downsampling)
        snap_std = np.std(snap_sizes) / 1e7

        # plot snap vs lc storage
        ax1 = plt.subplot2grid((1,1), (0,0))
        ax2 = ax1.twinx()
        
        try:
            llc = ax1.plot(self.snaps, lc_sizes, lw=2, 
                          label=r'$\mathrm{{{}\>lightcone}}$'.format(self.sim_name), c=self.cm[2])
        except ValueError:
            # this will happen if the input lightcone does not contain all of the steps
            # specified in self.snaps (lightcone is smaller than desired redshift range)
            llc = ax1.plot(sorted(lc_snaps[lc_snap_mask]), lc_sizes, lw=2, 
                          label=r'$\mathrm{{{}\>lightcone}}$'.format(self.sim_name), c=self.cm[2])
        orig_xlim = ax1.get_xlim()
        
        lsn = ax2.plot(orig_xlim, np.ones(2)*snap_mean, '--', 
                       c=self.cm[5], lw=1, label=r'$\mathrm{{{}\>snapshots}}$'.format(
                       self.sim_name))
        ax2.fill_between(orig_xlim, np.ones(2)*snap_mean - snap_std,
                                   np.ones(2)*snap_mean + snap_std,
                         color=self.cm[5], alpha=.3, lw=0)
        
        # formatting
        ax1.set_xticks(self.snaps[0::2])
        ax1.set_yticks(np.arange(ax1.get_ylim()[0], ax1.get_ylim()[1], 300))
        ax1.set_xticklabels([r'${}$'.format(s) for s in self.snaps[0::2].astype(str)], 
                            rotation= 270,fontsize=13)
        ax2.set_ylim([40000, 42000]) # <-- Outer Rim specific...
        for tick in ax1.yaxis.get_major_ticks(): tick.label.set_fontsize(13)
        
        ticks = ax1.get_yticks()
        for i in range(len(ticks)):
            if(i%2==0): continue
            ax1.fill_between([orig_xlim[0], orig_xlim[1]], 
                             np.ones(2)*(ticks[i] - (ticks[i]-ticks[i-1])), 
                             np.ones(2)*ticks[i], 
                             color='k', alpha=0.075)
        ax1.set_xlim(orig_xlim)

        lns = llc + lsn
        labs = [l.get_label() for l in lns]
        ax1.grid()
        ax1.set_xlabel(r'$\mathrm{snapshot}$', fontsize=16)
        ax1.set_ylabel(r'$\mathrm{Lightcone\>Storage\>Size\>[Gb]}$', fontsize=16, color=self.cm[2])
        ax2.set_ylabel(r'$\mathrm{Snapshot\>Storage\>Size\>[Gb]}$', fontsize=16, color=self.cm[5])
        ax2.legend(lns, labs, loc='upper right', fontsize=14)

        fig = plt.gcf()
        fig.set_size_inches(10.5, 6)
        fig.tight_layout()
        plt.show()   
    
    # =========================================================================


    def make_all(self):
        self.step_vs_z();
        self.step_vs_comv();
        self.step_vs_lcMem();

  
# =========================================================================


if __name__ == '__main__': 
    pl = ccp()
    pl.make_all()
