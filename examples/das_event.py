from obspy import read
from kvp import KVP


# Read traces
traces = read('OT/*', format='sac')

# Create a KVP instance with some parameters
kvp_picker = KVP(
    freqmax=20.0,
    octaves=4,
    voices=2,
    cf_cycles=90.0,
    jmp_cycles=2.0,
    jump=5.0,
    mingap=0.25,
    nbands=1,
    )


# Run KVP on all traces
picks = [kvp_picker.obspy(tr) for tr in traces]


# Print all outputs
print('Station name\tUTC time')
for channel in picks:
    print(channel)


# Zip traces and picks together
items = [(trace, pick) for trace, pick in zip(traces, picks)]



# Plot function
def plot_OT(items, fnum=0, show=True, save=False):

    
    import numpy as np
    from obspy import UTCDateTime
    import matplotlib.pyplot as plt
    from matplotlib.offsetbox import AnchoredText
    
    
    plt.rcParams["font.family"] = "monospace"
    
    labelsize = 8
    legendsize = 8
    textsize = 8
    ticksize = 8
    
    t_origin = UTCDateTime('2021-10-26T16:25:37')
    
    n_tr = len(items)
    
    fig, ax = plt.subplots(n_tr, 1, dpi=300, figsize=(5,0.8*n_tr), sharex=True)
    fig.subplots_adjust(hspace=0, wspace=0)
    
    for item, axis in zip(items, ax):
        
        tr_obspy, picks = item
        
        offset_t = tr_obspy.stats.starttime - t_origin
        t_axis = np.linspace(0,tr_obspy.stats.npts/tr_obspy.stats.sampling_rate, tr_obspy.stats.npts) + offset_t
        
        axis.tick_params(labelleft=False, labelsize=ticksize)
        axis.tick_params(axis='y', direction='in')
        
        ylim = 30
        axis.set_xlim([18.5,21.1])
        axis.set_ylim([-ylim,+ylim])
        
        axis.set_xticks(np.arange(18.6,21.1,0.3))
        axis.set_yticks([-20,0,+20])
        
        
        axis.plot(t_axis, tr_obspy, linewidth=0.25, color='gray', zorder=0)
        
        for pick in picks:
            ons = UTCDateTime(pick.onset())- t_origin
            color = 'k' if pick.nb>= 3 else '0.5'
            axis.axvline(ons, linewidth=0.75, color=color, zorder=4)
        
        sta = tr_obspy.stats.station
        
        prop = dict(size=textsize, alpha=0.8)
        fqtext = AnchoredText(f'ZI.{sta}.HS1', loc='center', bbox_to_anchor=(0.135,0.79),
                          bbox_transform=axis.transAxes, prop=prop, frameon=True)
        fqtext.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        fqtext.patch.set_alpha(0.5)
        axis.add_artist(fqtext)
            
        
        
    ax[0].tick_params(top=False, labelleft=True)
        
    ax[-1].set_xlabel('Time [s]', fontsize=labelsize)
    
    axis = ax[0]
    
    axis.set_ylabel('n$\\epsilon$ [nm/m]', fontsize=labelsize)
    
    # Ghost labels
    axis.axvline(-1, linewidth=0.75, color='k', label='KVP')
    axis.axvline(-1, linewidth=0.75, color='0.5', label='KVP < 3 bands')
    
    vpos = 0.93
    fig.text(0.5, vpos, 'OT :: 2021-10-26T16:25:37', fontsize=textsize, horizontalalignment='center')
    
    legend_vpos = 0.901
    handles, labels = axis.get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(0.511, legend_vpos), loc='center',
               ncol=4, frameon=False, fontsize=legendsize, labelspacing=0.1)
    
    if save:
        plt.savefig(f'KVP_OT_example_{fnum}.pdf', bbox_inches='tight')
    
    if show:
        plt.show()
    
    plt.close(fig)


# Plot in groups of 6 (manuscript figure)
plot_OT(items[0:6], fnum=1)
plot_OT(items[6:12], fnum=2)
plot_OT(items[12:18], fnum=3)
plot_OT(items[18:24], fnum=4)
plot_OT(items[24:30], fnum=5)
plot_OT(items[30:36], fnum=6)