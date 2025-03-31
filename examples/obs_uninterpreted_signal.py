from obspy import read
from kvp import KVP



# Read trace
trace = read('SAFE/*', format='sac')[0]

# Create a KVP instance with some parameters
kvp_picker = KVP(
    freqmax=20.0,
    octaves=3,
    voices=3,
    cf_cycles=90.0,
    jmp_cycles=2.0,
    jump=5.0,
    mingap=0.25,
    nbands=3,
    )



def plot_OBS02(trace, show=True, save=False):
    
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    plt.rcParams["font.family"] = "monospace"
    
    # Get picks with traces included
    picks = kvp_picker.obspy(trace, wt_traces=True, cf_traces=True)
    
    wt = picks.wt_traces()
    cf = picks.cf_traces()
    
    
    duration = trace.stats.endtime - trace.stats.starttime
    
    fig, ax = plt.subplots(len(cf)+1, 1, figsize=(5,7), dpi=300)
    fig.subplots_adjust(hspace=0, wspace=0)
    
    twins = []
    
    for axis, wt, cf in zip(ax[1:], wt, cf):
        
        tr = cf[0]
        offset_t = cf[1]['starttime'] - trace.stats.starttime.timestamp
        t_axis = np.linspace(0, duration, tr.size) + offset_t
        
        ft = wt[0]
        t_axis_ft = np.linspace(0, duration, ft.size)
        
        axis2 = axis.twinx()
        twins.append(axis2)
        
        axis.plot(t_axis, tr, linewidth=1, color='k', zorder=2)
        axis2.plot(t_axis_ft, ft, linewidth=0.25, color='0.8', zorder=0)
        axis.axhline(3+5, linewidth=0.25, color='0.3', linestyle='dashdot', zorder=1)
        axis.axhline(3+3, linewidth=0.25, color='0.3', linestyle='dashed', zorder=1)
        axis.axhline(3, linewidth=0.25, color='0.3', linestyle='dotted', zorder=1)
        
        axis.set_zorder(axis2.get_zorder()+1)
        axis.patch.set_visible(False)
        
        axis.set_xlim([20,80])
        axis.set_ylim([0,11])
        axis.set_yticks([3,6,9])
        
        axis.tick_params(labelbottom=False, bottom=False, labelleft=False)
        axis2.tick_params(right=False, labelright=False)
        
        ylim_ = np.max(np.abs(ft))*1.1
        
        axis2.set_ylim([-ylim_,+ylim_])
        
        fq = round(wt[2],1)
        
        box_text = f'{fq} Hz'
        box_props = dict(boxstyle='round', facecolor='white', alpha=0)
        axis.text(0.02, 0.95, box_text, fontsize=8, transform=axis.transAxes, verticalalignment='top', bbox=box_props,zorder=3)
    
    ax[-1].set_xlabel('Time [s]', fontsize=8)
    ax[1].tick_params(labelleft=True, labelsize=8)
    ax[1].set_ylabel('Kurtosis', fontsize=6)
    
    vpos = 0.915
    fig.text(0.5, vpos, 'SAFE :: 7M.OBS02.EH1 signal (uninterpreted)', fontsize=8, horizontalalignment='center')
    
    # Ghost labels
    ax[0].axvline(-1, linewidth=0.25, color='0.3', linestyle='dotted', zorder=0, label='Gaussian')
    ax[0].axvline(-1, linewidth=0.25, color='0.3', linestyle='dashed', zorder=0, label='Gaussian + 3')
    ax[0].axvline(-1, linewidth=0.25, color='0.3', linestyle='dashdot', zorder=0, label='Gaussian + 5')
    
    legend_vpos = 0.895
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, bbox_to_anchor=(0.511, legend_vpos), loc='center',
               ncol=4, frameon=False, fontsize=8, labelspacing=0.1)
    
    s_t_axis = np.linspace(0,duration,trace.stats.npts)
    ax[0].plot(s_t_axis, trace, linewidth=0.2, color='0.5')
    
    ylim = np.max(np.abs(trace))*1.1
    
    ax[0].tick_params(bottom=False, labelbottom=False, labelsize=8)
    
    ax[0].set_xlim([20,80])
    ax[0].set_ylim([-ylim,+ylim])
    ax[0].set_yticks([-20000,0,+20000])
    ax[0].set_yticklabels([-2,0,+2])
    
    ax[0].set_ylabel('X $10^{4}$\n[counts]', fontsize=6)
    
    twins[0].set_ylabel('Normalized\namplitude', fontsize=6)
    
    ax[-1].tick_params(bottom=True, labelbottom=True, labelsize=6)
    
    ax[-1].set_xticks([20,25,30,35,40,45,50,55,60,65,70,75,80])
    ax[-1].set_xticklabels([0,5,10,15,20,25,30,35,40,45,50,55,60])
    
    box_text = 'Raw HP 1 Hz'
    box_props = dict(boxstyle='round', facecolor='white', alpha=0)
    ax[0].text(0.0105, 0.93, box_text, fontsize=8, transform=ax[0].transAxes, verticalalignment='top', bbox=box_props,zorder=3)
    
    
    # Pick onsets manually
    manual_ons = [33.3, 39.6, 50.0, 52.6]
    
    for ons in manual_ons:
        for axis in ax:
            axis.axvline(ons, linewidth=0.25, color='k')
    
    
    # Automatic picks, change picking parameters (e.g. jump=2.0, nbands=1)
    for pick in picks:
        ons = pick.onset(posix=False)
        ax[0].axvline(ons, linewidth=0.25, color='r')
    
    
    if save:
        plt.savefig('KVP_SAFE_example.pdf', bbox_inches='tight')
    
    if show:
        plt.show()
    plt.close(fig)



plot_OBS02(trace)
