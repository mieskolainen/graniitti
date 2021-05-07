# Advanced histogramming & automated plotting functions
# 
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import pathlib
import matplotlib
matplotlib.use('Agg') # Important for multithreaded applications
from matplotlib import pyplot as plt

import numpy as np
import math
import copy


def chi2_cost(h_mc, h_data):
    """
    Chi2 cost function between two histograms
    """
    counts_mc   = h_mc.counts * h_mc.binscale
    err_mc      = h_mc.errs   * h_mc.binscale

    counts_data = h_data.counts * h_data.binscale
    err_data    = h_data.errs   * h_data.binscale

    return np.sum((counts_mc - counts_data)**2 / (err_mc**2 + err_data**2))


def set_global_style(dpi=120, figsize=(4,3.75), font='serif', font_size=8, legend_fontsize=7, legend_handlelength=1):
    """ Set global plot style.
    """
    plt.rcParams['legend.fontsize']     = legend_fontsize
    plt.rcParams['legend.handlelength'] = legend_handlelength

    plt.rcParams['figure.dpi']     = dpi
    plt.rcParams['figure.figsize'] = figsize

    plt.rcParams['font.family']    = font
    plt.rcParams['font.size']      = font_size


# Colors
imperial_dark_blue  = (0, 0.24, 0.45)
imperial_light_blue = (0, 0.43, 0.69)
imperial_dark_red   = (0.75, 0.10, 0.0)
imperial_green      = (0.0, 0.54, 0.23)


def colors(i, power=0.34):

    c = [imperial_dark_red, imperial_dark_blue, imperial_green, imperial_light_blue]

    if i < len(c):
        return c[i]
    else:
        return c[i%len(c)] * (1.0/power)


""" Global marker styles

zorder : approximate plotting order
lw     : linewidth
ls     : linestyle
"""
errorbar_style  = {'zorder': 3, 'ls': ' ', 'lw': 1, 'marker': 'o', 'markersize': 2.5}
plot_style      = {'zorder': 2, 'ls': '-', 'lw': 1}
hist_style_step = {'zorder': 0, 'ls': '-', 'lw': 1, 'histtype': 'step'}
hist_style_fill = {'zorder': 0, 'ls': '-', 'lw': 1, 'histtype': 'stepfilled'}
hist_style_bar  = {'zorder': 0, 'ls': '-', 'lw': 1, 'histtype': 'bar'}


class hobj:
    """ Minimal histogram data object.
    """
    def __init__(self, counts = 0, errs = 0, bins = 0, cbins = 0, binscale=1.0):
        self.counts   = counts
        self.errs     = errs
        self.bins     = bins
        self.cbins    = cbins
        self.binscale = binscale

        if (np.sum(counts) == 0):
            self.is_empty = True
        else:
            self.is_empty = False

    # + operator
    def __add__(self, other):

        if (self.is_empty == True): # Take the rhs
            return other

        if ((self.bins == other.bins).all() == False):
            raise(__name__ + ' + operator: cannot operate on different sized histograms')

        # Harmonic sum
        binscale = 1/(1/self.binscale + 1/other.binscale)

        counts = self.counts + other.counts
        errs   = np.sqrt(self.errs**2   + other.errs**2)

        return hobj(counts, errs, bins, cbins, binscale)
    
    # += operator
    def __iadd__(self, other):

        if (self.is_empty == True): # Still empty
            return other

        if ((self.bins == other.bins).all() == False):
            raise(__name__ + ' += operator: cannot operate on different sized histograms')

        self.counts = self.counts + other.counts
        self.errs   = np.sqrt(self.errs**2 + other.errs**2)
        
        # Harmonic sum
        self.binscale = 1/(1/self.binscale + 1/other.binscale)

        return self


def stepspace(start, stop, step):
    """ Linear binning edges between [start, stop]
    """
    return np.arange(start, stop + step, step)


def plot_horizontal_line(ax, color=(0.5,0.5,0.5), linewidth=0.9):
    """ For the ratio plot
    """
    xlim = ax.get_xlim()
    ax.plot(np.linspace(xlim[0], xlim[1], 2), np.array([1,1]), color=color, linewidth=linewidth)


def tick_calc(lim, step, N=6):
    """ Tick spacing calculator.
    """
    return [np.round(lim[0] + i*step, N) for i in range(1+math.floor((lim[1]-lim[0])/step))]

def set_axis_ticks(ax, ticks, dim='x'):
    """ Set ticks of the axis.
    """
    if   (dim == 'x'):
        ax.set_xticks(ticks)
        ax.set_xticklabels(list(map(str, ticks)))
    elif (dim == 'y'):
        ax.set_yticks(ticks)
        ax.set_yticklabels(list(map(str, ticks)))

def tick_creator(ax, xtick_step=None, ytick_step=None, ylim_ratio=(0.7, 1.3),
        ratio_plot=True, minorticks_on=True, ytick_ratio_step=0.15, labelsize=9,
        labelsize_ratio=8, **kwargs) :
    """ Axis tick constructor.
    """

    # Get limits
    xlim = ax[0].get_xlim()
    ylim = ax[0].get_ylim()

    # X-axis
    if (xtick_step is not None):
        ticks = tick_calc(lim=xlim, step=xtick_step)
        set_axis_ticks(ax[-1], ticks, 'x')

    # Y-axis
    if (ytick_step is not None):    
        ticks = tick_calc(lim=ylim, step=ytick_step)
        set_axis_ticks(ax[0], ticks, 'y')

    # Y-ratio-axis
    if ratio_plot:
        ax[0].tick_params(labelbottom=False)
        ax[1].tick_params(axis='y', labelsize=labelsize_ratio)

        ticks = tick_calc(lim=ylim_ratio, step=ytick_ratio_step)
        ticks = ticks[1:-1] # Remove the first and the last
        set_axis_ticks(ax[1], ticks, 'y')

    # Tick settings
    for a in ax:
        if minorticks_on: a.minorticks_on()
        a.tick_params(top=True, bottom=True, right=True, left=True, which='both', direction='in', labelsize=labelsize)

    return ax

def create_axes(xlabel='$x$', ylabel=r'Counts', ylabel_ratio='Ratio',
    xlim=(0,1), ylim=None, ylim_ratio=(0.7, 1.3),
    ratio_plot=True, figsize=(5,4), fontsize=9, units={'x': '', 'y': ''}, **kwargs):
    """ Axes creator.
    """
    
    # Create subplots
    N = 2 if ratio_plot else 1
    gridspec_kw = {'height_ratios': (3.333, 1) if ratio_plot else (1,), 'hspace': 0.0}
    fig, ax = plt.subplots(N,  figsize=figsize, gridspec_kw=gridspec_kw)
    ax = [ax] if (N == 1) else ax

    # Axes limits
    for a in ax: a.set_xlim(*xlim)

    if ylim is not None:
        ax[0].set_ylim(*ylim)

    # Axes labels
    if kwargs['density']:
        ylabel = f'$1/N$  {ylabel} / [{units["x"]}]'
    else:
        ylabel = f'{ylabel}  [{units["y"]} / {units["x"]}]'
    xlabel = f'{xlabel} [{units["x"]}]'
    
    ax[0].set_ylabel(ylabel, fontsize=fontsize)
    ax[-1].set_xlabel(xlabel, fontsize=fontsize)

    # Ratio plot
    if ratio_plot:
        ax[1].set_ylabel(ylabel_ratio, fontsize=fontsize)
        ax[1].set_ylim(*ylim_ratio)

    # Setup ticks
    ax = tick_creator(ax=ax, ratio_plot=ratio_plot, **kwargs)

    return fig, ax


def ordered_legend(ax=None, order=None, frameon=False, unique=False, **kwargs):
    """ Ordered legends.
    """

    def unique_everseen(seq, key=None):
        seen = set()
        seen_add = seen.add
        return [x for x,k in zip(seq,key) if not (k in seen or seen_add(k))]

    if ax is None: ax=plt.gca()
    handles, labels = ax.get_legend_handles_labels()

    # Sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    
    # Sort according to a given list, which may be incomplete
    if order is not None: 
        keys=dict(zip(order,range(len(order))))
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t, keys=keys: keys.get(t[0],np.inf)))

    # Keep only the first of each handle
    if unique:  labels, handles= zip(*unique_everseen(zip(labels,handles), key = labels)) 
    ax.legend(handles, labels, frameon=frameon, **kwargs)

    return (handles, labels)


def binwidth(bins):
    """ Return binwidth from a linear array """
    return (bins[1:] - bins[0:-1])


def edge2centerbins(bins) :
    """ Get centerbins from edgebins.
    """
    return (bins[1:] + bins[0:-1])/2


def ratioerr(A, B, sigma_A, sigma_B, sigma_AB = 0, EPS = 1E-15):
    """ Ratio f(A,B) = A/B error, by Taylor expansion of f.
    """
    A[np.abs(A) < EPS] = EPS
    B[np.abs(B) < EPS] = EPS
    return np.abs(A/B) * np.sqrt((sigma_A/A)**2 + (sigma_B/B)**2 - 2*sigma_AB/(A*B))


def hist_to_density(counts, errs, bins):
    """ Normalize to unit integral density function over the visible histogram range """

    norm    = binwidth(bins) * counts.sum()
    return counts/norm, errs/norm


def hist_to_density_fullspace(counts, errs, bins, totalweight):
    """ Normalize histogram to a unit integral density function
        over total sum of event weights (not just the visible histogram range mass)
    """
    norm    = binwidth(bins) * totalweight
    return counts/norm, errs/norm


def hist(x, bins=30, density=False, weights=None):
    """ Calculate a histogram.
    """
    x = np.array(x)

    # Calculate histogram
    if weights is None:
        weights = np.ones(x.shape)

    weights = np.array(weights)

    counts, bins = np.histogram(x, bins=bins, weights=weights)
    cbins = edge2centerbins(bins)

    # Input data to histogram bins
    # Weighted error on bin counts given by (square root of) sum of squared weights
    inds = np.digitize(x, bins)
    errs = np.asarray([np.linalg.norm(weights[inds==k],2) for k in range(1, len(bins))])

    # Density integral 1 over the histogram bins range
    if density:
        counts, errs = hist_to_density(counts=counts, errs=errs, bins=bins)

    return counts, errs, bins, cbins


def hist_obj(x, bins=30, weights=None):
    """ A wrapper to return a histogram object.
    """
    counts, errs, bins, cbins = hist(x, bins=bins, weights=weights)
    return hobj(counts, errs, bins, cbins)


def generate_colormap():
    """ Default colormap.
    """
    # Take colors
    color = plt.cm.Set1(np.linspace(0,1,10))

    # Add black 
    black = np.ones((1,4))
    black[:,0:3] = 0.0
    color = np.concatenate((black, color))    

    return color


def hist_filled_error(ax, bins, cbins, y, err, color, **kwargs):
    """ Stephist style error.
    """
    new_args = kwargs.copy()
    new_args['lw'] = 0
    new_args.pop('histtype', None) # Remove

    ax.fill_between(bins[0:-1], y-err, y+err, step='post', alpha=0.3, color=color, **new_args)

    # The last bin
    ax.fill_between(bins[-2:],  (y-err)[-2:], (y+err)[-2:], step='pre', alpha=0.3, color=color, **new_args)


def superplot(data, observable=None, ratio_plot=True, yscale='linear', ratio_error_plot=True, \
    legend_counts=False, color=None, legend_properties={'fontsize': 9}, bottom_PRC=5, EPS=1E-12):
    """ Superposition (overlaid) plotting
    """
    if observable == None:
        observable = data[0]['obs']
    
    fig, ax = create_axes(**observable, ratio_plot=ratio_plot)

    if color == None:
        color = generate_colormap()
    
    legend_labels  = []

    # y-axis limit
    bottom_count  = 1e32
    ceiling_count = 0

    # Plot histograms
    for i in range(len(data)):

        if data[i]['hdata'].is_empty:
            print(__name__ + f'.superplot: Skipping empty histogram for entry {i}')
            continue

        c = data[i]['color']
        if c is None: c = color[i]

        counts = data[i]['hdata'].counts * data[i]['hdata'].binscale
        errs   = data[i]['hdata'].errs   * data[i]['hdata'].binscale
        bins   = data[i]['hdata'].bins
        cbins  = data[i]['hdata'].cbins

        # -----------------------------------------------
        # ** For visualization autolimits **
        # Use percentile for the bottom (~ handle noisy small bins)
        bottom_count  = np.min([bottom_count, np.percentile(counts[counts > EPS], bottom_PRC)])
        ceiling_count = np.max([ceiling_count, np.max(counts[counts > 0])])
        # -----------------------------------------------
        
        label = data[i]['label']
        if legend_counts == True:
            label += f' $N={np.sum(data[i]["hdata"].counts):.1f}$'

        legend_labels.append(label)

        if   data[i]['hfunc'] == 'hist' :
            ax[0].hist(x=cbins, bins=bins, weights=counts, color=c, label=label, **data[i]['style'])
            hist_filled_error(ax=ax[0], bins=bins, cbins=cbins, y=counts, err=errs, color=c, **data[i]['style'])

        elif data[i]['hfunc'] == 'errorbar' :
            ax[0].errorbar(x=cbins, y=counts, yerr=errs, color=c, label=label, **data[i]['style'])

        elif data[i]['hfunc'] == 'plot' :
            ax[0].plot(cbins, counts, color=c, label=label, **data[i]['style'])
            
            new_args = data[i]['style'].copy()
            new_args['lw'] = 0
            ax[0].fill_between(cbins, counts-errs, counts+errs, alpha=0.3, color=c, **new_args)

    # Plot ratiohistograms
    if ratio_plot:

        plot_horizontal_line(ax[1])
        
        for i in range(len(data)):

            if data[i]['hdata'].is_empty:
                print(__name__ + f'.superplot: Skipping empty histogram for entry {i} (ratioplot)')
                continue

            c = data[i]['color']
            if c is None: c = color[i]

            A        = data[i]['hdata'].counts * data[i]['hdata'].binscale
            B        = data[0]['hdata'].counts * data[0]['hdata'].binscale
            sigma_A  = data[i]['hdata'].errs   * data[i]['hdata'].binscale
            sigma_B  = data[0]['hdata'].errs   * data[0]['hdata'].binscale
            sigma_AB = 0
            ratio_errs = ratioerr(A=A, B=B, sigma_A=sigma_A, sigma_B=sigma_B, sigma_AB=sigma_AB)

            EPS      = 1E-30
            ratio    = A / (B + EPS)
            bins     = data[i]['hdata'].bins
            cbins    = data[i]['hdata'].cbins

            # If no errors turned on
            if ratio_error_plot == False:
                ratio_errs = np.zeros(ratio_errs.shape)

            if   data[i]['hfunc'] == 'hist':
                ax[1].hist(x=cbins, bins=bins, weights=ratio, color=c, **data[i]['style'])        
                hist_filled_error(ax=ax[1], bins=bins, cbins=cbins, y=ratio, err=ratio_errs, color=c, **data[i]['style'])

            elif data[i]['hfunc'] == 'errorbar':
                ax[1].errorbar(x=cbins, y=ratio, yerr=ratio_errs, color=c, **data[i]['style'])

            elif data[i]['hfunc'] == 'plot':
                ax[1].plot(cbins, ratio, color=c, **data[i]['style'])

                new_args = data[i]['style'].copy()
                new_args['lw'] = 0
                ax[1].fill_between(cbins, ratio-ratio_errs, ratio+ratio_errs, alpha=0.3, color=c, **new_args)
    
    # Legend
    if legend_labels != []:
        ordered_legend(ax = ax[0], order=legend_labels, **legend_properties)

    # --------------------------------------------------------------------
    # Upper figure

    # Log y-scale
    ax[0].set_yscale(yscale)

    # y-limits
    if observable['ylim'] is None:
        ylim_now = ax[0].get_ylim()
        if yscale == 'log':
            ax[0].set_ylim([bottom_count / 4, ceiling_count * 10])
        else:
            ax[0].set_ylim([0, ceiling_count * 1.5])
    else:
        ax[0].set_ylim(observables.ylim)
    # --------------------------------------------------------------------    

    return fig, ax


def change2density_labels(all_obs):
    """ Change to density ~ 1/N dN/dx [1/xdim] type label to y-axis """

    for key in all_obs.keys():
        xlabel = all_obs[key]['xlabel'].replace('$', '')
        all_obs[key]['ylabel'] = '$\\frac{1}{N} \\; ' + f'dN/d{xlabel}$'
        all_obs[key]['units']['y'] = '1'

    return all_obs


def histmc(mcdata, all_obs, density=False, scale=None, color=(0,0,1), label='none', style=hist_style_step):
    """ Over all observables of an MC sample """

    obj       = {}

    for OBS in all_obs.keys():

        # Histogram it
        counts, errs, bins, cbins = hist(x=mcdata['data'][OBS], bins=all_obs[OBS]['bins'], weights=mcdata['weights'])

        # Compute differential cross section within histogram range
        # Note that division by sum(weights) handles the histogram range integral (overflow) properly
        binscale = mcdata['xsection_pb'] / binwidth(bins) / np.sum(mcdata['weights']) 

        # Additional scale factor
        if scale is not None:
            binscale *= scale
        
        # Density integral 1 over the histogram bins
        if density:
            counts,errs = hist_to_density(counts=counts, errs=errs, bins=bins)
            binscale    = 1.0
        
        obj[OBS] = {'hdata': hobj(counts, errs, bins, cbins, binscale), 'hfunc' : 'hist', 'color': color, 'label': label, 'style' : style}

    return obj


def histhepdata(hepdata, all_obs, scale=None, density=False, MC_XS_SCALE=1E12, label='Data', style=hist_style_step):

    # Over all observables
    obj = {}

    for OBS in all_obs.keys():

        # Over all DATA files (now fixed to one)
        data_obj = []

        y        = hepdata[OBS]['y']    
        yerr     = hepdata[OBS]['y_err']
        bins     = hepdata[OBS]['bins']
        cbins    = hepdata[OBS]['x']

        binscale = hepdata[OBS]['scale'] * MC_XS_SCALE

        # Additional scale factor
        if scale is not None:
            binscale *= scale
        
        # Density integral 1 over the histogram bins
        if density:
            norm  = hepdata[OBS]['binwidth'] * y.sum()
            y    /= norm
            yerr /= norm
            binscale = 1.0

        obj[OBS] = {'hdata': hobj(y, yerr, bins, cbins, binscale), 'hfunc' : 'hist', 'color': (0,0,0), 'label': label, 'style' : style}
    
    return obj


def fuse_histograms(hist_list):
    """
    Fuse a list of count histogram objects
    """
    hsum = copy.deepcopy(hist_list[0])
    for c in range(1, len(hist_list)):
        for OBS in hist_list[0].keys():
            hsum[OBS]['hdata'] += hist_list[c][OBS]['hdata']
    
    return hsum


def test_iceplot():
    """ Visual unit tests """
    
    import pytest
    import pathlib

    pathlib.Path("./testfigs").mkdir(parents=True, exist_ok=True)


    # ------------------------------------------------------------------------
    set_global_style()


    # Synthetic input data
    r1 = np.random.randn(25000) * 0.8
    r2 = np.random.randn(25000) * 1
    r3 = np.random.randn(25000) * 1.2
    r4 = np.random.randn(25000) * 1.5


    # ------------------------------------------------------------------------
    # Mathematical definitions

    # Momentum squared
    def pt2(x):
        return np.power(x,2);


    # ------------------------------------------------------------------------
    # Observables containers

    obs_pt2 = {

    # Axis limits
    'xlim'    : (0, 1.5),
    'ylim'    : None,
    'xlabel'  : r'$p_t^2$',
    'ylabel'  : r'Counts',
    'units'   : {'x': r'GeV$^2$', 'y' : r'counts'},
    'label'   : r'Transverse momentum squared',
    'figsize' : (4, 3.75),

    # Ratio
    'ylim_ratio' : (0.7, 1.3),

    # Histogramming
    'bins'    : np.linspace(0, 1.5, 60),
    'density' : False,
     
    # Function to calculate
    'func'    : pt2
    }


    # ------------------------------------------------------------------------
    # ** Example **

    fig1, ax1 = create_axes(**obs_pt2, ratio_plot=False)
    counts, errs, bins, cbins = hist(obs_pt2['func'](r1), bins=obs_pt2['bins'], density=obs_pt2['density'])
    ax1[0].errorbar(x=cbins, y=counts, yerr=errs, color=(0,0,0), label='Data $\\alpha$', **errorbar_style)
    ax1[0].legend(frameon=False)
    fig1.savefig('./testfigs/testplot_1.pdf', bbox_inches='tight')


    # ------------------------------------------------------------------------
    # ** Example **

    fig2, ax2 = create_axes(**obs_pt2, ratio_plot=False)
    counts, errs, bins, cbins = hist(obs_pt2['func'](r1), bins=obs_pt2['bins'], density=obs_pt2['density'])
    ax2[0].hist(x=cbins, bins=bins, weights=counts, color=(0.5, 0.2, 0.1), label='Data $\\alpha$', **hist_style_step)
    ax2[0].legend(frameon=False)
    fig2.savefig('./testfigs/testplot_2.pdf', bbox_inches='tight')


    # ------------------------------------------------------------------------
    # ** Example **

    fig3, ax3 = create_axes(**obs_pt2, ratio_plot=True)

    counts1, errs, bins, cbins = hist(obs_pt2['func'](r1), bins=obs_pt2['bins'], density=obs_pt2['density'])
    ax3[0].hist(x=cbins, bins=bins, weights=counts1, color=(0,0,0), label='Data 1', **hist_style_step)

    counts2, errs, bins, cbins = hist(obs_pt2['func'](r2), bins=obs_pt2['bins'], density=obs_pt2['density'])
    ax3[0].hist(x=cbins, bins=bins, weights=counts2, color=(1,0,0), alpha=0.5, label='Data 2', **hist_style_step)

    ordered_legend(ax = ax3[0], order=['Data 1', 'Data 2'])

    # Ratio
    plot_horizontal_line(ax3[1])
    ax3[1].hist(x=cbins, bins=bins, weights=counts2 / (counts1 + 1E-30), color=(1,0,0), alpha=0.5, label='Data $\\beta$', **hist_style_step)

    fig3.savefig('./testfigs/testplot_3.pdf', bbox_inches='tight')


    # ------------------------------------------------------------------------
    # ** Example **

    data_template = {
        'data'   : None,
        'weights': None,
        'label'  : 'Data',
        'hfunc'  : 'errorbar',
        'style'  : errorbar_style,
        'obs'    : obs_pt2,
        'hdata'  : None,
        'color'  : None
    }

    # Data source <-> Observable collections
    data1 = data_template.copy() # Deep copies
    data2 = data_template.copy()
    data3 = data_template.copy()
    data4 = data_template.copy()


    data1.update({
        'data'   : r1,
        'label'  : 'Data $\\alpha$',
        'hfunc'  : 'errorbar',
        'style'  : errorbar_style,
    })
    data2.update({
        'data'   : r2,
        'label'  : 'Data $\\beta$',
        'hfunc'  : 'hist',
        'style'  : hist_style_step,
    })
    data3.update({
        'data'   : r3,
        'label'  : 'Data $\\gamma$',
        'hfunc'  : 'hist',
        'style'  : hist_style_step,
    })
    data4.update({
        'data'   : r4,
        'label'  : 'Data $\\delta$',
        'hfunc'  : 'plot',
        'style'  : plot_style,
    })

    data = [data1, data2, data3, data4]


    # Calculate histograms
    for i in range(len(data)):
        data[i]['hdata'] = hist_obj(data[i]['obs']['func'](data[i]['data']), bins=data[i]['obs']['bins'])

    # Plot it
    fig4, ax4 = superplot(data, ratio_plot=True, yscale='log')
    fig5, ax5 = superplot(data, ratio_plot=True, yscale='linear', ratio_error_plot=False)

    fig4.savefig('./testfigs/testplot_4.pdf', bbox_inches='tight')
    fig5.savefig('./testfigs/testplot_5.pdf', bbox_inches='tight')

