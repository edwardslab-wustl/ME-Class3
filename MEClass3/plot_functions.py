
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

def pull_up_dn_bounds(param, features, args):
    dn = -param.region_size
    up = param.region_size
    if features == 'tss':
        if args.tss_upperBound < up:
            up = args.tss_upperBound
        if args.tss_lowerBound > dn:
            dn = args.tss_lowerBound
    elif features == 'enh':
        if args.enh_upperBound < up:
            up = args.enh_upperBound
        if args.enh_lowerBound > dn:
            dn = args.enh_lowerBound
    if param.lowerBound !=0 or param.upperBound != 0:
        if param.upperBound < up:
            up = param.uppperBound
        if param.lowerBound > dn:
            dn = param.lowerBound
    return dn, up

def setup_x_axis(figure, axis, param, anno_type, args):
    x_range = []
    if param.lowerBound !=0 or param.upperBound !=0:
        x_range = [param.lowerBound,param.upperBound]
    else:
        x_range = [-param.region_size,param.region_size]
    figure.xlim(x_range)
    if anno_type == 'tss':
        feat_label = "TSS"
        x_tick_minorLocator = MultipleLocator(1000)
        axis.xaxis.set_minor_locator(x_tick_minorLocator)
    elif anno_type == 'enh':
        feat_label = "Enhancer"
        plt.legend(loc='right', title='',bbox_to_anchor=(1.4,0.5),handlelength=1,handletextpad=0.5,frameon=False)
        x_tick_minorLocator = MultipleLocator(100)
        axis.xaxis.set_minor_locator(x_tick_minorLocator)
    else:
        feat_label = anno_type
        x_tick_minorLocator = MultipleLocator(100)
        axis.xaxis.set_minor_locator(x_tick_minorLocator)
    axis.set_xlabel(f"Distance to {feat_label} (bp)")
    figure.xticks((x_range[0],0,x_range[1]))
    return