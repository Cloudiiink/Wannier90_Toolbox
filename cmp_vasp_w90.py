#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib as mpl
import logging
from scipy import interpolate
from scipy.signal import savgol_filter
mpl.use("Agg")

logging.basicConfig(level = logging.INFO)
logger = logging.getLogger("cmp")

def parse_dat(datfile):
    data = pd.read_csv(datfile, header=None, sep=r'\s+', comment='#')
    k_idxs = np.where(data[0]==np.min(data[0]))[0]
    nk = k_idxs[1] - k_idxs[0]
    kk = np.array(data[0])[:nk]
    EE = np.array(data[1]).reshape(-1, nk)
    return kk, EE

# Used in wannier_fit evaluation
def gaussian(x, mid=0, width=3):
    mu, sig = mid, width
    return np.exp(-np.power((x - mu)/sig, 2)/2)
    # return 1/(np.sqrt(2*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2)/2)

# Used in wannier_fit evaluation
def unit(x, mid=0, width=3):
    res = np.abs(x-mid) < width
    return res.astype(float)

def plot_cmp_vasp_w90(vkk, vee, wkk, wee, 
                      ylim=None,
                      efermi=0,
                      font="Open Sans", size=18):
    w2v_ratio = np.max(vkk) / np.max(wkk)   # Theoretically, it should be 2 * pi

    # general options for plot
    plt.rcParams['font.family'] = font
    plt.rcParams['font.weight'] = "regular"
    plt.rcParams['font.size']   = size

    # plot
    fig, ax = plt.subplots(figsize=(8, 6))

    for idx, EE in enumerate(vee):
        ax.plot(vkk / w2v_ratio,
                EE - efermi,
                'k-',
                label='VASP' if idx==0 else None)
    for idx, EE in enumerate(wee):
        ax.plot(wkk[::1],
                EE[::1] - efermi,
                'r--',
                label='Wannier90' if idx==0 else None)

    # kpoints labels: from `wannier90_band.labelinfo.dat` file
    with open('wannier90_band.labelinfo.dat', 'r') as f:
        lines = f.readlines()
    label = [l.split()[0] for l in lines]
    logger.info(f"k label: {label}")
    k_node = [eval(l.split()[2]) for l in lines]
    for i, lab in enumerate(label):
        if lab[-1].isdigit():
            label[i] = lab[:-1] + r'$_' + lab[-1] + r'$'
   
    ax.set_xlim(k_node[0], k_node[-1])
    # put tickmarks and labels at node positions
    ax.set_xticks(k_node)
    ax.set_xticklabels(label)
    # add vertical lines at node positions
    ax.grid(axis='x')
    ax.hlines(y=0, xmin=k_node[0], xmax= k_node[-1], color="grey", linestyles="dashed", lw=0.5)
    ax.set_ylabel("Energy / eV")

    # get plot bound
    if not ylim:
        ylim = (wee.min()-efermi-1, wee.max()-efermi+1)
    ax.set_ylim(ylim)
    logger.info(f"Energy range: {ylim}")

    # legend
    # ------
    ax.legend(fancybox=False,
              shadow=False,
              frameon=True, #False,
              framealpha=1.0,
              facecolor='white',
              edgecolor='black',
              loc='upper right',
              prop={'size': 14})

    # plt.show()
    plt.savefig(output_figure,  bbox_inches='tight', transparent=True, dpi=300)
    logger.info(f"Output figure: {output_figure}")

def evaluate_cmp_vasp_w90(vkk, vee, wkk, wee, kernel='unit', mid=0, width=3):

    if kernel[0].lower() == 'u':
        kernel = lambda x: unit(x, mid=mid, width=width)
    elif kernel[0].lower() == 'g':
        kernel = lambda x: gaussian(x, mid=mid, width=width)

    nbnds, _ = wee.shape  # num of bands in wannier90
    w2v_ratio = np.max(vkk) / np.max(wkk)   # Theoretically, it should be 2 * pi

    ve, we = vee[:, 0], wee[:, 0]
    Nv, Nw = len(ve), len(we)
    diff = np.array([np.sum(np.abs(ve[i:i+Nw] - we)) for i in range(Nv-Nw+1)])
    nbnds_excl = np.argmin(diff)
    logger.info(f"nbnds_excl: {nbnds_excl}")

    dEs, wgts = [], []

    # mask of VASP data
    diff_vkk = vkk[1:] - vkk[:-1]
    vmask = [True] + list(np.logical_not(diff_vkk < 1e-7))
    # mask of W90 data
    diff_wkk = wkk[1:] - wkk[:-1]
    wmask = list(np.logical_not(diff_wkk < 1e-7)) + [True]

    for i in range(nbnds):
        # Interpolation need to remove the duplicates
        # REF: https://stackoverflow.com/questions/12054060/scipys-splrep-splev-for-python-interpolation-returns-nan

        tck = interpolate.splrep(wkk[wmask], wee[i][wmask])
        fit_wee_i = interpolate.splev(vkk[vmask] / w2v_ratio, tck)
        vee_i = vee[i + nbnds_excl][vmask]

        dEi = fit_wee_i - vee_i
        # Using filter to smooth the data and remove
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.savgol_filter.html
        dEi =  np.abs(savgol_filter(dEi, 9, 2))

        # ! AVERAGE DISTANCE
        # dEs.append(np.sum(kernel(vee[i + nbnds_excl][vmask]) * dEi) / len(dEi) * 1000)
        # wgts.append(np.sum(kernel(vee[i + nbnds_excl][vmask])) / len(dEi))

        # ! MAX DISTANCE
        dEs.append(np.max(kernel(vee_i) * dEi) * 1000)
        wgts.append(np.max(kernel(vee_i)))

    dEs, wgts = np.array(dEs), np.array(wgts)
    return dEs, wgts

def show_vasp_w90_diff(dEs, wgts):
    import gnuplotlib as gp

    # Normalization: dE = dEs / wgts
    N = len(dEs)
    logger.info('=== MAX DIFF of VASP vs W90 with each bands (meV) ===')
    gp.plot(np.linspace(1, N, N), dEs,
            _with    = 'lines',
            terminal = 'dumb 100, 30',
            unset    = 'grid')
    logger.info(f'Average dE (meV): {np.sum(dEs) / np.sum(wgts)}')

def show_spreading(path):
    import gnuplotlib as gp

    if 'wannier90.wout' in os.listdir(path):
        conv_str = os.popen(f'grep CONV {path}/wannier90.wout').read().split('\n')
        if len(conv_str) > 4:
            conv = conv_str[3:-1]
            spread = np.array([eval(re.split('\s+', c)[4]) for c in conv])

            x = np.linspace(0, len(spread)-1, len(spread))
            spread_min = spread.min()
            gp.plot( x, spread,
                    _with    = 'lines',
                    _yrange  = [spread_min * 0.7, spread_min * 2],
                    terminal = 'dumb 100, 30',
                    unset    = 'grid')
            idx = np.argmin(spread)
            logger.info(f'MIN_NUM_ITER: {idx}   SPREAD: {spread[idx]}')
    else:
        logger.info(f"There if no `{path}/wannier90.wout`.")

def show_all_fonts():
    # Ref: python - How to get a list of all the fonts currently available for Matplotlib? - Stack Overflow
    # https://stackoverflow.com/questions/8753835/how-to-get-a-list-of-all-the-fonts-currently-available-for-matplotlib
    from matplotlib import font_manager as fm
    fpaths = fm.findSystemFonts()
    family_name = set([fm.get_font(i).family_name for i in fpaths])
    logger.info(sorted(family_name))

def get_efermi(args):
    if args.efermi:
        efermi = float(args.efermi)
    else:
        efermi_str = os.popen(f'grep fermi {args.path}/vasprun.xml').read().strip()
        m = re.match('.+ ((\-|\+)?\d+(\.\d+)?) .+', efermi_str)
        efermi = float(m.groups()[0])
    return efermi

def get_args():
    '''
    CML parser.
    '''
    parser = argparse.ArgumentParser(description='Comparison between VASP band and Wannier90 band. `bnd.dat` for VASP band data in p4vasp format and `wannier90_band.dat`, `wannier90_band.labelinfo.dat`, and `wannier90.wout` are required for plotting and analysis.')
    parser.add_argument('name',
                        help='name of system')
    parser.add_argument('--efermi', default=None, 
                        help='Fermi level. Default value is generated from `vasprun.xml`.')
    parser.add_argument('--path', default='.',
                        help='Default: .')
    parser.add_argument('--vasp', default='bnd.dat',
                        help="location of VASP band file in p4vasp format. Default: bnd.dat")
    parser.add_argument('--ylim', default=None, nargs=2, type=float,
                        help="Energy bound for plot. Default: [E_w90.min - 1, E_w90.max + 1]")
    parser.add_argument('--kernel', default='unit,2,5',
                        help="kernel function for evaluating diff: type, middle, width. There are two type of kernel function: `unit` and `gaussian`. Defalut: unit,2,5")
    parser.add_argument('--show-fonts', default=False, action="store_true",
                         help="Show all availabel font families can be used in `rcParams`")
    parser.add_argument('--fontfamily', default='Open Sans',
                        help="Set font family manually. Default: Open Sans")
    parser.add_argument('--fontsize', default=18, type=int,
                        help="Set font size manually. Default: 18")
    parser.add_argument("--no-spread", default=False, action="store_true",
                         help="Don't plot spreading")
    parser.add_argument("--no-quality", default=False, action="store_true",
                         help="Don't show quality of fitting")
    parser.add_argument("--quiet", default=False, action="store_true",
                         help="Equal to --no-spreading --no-quality")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()

    name = args.name
    output_figure = f'{args.path}/{name}_VASP_W90_cmp.png'
    efermi = get_efermi(args)
    
    logger.info(f'Reading Data from {args.path}/{args.vasp}')
    vkk, vee = parse_dat(f'{args.path}/{args.vasp}')
    wkk, wee = parse_dat(f'{args.path}/wannier90_band.dat')

    plot_cmp_vasp_w90(vkk, vee, wkk, wee, 
                      ylim=args.ylim,
                      efermi=efermi,
                      font=args.fontfamily, size=args.fontsize)

    if not args.no_quality and not args.quiet:
        logger.info('Evaluating Band Quality:')
        l = args.kernel.split(',')
        kernel, mid, width = l[0], float(l[1]), float(l[2])
        dEs, wgts = evaluate_cmp_vasp_w90(vkk, vee, wkk, wee,
                                          kernel=kernel, mid=mid, width=width)
        show_vasp_w90_diff(dEs, wgts)

    if not args.no_spread and not args.quiet:
        logger.info('Show spreading convergence:')
        show_spreading(args.path)

    if args.show_fonts:
        show_all_fonts()
