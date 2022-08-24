import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import os

# pymatgen
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen import Orbital
# spin
from pymatgen import Spin

# 计算在某一个能量区间范围内的占比
# from scipy.interpolate import interp1d
from scipy.interpolate import splrep, splev
from scipy.integrate import fixed_quad

from collections import Counter
from functools import reduce

import matplotlib
matplotlib.use('Agg')

def dos_distribute(e, dos, window):
    # dos_func = interp1d(e, dos, kind='nearest')
    dos_cubic_fitting = splrep(e, dos, s=0)
    dos_func = lambda x: splev(x, dos_cubic_fitting, der=0)
    res = fixed_quad(dos_func, window[0], window[-1])
    return res[0]

def gen_dos_df(dos_data_total, left, right, orb_names=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']):
    structure = dos_data_total.structure
    num_sites = structure.num_sites
    dos_df = pd.DataFrame(columns=["species", "structure_id", "orb_id", "orb_name", "key_string", "dos"])
    for i in range(num_sites):
        for j in range(len(orb_names)):
            site, orb = structure[i], Orbital(j)
            dos = dos_data_total.get_site_orbital_dos(site, orb)
            dis = dos_distribute(dos.energies - dos.efermi, dos.densities[Spin.up], [left, right])
            key = site.species_string + '_' + str(i) + '_' + orb_names[j]
            new = pd.DataFrame({"species": site.species_string,
                                "structure_id": i,
                                "orb_id": j,
                                "orb_name": orb_names[j],
                                "key_string": key,
                                "dos": dis}, index=[1])
            dos_df = dos_df.append(new, ignore_index=True)
    return dos_df

def plot_dos_dis(dos_df, pct=1, selected=set(), path='.', filename='dos_analysis.png', colors=['brown', 'orange']):
    threshold = pct / 100 * dos_df['dos'].max()
    df = dos_df[dos_df['dos'] > threshold]
    df = df.sort_values(by='dos', ascending=False)
    # color = [colors[0] if row['key_string'] in selected else colors[-1]
    #              for _, row in df.iterrows()]
    color = [colors[0] if row in selected else colors[-1] for row in df['key_string']]
    print(f'Plot with selected orbitals in `{colors[0]}` and not selected orbitals in `{colors[-1]}`.')
    print(f'{len(selected)} orbitals selected')
    mask = np.array(color) == colors[0]
    print(f'Plot {len(df)} orbitals with {sum(mask)} selected')
    df.plot.barh(x="key_string", y="dos", color=color, figsize=(8, 20))
    plt.savefig(filename, dpi=200, bbox_inches='tight', transparent=True)

# 返回最终选择的 轨道数量 和 轨道列表, 以 pandas.DataFrame 的形式
# 其中 site 如果为 -1, 意味着选择该元素的全体轨道
# pick_rate 指选择的轨道成分为 0.1*max_dos
def dos_analysis_df(dos_df, pick_rate=0.1):
    # 选择出满足条件的轨道
    threshold = pick_rate * dos_df['dos'].max()
    df = dos_df[dos_df['dos'] > threshold]

    # print(df)
    n_orb = len(df)

    # 存储筛选结果数据表, 分为两个部分, site 和 orb
    res_df = pd.DataFrame(columns=["species", "site", "orb"])

    # 判断一个格点上是否同一类型的轨道都被选择了
    d_orb = ['dxy', 'dyz', 'dxz', 'dx2', 'dz2']
    p_orb = ['px', 'py', 'pz']
    s_orb = ['s']

    # 循环遍历所有元素进行检查
    for species in df.species.unique():
        species_sites = df[df["species"]==species].structure_id.unique()

        # 如果轨道集合中包含 spd 轨道或者和其它轨道的组成, 需要进行简化
        # 比如将 [s, px, py, pz, dx2] 简化为 [s, p, dx2]
        # 检查各格点上的轨道是否可以统一为 spd
        for idx in species_sites:
            site_df = df[df["structure_id"]==idx]
            site_orbs = site_df.orb_name.unique()
            for orb in ["s", "p", "d"]:
                if orb == "s":
                    orb_list = s_orb
                elif orb == "p":
                    orb_list = p_orb
                elif orb == "d":
                    orb_list = d_orb
                if set(orb_list).issubset(set(site_orbs)):
                    site_orbs = list(set(site_orbs) - set(orb_list)) + [orb]
            new = pd.DataFrame({"species": species,
                                "site": idx,
                                "orb": [site_orbs]}, index=[1])
            res_df = res_df.append(new, ignore_index=True)

    # display(res_df)

    # 合并具有同样元素同一轨道的格点
    for species in res_df.species.unique():
        # 获取这个元素对应的所有 structure 的 id
        species_site_id = dos_df[dos_df["species"]==species].structure_id.unique()
        species_site = res_df[res_df["species"]==species].site.unique()
        species_site_orbs = res_df[res_df["species"]==species].orb
        if len(species_site_id) == len(species_site):
            ref_orb = species_site_orbs.to_list()[0]
            is_orb_united = reduce(lambda x, y: x and y, [Counter(ref_orb)==Counter(orb) for orb in species_site_orbs])
            if is_orb_united:
                # 删除所有具有这一元素的数据
                res_df = res_df[res_df['species']!=species]
                # 添加一行合并后的数据, -1 表示合并后的格点
                new = pd.DataFrame({"species": species,
                                    "site": -1,
                                    "orb": [list(ref_orb)]})
                res_df = res_df.append(new, ignore_index=True)

    # display(res_df)
    return n_orb, res_df, df

def orb_string(orb_list):
    orb_string_list = []
    for o in orb_list:
        if o == "s":
            orb_string_list.append('l=0')
        elif o == "p":
            orb_string_list.append('l=1')
        elif o == "d":
            orb_string_list.append('l=2')
        else:
            orb_string_list.append(o)
    return ';'.join([s for s in orb_string_list])

# 生成 wannier90 对应的 projection 选择格式
def w90_string(res_df, structure):
    w90_string_list = []
    for _, row in res_df.iterrows():
        if row['site'] == -1:
            w90_string_list.append(f"{row['species']}:{orb_string(row['orb'])}")
        else:
            site_string = ','.join([f"{i:.6f}" for i in structure[row['site']].frac_coords])
            w90_string_list.append(f"c={site_string}:{orb_string(row['orb'])}")
    print('\n'.join([s for s in w90_string_list]))

def dos_given_orb(dos_data_total, left, right, full_res_df):
    structure = dos_data_total.structure
    res = 0
    for _, row in full_res_df.iterrows():
        site = structure[row["structure_id"]]
        orb = Orbital(row["orb_id"])
        dos = dos_data_total.get_site_orbital_dos(site, orb)
        dis = dos_distribute(dos.energies - dos.efermi, dos.densities[Spin.up], [left, right])
        res += dis
    return res

# from dis_win_suggest file
def get_dis_win_table(band_data, num_wann, n_excl=0, epsilon=0.02):
    if n_excl < 0:
        raise ValueError('Wrong Input!')

    nbands = band_data.shape[0]
    n_choices = nbands-num_wann-n_excl
    dis_win_table = np.zeros((n_choices, 4))

    # 获取所有可能的能带选择
    for idx in range(n_choices):
        dis_win_max  = band_data[n_excl + idx + num_wann - 1].max() + epsilon
        dis_froz_max = band_data[n_excl + idx + num_wann].min() - epsilon
        dis_win_min  = band_data[n_excl + idx].min() - epsilon
        dis_froz_min = band_data[n_excl + idx - 1].max() + epsilon if n_excl + idx > 0 else dis_win_min

        dis_win_table[idx] = np.array((dis_win_max, dis_froz_max, dis_froz_min, dis_win_min))

    return dis_win_table

def print_dis_win_table(dis_win_table):
    # output the result
    # print('dis_win_max  : larger than  {0}'.format(dis_win_max))
    # print('dis_froz_max : smaller than {0}'.format(dis_froz_max))
    # print('dis_froz_min : smaller than {0}'.format(dis_froz_min))
    # print('dis_win_min  : smaller than {0}'.format(dis_win_min))
    print('dis_win_max  dis_froz_max  dis_froz_min  dis_win_min')
    for line in dis_win_table:
        dis_win_max, dis_froz_max, dis_froz_min, dis_win_min = line
        print(f'{dis_win_max:4.2f}  {dis_froz_max:4.2f}  {dis_froz_min:4.2f}  {dis_win_min:4.2f}')

def select_str2list(s, orb_names=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']):
    def num2list(s):
        l, r = map(int, s.split('-')) if '-' in s else [int(s)] * 2
        return list(range(l, r + 1, 1))

    if len(s) == 0:
        return []

    res = []
    for l in s.split(';'):
        spec, pos, orb = l.split(',')
        pos, orb = num2list(pos), num2list(orb)
        for p in pos:
            for o in orb:
                res.append(f'{spec}_{p}_{orb_names[o]}')
    return set(res)

def kpath(kpoints, delimeter=None):
    delimeter = delimeter if delimeter else '!'
    with open(kpoints, 'r') as f:
        lines = f.readlines()
    lines    = [l.strip() for l in lines if len(l.strip())!=0 ][4:]
    labels   = [s.split(delimeter)[-1].strip() for s in lines]
    nums_str = [s.split(delimeter)[0].strip() for s in lines]
    nums     = [[eval(n) for n in s.split()] for s in nums_str]
    nice_string = lambda p, label: f"{label:5s} {p[0]: 02.5f}    {p[1]: 02.5f}    {p[2]: 02.5f}"
    for i in range(0, len(labels), 2):
        s1 = nice_string(nums[i], labels[i])
        s2 = nice_string(nums[i+1], labels[i+1])
        print(f'{s1}  {s2}')

def export_vasp_band(path):
    def export2dat(kk, EE, filename):
        lines = []
        nbnds, nk = EE.shape
        for i in range(nbnds):
            for j in range(nk):
                lines.append(f'{kk[j]: 12.6f} {EE[i, j]: 12.6f}\n')
            lines.append('\n')
        with open(filename, 'w') as f:
            f.writelines(lines)

    print(f"Generate p4vasp format bnd.dat file from {path}/vasprun.xml")
    run = Vasprun(f'{path}/vasprun.xml')
    bands = run.get_band_structure(f'{path}/KPOINTS', line_mode=True)
    nspin = len(bands.bands.keys())
    kk = bands.distance
    if nspin == 1:
        print('export band data to `bnd.dat`.')
        export2dat(kk, bands.bands[Spin.up], f'{path}/bnd.dat')
    elif nspin == 2:
        print('NSPIN = 2')
        print('export band data to `bnd_up.dat` and `bnd.down.dat` seperatly.')
        export2dat(kk, bands.bands[Spin.up], f'{path}/bnd_up.dat')
        export2dat(kk, bands.bands[Spin.down], f'{path}/bnd_down.dat')

def template(flag):
    if flag == 'wann':
        print('###############\n' \
              '#     W90     #\n' \
              '###############\n' \
              'dis_win_max       =  0.0\n' \
              'dis_froz_max      =  0.0\n' \
              'dis_froz_min      = -5.0\n' \
              'dis_win_min       = -5.0\n' \
              '\n'                         \
              'num_iter          = 1000\n' \
              'num_print_cycles  =   40\n' \
              'dis_num_iter      = 5000\n' \
              'dis_mix_ratio     =  1.0\n')
    elif flag == 'band':
        print('###############\n' \
              '#  Band Plot  #\n' \
              '###############   # restart = plot\n' \
              'write_hr          = true\n' \
              'bands_plot        = true\n' \
              'bands_num_points  = 151 \n' \
              'bands_plot_format = gnuplot\n' \
              '\n' \
              'Begin Kpoint_Path\n' \
              'End Kpoint_Path\n')
    else:
        print('Unknown flag: ', flag)

def get_args():
    '''
    CML parser.
    '''
    parser = argparse.ArgumentParser(description="Pre analysis before Wannier90 Interpolation.\n    dos: python pre_w90_tool.py dos --plot -e -4 7 --extra 'Bi,4-7,0-3;F,8-23,1-3'",     add_help=True)

    parser.add_argument('mode', help='Mode: kpath, band, template, dos')
    parser.add_argument('--path', default='.',
                        help='Default: .')
    parser.add_argument('--no-soc', action='store_true', default=False,
                        help='with SOC or not. Default: False')
    parser.add_argument('--pick', action='store', type=float, default=0.1,
                        help='Minimum of selected orbital / max single orbital. default: 0.1')
    parser.add_argument('-e', dest='erange', action='store', type=float,
                        default=None, nargs=2,
                        help='Energy range.')
    parser.add_argument('--plot', default=False, action="store_true",
                        help='plot the dos distribution')
    parser.add_argument('--extra', action='store', type=str, default='',
                        help='Extra input.')
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()

    if args.mode[0].lower() == 'k': # generate kpath
        if 'KPOINTS' in os.listdir(args.path):
            print(f'Generate Kpoint_Path from {args.path}/KPOINTS:\n')
            kpath(f'{args.path}/KPOINTS', delimeter=args.extra)
        else:
            print(f'There is no KPOINTS file in {args.path}. Please check your input!')
    elif args.mode[0].lower() == 't':   # print W90 parameter template
        if len(args.extra) == 0:
            template('wann')
            template('band')
        else:
            template(args.extra)
    elif args.mode[0].lower() == 'b':   # generate p4vasp-like band data file
        export_vasp_band(args.path)
    elif args.mode[0].lower() == 'd':   # DOS Analysis
        left, right = args.erange #-4, 8
        pick_rate = args.pick #0 .08

        vasprun = Vasprun(f"{args.path}/vasprun.xml")
        print(f"Reading vasprun.xml file from `{args.path}/vasprun.xml` for dos analysis")
        dos_data_total = vasprun.complete_dos       # get dos data
        structure = dos_data_total.structure
        dos_df = gen_dos_df(dos_data_total, left, right)
        print(f"\nCalculated DOS Energy Range: {left}, {right}")
        print(dos_df[:30])

        norb, simple_res_df, full_res_df = dos_analysis_df(dos_df, pick_rate=pick_rate)
        nwann = norb if args.no_soc else 2 * norb
        print(f"\nNumber of Selected Orbitals: {norb}")
        print(f"\nNumber of Selected WFs: {nwann}")
        print("\nSelected Orbitals: ")
        print(simple_res_df)

        if args.plot:
            selected = select_str2list(args.extra)
            plot_dos_dis(dos_df, selected=selected, filename=f'{args.path}/dos_analysis.png')

        # Manually Input the Orbital Selection
        # ------------------------------------
        # orb id order: s = 0 py = 1 pz = 2 px = 3 dxy = 4 dyz = 5 dz2 = 6 dxz = 7 dx2 = 8
        # see orb id in pymatgen Orbital https://pymatgen.org/pymatgen.electronic_structure.core.html#pymatgen.electronic_structure.core.Orbital
        # data = [['Ga', 0, 6],
        #         ['As', 1, 1],
        #         ['As', 1, 2],
        #         ['As', 1, 3]]
        # full_res_df = pd.DataFrame(data, columns = ['species', 'structure_id', 'orb_id'])
        # print("\nManually Selected Orbitals: ")
        # print(full_res_df)
        # n_orb = 4

        print("\nWannier90 Projection:")
        w90_string(simple_res_df, structure)

        # with given dis_froz_min & dis_froz_max value, calculate the dos distribution
        # parse band data

        # v = Vasprun(band_path + '/' + band_vasprun_file)
        # bands = v.get_band_structure(band_path + '/' + 'KPOINTS', line_mode=True)
        # band_data = bands.bands[Spin.up]

        # num_wann = 2 * n_orb if is_soc else n_orb
        # # n_excl = 0
        # dis_win_table = get_dis_win_table(band_data, num_wann, n_excl=0, epsilon=0.02)
        # # print_dis_win_table(dis_win_table)
        # # 存储为 pandas.DataFrame
        # dis_win_df = pd.DataFrame(columns=["dis_win_max", "dis_froz_max", "dis_froz_min", "dis_win_min", "dos"])
        # for line in dis_win_table:
        #     dis_win_max, dis_froz_max, dis_froz_min, dis_win_min = line
        #     e_dos = dos_given_orb(dos_data_total, dis_froz_min, dis_froz_max, full_res_df)
        #     new = pd.DataFrame({"dis_win_max": dis_win_max,
        #                         "dis_froz_max": dis_froz_max,
        #                         "dis_froz_min": dis_froz_min,
        #                         "dis_win_min": dis_win_min,
        #                         "dos": e_dos}, index=[1])
        #     dis_win_df = dis_win_df.append(new, ignore_index=True)

        # # 展示前 70% 的数据
        # threshold = 0.7 * dis_win_df['dos'].max()
        # print(dis_win_df[dis_win_df['dos'] > threshold])
