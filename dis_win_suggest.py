import numpy as np
import pandas as pd
import os, re
import matplotlib.pyplot as plt
import sys
import argparse

class W90():
    def __init__(self, eig='EIGENVAL', path='.', win='wannier90.win', efermi=None, nbnds_excl=None, nwann=None, ndeg=1):
        '''
        Init
        '''
        self._fname = eig
        self._win = win
        # the directory containing the input file
        self._dname = path
        # if poscar is None:
        #     self.poscar = self._dname + '/POSCAR'
        self.efermi = efermi
        self.nbnds_excl = nbnds_excl
        self.nwann = nwann
        self.ndeg = ndeg # denegeracy of bands, actually only Kramers degeneracy counts

        self.read_eigenval()

        # if os.path.isfile(self._win):
        #     self.read_wannier90_win()

    def read_eigenval(self):
        '''
        Read band energies from VASP EIGENVAL file or Wannier90 .eig file.
        '''
        if self._fname[-3:] == 'eig':
            self.nspin = 1
            self.nelect = None
            data = np.loadtxt(f'{self._dname}/{self._fname}')
            self.ir_kpath = None
            self.ir_kptwt = None
            self.ir_nkpts = int(data[:, 1].max())
            self.nbnds = int(data[:, 0].max())
            self.ir_ebands = data[:, 2].reshape((self.ir_nkpts, 1, -1)).swapaxes(0, 1)
        else:
            with open(f'{self._dname}/{self._fname}') as inp:
                # read all the data.
                dat = np.array([line.strip() for line in inp if line.strip()])

            # extract the needed info in the header.
            self.nspin = int(dat[0].split()[-1])
            self.nelect, self.ir_nkpts, self.nbnds = map(int, dat[5].split())

            # remove the header
            dat = dat[6:]

            # extract the k-points info
            dump = np.array(
                [xx.split() for xx in dat[::self.nspin * self.nbnds + 1]], dtype=float)
            self.ir_kpath = dump[:self.ir_nkpts, :3]
            self.ir_kptwt = dump[:self.ir_nkpts, -1]

            # extract the ebands info
            ebands_flag = np.ones(dat.size, dtype=bool)
            ebands_flag[::self.nspin * self.nbnds + 1] = 0
            if self.nspin == 1:
                ebands = np.array([xx.split()[1] for xx in dat[ebands_flag]],
                                    dtype=float)
            else:
                ebands = np.array([xx.split()[1:3] for xx in dat[ebands_flag]],
                                    dtype=float)
            ebands.shape = (self.ir_nkpts, self.nspin, self.nbnds)
            self.ir_ebands = ebands.swapaxes(0, 1)

        self.emax = self.ir_ebands.max()
        self.emin = self.ir_ebands.min()
        self.eband_max = np.max(self.ir_ebands, axis=1)[0]
        self.eband_min = np.min(self.ir_ebands, axis=1)[0]

    def read_wannier90_win(self):
        '''
        Read parameters from wannier90.win file.
        '''
        with open(self._win) as f:
            dat = [line.strip() for line in f if line.strip()]

        for l in dat:
            r = re.match('num_wann\s*=\s*(\d+)', l)
            if r:
                self.nwann = eval(r.group(1))
            r = re.match('dis_win_max\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.win_max = eval(r.group(1))
            r = re.match('dis_froz_max\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.froz_max = eval(r.group(1))
            r = re.match('dis_froz_min\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.froz_min = eval(r.group(1))
            r = re.match('dis_win_min\s*=\s*(-?\d+\.\d+)', l)
            if r:
                self.win_min = eval(r.group(1))

    def plot_eigenval(self, erange=None, separate=False, savefig='eigenval_dis.png'):
        fig, ax = plt.subplots(figsize=(8, 6))

        def label_bar(string, height, rect):
            """Attach a text label on top of bar."""
            ax.annotate(f'{string}',
                        xy=(rect.get_x() + rect.get_width()/2, height),
                        xytext=(0, 0),  # 4 points vertical offset.
                        textcoords='offset points',
                        ha='center', va='bottom')

        if not separate:     # 不区分单独的能带绘制分布图
            idx = self.eband_min[1:] > self.eband_max[:-1]
            idx = np.nonzero(idx)[0]

            if self.nbnds_excl:
                idx = idx[idx > (self.nbnds_excl - 1)]
                idx_min, idx_max = [self.nbnds_excl]+list(idx+1), list(idx)+[self.nbnds-1]
            else:
                idx_min, idx_max = [0]+list(idx+1), list(idx)+[self.nbnds-1]
            
            eplot_min, eplot_max = self.eband_min[idx_min], self.eband_max[idx_max]

            ymin, ymax = erange if erange else (-1e4, 1e4)

            for i, (emin, emax) in enumerate(zip(eplot_min, eplot_max)):
                if emax > ymin and emin < ymax:
                    rect = ax.bar(i, emax-emin, width=0.6, bottom=emin-self.efermi, color='b')
                    label_bar(f'{idx_min[i]-self.nbnds_excl}', emax, rect[0])

        else:                   # 区分单独的能带绘制分布
            idx = np.arange(self.nbnds_excl, self.nbnds-1, self.ndeg) if self.nbnds_excl else np.arange(0, self.nbnds-1, self.ndeg)
            ymin, ymax = erange if erange else (-1e4, 1e4)
            for i in idx:
                emin, emax = self.eband_min[i], self.eband_max[i]
                if emax > ymin and emin < ymax:
                    rect = ax.bar(i, emax-emin, width=0.6, bottom=emin-self.efermi, color='b')
                    label_bar(f'{i-self.nbnds_excl}', emax, rect[0])

        ax.set(ylabel='Energy / eV')
        ax.grid()

        plt.savefig(savefig, dpi=200, bbox_inches='tight', transparent=True)
        # plt.show()

    def report_eigenval(self, erange=None, separate=False):
        print(f'EFERMI: {self.efermi: 2.6f}')
        print('--------------------------------')
        print('Band No.     EMIN        EMAX')
        print('--------------------------------')
        if not separate:     # 不区分单独的能带绘制分布图
            idx = self.eband_min[1:] > self.eband_max[:-1]
            idx = np.nonzero(idx)[0]

            if self.nbnds_excl:
                idx = idx[idx > (self.nbnds_excl - 1)]
                idx_min, idx_max = [self.nbnds_excl]+list(idx+1), list(idx)+[self.nbnds-1]
            else:
                idx_min, idx_max = [0]+list(idx+1), list(idx)+[self.nbnds-1]
            
            eplot_min, eplot_max = self.eband_min[idx_min], self.eband_max[idx_max]

            ymin, ymax = erange if erange else (-1e4, 1e4)

            for i, (emin, emax) in enumerate(zip(eplot_min, eplot_max)):
                if emax > ymin and emin < ymax:
                    print(f'{idx_min[i]:3d}~{idx_max[i]:3d}  {emin:+10.5f}  {emax:+10.5f}')

        else:                   # 区分单独的能带绘制分布
            idx = np.arange(self.nbnds_excl, self.nbnds-1, self.ndeg) if self.nbnds_excl else np.arange(0, self.nbnds-1, self.ndeg)
            ymin, ymax = erange if erange else (-1e4, 1e4)
            for i in idx:
                emin, emax = self.eband_min[i], self.eband_max[i]
                if emax > ymin and emin < ymax:
                    print(f'  {i:3d}    {emin:+10.5f}  {emax:+10.5f}')
        print('--------------------------------')
    
    def count_states(self, erange):
        emin, emax = erange
        mask = np.logical_and(self.eband_min <= emax, self.eband_max >= emin)
        return sum(mask)

    def suggest_win_max(self, emin, nwann=None, eps=4e-3):
        '''
        Lower bound of dis_win_max for given froz_min and nwann
        '''
        nwann = nwann if nwann else self.nwann
        mask = self.eband_min >= emin
        idx = np.argmax(mask) + nwann - 1
        if idx >= self.nbnds:
            # TODO handle error case for no enough states
            print(f'There is no enough states for {emin} with {nwann} WFs!')
            res = int(self.emax) + 1.
        else:
            res = self.eband_max[idx] + eps
        return res

    def suggest_froz_min(self, emax, nwann=None, eps=4e-3):
        '''
        Lower bound of froz_min for given froz_max and nwann
        '''
        nwann = nwann if nwann else self.nwann
        mask_emax = self.eband_min <= emax
        idx = np.argmin(mask_emax) - nwann - 1
        res = int(self.emin) - 1. if idx < 0 else self.eband_max[idx] + eps
        return res

    def suggest_froz_max(self, emin, nwann=None, eps=4e-3):
        '''
        Upper bound of froz_max for given froz_min and nwann
        '''
        nwann = nwann if nwann else self.nwann
        mask_emin = self.eband_max >= emin
        idx = np.argmax(mask_emin) + nwann
        res = int(self.emax) + 1. if idx >= self.nbnds else self.eband_min[idx] - eps
        return res

    def get_dis_froz_df(self, erange, eps=4e-3):
        # suggest frozen window with given energy interval
        N = self.count_states(erange)
        print(f'There are {N} states in {erange} with Fermi level at {self.efermi}.')
        emin, emax = erange
        dN = N - self.nwann
        froz_min_list, froz_max_list = [], []

        if self.nwann <= 0:
            print(f'Please input vaild number of WF, now is {self.nwann}.')
            return pd.DataFrame(columns=['dis_froz_min', 'dis_froz_max'])

        elif dN > 0:
            print('Suggest dis_froz_min & dis_froz_max as following:')
            print(f'nwann: {self.nwann}    degenercy: {self.ndeg}    Fermi: {self.efermi:12.6f}')

            for i in range(1, dN + 1, self.ndeg):
                # First get froz_max for nwann = nwann_input + i, then get froz_min for nwann = nwann_input and froz_max
                froz_max = self.suggest_froz_max(emin, nwann=self.nwann+i, eps=eps)
                froz_min = self.suggest_froz_min(froz_max, eps=eps)
                froz_max_list.append(froz_max)
                froz_min_list.append(froz_min)

            # number of missing states between `emin` and lowest `dis_froz_min`
            num_missing = self.count_states((emin, min(froz_min_list)))
            df = pd.DataFrame(zip(froz_min_list, froz_max_list), columns=['dis_froz_min', 'dis_froz_max'])

            if num_missing > 0:
                print(f'\nWANRING: There are states between given `emin`: {emin} and lowest `dis_froz_min`: {min(froz_min_list)}. Please carefully treat the suggestion of dis_froz_min / dis_froz_max and check energy range of each bands again. This situation usually happens in no-SOC system with many denegeracy point. But we still want to give you some useful energy window information with states less than number of WFs.\n')
                for i in range(self.nwann - num_missing + 1, self.nwann, 1):
                    new = pd.DataFrame({"dis_froz_min" : emin,
                                        "dis_froz_max" : self.suggest_froz_max(emin, nwann=i)}, index=[1])
                    df = df.append(new, ignore_index=True)

            return df
        else:
            return pd.DataFrame(columns=['dis_froz_min', 'dis_froz_max'])

def get_efermi(args, direct=False):
    if not direct and args.efermi:
        efermi = float(args.efemri)
    else:
        efermi_str = os.popen(f'grep fermi {args.path}/vasprun.xml').read().strip()
        m = re.match('.+ ((\-|\+)?\d+(\.\d+)?) .+', efermi_str)
        efermi = float(m.groups()[0])
    return efermi

def get_args():
    '''
    CML parser.
    '''
    parser = argparse.ArgumentParser(description='CLI Tool for W90 energy windows.', add_help=True)

    parser.add_argument('mode', help='Mode: report, plot, count, suggest')
    parser.add_argument('-i', dest='eig', action='store', type=str,
                        default='EIGENVAL',
                        help='Select wannier90.eig file or EIGENVAL file. Default: EIGENVAL')
    parser.add_argument('--path', default='.',
                        help='Default: .')
    parser.add_argument('--efermi', dest='efermi', action='store',
                        default=None,
                        help='Fermi level. Default value is generated from `vasprun.xml`.')
    parser.add_argument('-w', dest='nwann', action='store', type=int,
                        default=0,
                        help='Number of Wannier Functions. Default: 0')
    parser.add_argument('-n', dest='nbnds_excl', action='store', type=int,
                        default=0,
                        help='Number of bands excluded')
    parser.add_argument('-d', dest='ndeg', action='store', type=int,
                        default=2,
                        help='Number of degeneracy')
    parser.add_argument('-e', dest='erange', action='store', type=float,
                        default=None, nargs=2,
                        help='Energy range.')
    parser.add_argument('--separate', default=False, action="store_true",
                        help='Calculate bands not separately.')
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()

    w90 = W90(eig=args.eig,
              path=args.path,
              efermi=get_efermi(args), 
              nbnds_excl=args.nbnds_excl, 
              nwann=args.nwann, 
              ndeg=args.ndeg)

    if args.mode[0].lower() == 'p': # plot
        w90.plot_eigenval(erange=args.erange, separate=args.separate)
    elif args.mode[0].lower() == 'r': # report
        w90.report_eigenval(erange=args.erange, separate=args.separate)
    elif args.mode[0].lower() == 'c': # count
        # Count how many states inside the energy interval
        print(f'There are {w90.count_states(args.erange)} states in {args.erange}.')
    elif args.mode[0].lower() == 's': # suggest
        # suggest frozen window with given energy interval
        print(f'`dis_froz_min` and `dis_froz_max` Table:')
        df = w90.get_dis_froz_df(args.erange, eps=4e-3)
        if len(df) > 0:
            print(df)
        # dis_windows require energy window containing states larger than number of target WFs. This will also generate some constraint for dis_windows
        dis_win_max = w90.suggest_win_max(args.erange[0])
        print(f'\nLowest `dis_win_max` for {args.erange[0]}: {dis_win_max}')

    else:
        print(f'Unsupported mode: {args.mode}')