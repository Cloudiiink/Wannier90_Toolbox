#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import numpy as np

import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, Orbital

from scipy.ndimage import gaussian_filter1d as gaussian

import matplotlib
matplotlib.use('Agg')

def hex2rgb(hexcode):
    return tuple([int(hexcode[2*i+1:2*i+3], 16) for i in range(3)])

def orbstring2py(orb):
    orb_names = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2']
    idx = 0
    for i, orb_name in enumerate(orb_names):
        if orb_name == orb:
            idx = i
    return Orbital(idx)

if __name__ == "__main__":
    # TODO Need to polish
    # set up matplotlib plot
    # ----------------------
    plt.rcParams['font.family'] = "Open Sans"
    plt.rcParams['font.weight'] = "regular"
    plt.rcParams['font.size'] = 24

    fig, ax = plt.subplots(figsize=(10, 5))

    # set ylim for the plot
    # ---------------------
    emin  = -6
    emax  =  10
    sigma =  0.02
    ymax  =  5

    ax.set_xlim(emin, emax)

    # Basic Orbital Information
    # -------------------------
    # test
    is_index_given = True
    el_orb_dict = { "0" : ["p"],        # separate with comma, e.g. "1,2,3,5": ["pz"]
                    "1" : ["p"]}
    labels = ["Ga-p", "As-p"]
    base_colors = ['#E0AC09', '#F72916']
    #
    # is_index_given = False
    # el_orb_dict = {"Ga": ["s", "p"], "As":["p"]}
    # base_colors = ['#E0AC09', '#F72916', '#6938E0', '#1C00BF', '#00FF00']
    # # color from https://color.adobe.com/zh/create/color-wheel

    # Density of states
    # -----------------

    # read dos data
    # -------------
    dospath = '.'
    dosrun = Vasprun(dospath + "/vasprun.xml")

    ee = dosrun.tdos.energies
    # get element_orbital_dos from site_orbital_dos
    spec_orb_dos = {}
    structure = dosrun.structures[0]
    for species, orbs in el_orb_dict.items():
        orb_dos = {}
        for orb in orbs:
            orb_py = orbstring2py(orb)
            densities = np.zeros_like(ee)

            if is_index_given == False:
                for site in structure.sites:
                    if site.species_string == species:
                        if orb == 'p':
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('px')][Spin.up]
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('py')][Spin.up]
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('pz')][Spin.up]
                        elif orb == 'd':
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('dxy')][Spin.up]
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('dyz')][Spin.up]
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('dz2')][Spin.up]
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('dxz')][Spin.up]
                            densities += dosrun.complete_dos.pdos[site][orbstring2py('dx2')][Spin.up]
                        else:
                            densities += dosrun.complete_dos.pdos[site][orb_py][Spin.up]
            else:
                site_idx_list = [int(i) for i in species.split(',')]
                for site_idx in site_idx_list:
                    site = structure.sites[site_idx]
                    if orb == 'p':
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('px')][Spin.up]
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('py')][Spin.up]
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('pz')][Spin.up]
                    elif orb == 'd':
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('dxy')][Spin.up]
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('dyz')][Spin.up]
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('dz2')][Spin.up]
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('dxz')][Spin.up]
                        densities += dosrun.complete_dos.pdos[site][orbstring2py('dx2')][Spin.up]
                    else:
                        densities += dosrun.complete_dos.pdos[site][orb_py][Spin.up]

            orb_dos[orb] = densities
        spec_orb_dos[species] = orb_dos

    # ax.set_xticklabels([])
    # ax.grid()
    ax.set_ylim(1e-6, ymax)
    # ax.set_xticklabels([])
    ax.vlines(x=dosrun.efermi, ymin=0, ymax=ymax, color="grey", linestyles="dashed", lw=1)
    ax.set_ylabel("DOS / a.u.")
    ax.set_xlabel("Energy / eV")
    #ax.set_title("Density of States")

    # spd contribution
    color_idx = 0
    for species, orbs_dos in spec_orb_dos.items():
        for orb, dos in orbs_dos.items():
            # if gaussian
            dos = gaussian(dos, sigma*len(ee)/(ee.max()-ee.min()))
            if is_index_given == False:
                ax.plot(ee,
                        dos,
                        color=base_colors[color_idx],
                        label=f'{species}-{orb}', lw=2, zorder=30)
            else:
                ax.plot(ee,
                        dos,
                        color=base_colors[color_idx],
                        label=labels[color_idx], lw=2, zorder=30)
            color_idx += 1

    # total dos
    dos = dosrun.tdos.densities[Spin.up]
    # if gaussian
    dos = gaussian(dos, sigma*len(ee)/(ee.max()-ee.min()))
    ax.fill_between(ee,
                    dos,
                    color=(0.7, 0.7, 0.7),
                    facecolor=(0.7, 0.7, 0.7))

    ax.plot(ee,
            dos,
            color=(0.3, 0.3, 0.3),
            label="total", zorder=10)

    # plot format style
    # -----------------
    ax.tick_params(axis='both', direction='in', labelsize=20)
    legend = ax.legend(fancybox=False,
                        shadow=False,
                        facecolor='white',
                        edgecolor='black',
                        loc='upper right', 
                        framealpha=1.0, 
                        frameon=True, #False,
                        prop={'size': 20})
    legend.get_frame().set_linewidth(0)

    # plt.show()
    plt.savefig(sys.argv[0].strip(".py") + ".png", format="png", bbox_inches='tight', dpi=150)
