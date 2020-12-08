from ase.io import read, write
from gpaw import GPAW, PW, FermiDirac
import os
import sys
import re
import pathlib
import numpy as np
import pickle

import pymatgen as p
from pymatgen.io import ase
from pymatgen.io.cif import CifParser

from ase.parallel import parprint, world

from gpaw import GPAW, restart
from gpaw.unfold import Unfold, find_K_from_k
from gpaw import GPAW, FermiDirac


def get_ase(system="Li"):
    parser = CifParser("small.cif")
    lco = parser.get_structures()[0]
    if system == 'Li':
        return ase.AseAtomsAdaptor().get_atoms(lco)
    if system == "H":
        lco.replace(0, "H")
        return ase.AseAtomsAdaptor().get_atoms(lco)


def gs_calc(licoo2):
    # check if the calcualtion has been already done
    if os.path.isfile('licoo2_gs.gpw') == False:
        parprint("restart file not found :(")
        if world.rank == 0:
            pathlib.Path('outputs').mkdir(parents=True, exist_ok=True)
        calc = GPAW(mode=PW(500),
                    xc='PBE',
                    kpts=(12, 12, 12),
                    random=True,
                    occupations=FermiDirac(0.01),
                    nbands=-90,
                    txt='outputs/licoo2_gs.txt')
        licoo2.calc = calc
        licoo2.get_potential_energy()
        calc.write('licoo2_gs.gpw', 'all')


def band_calc():
    if world.rank == 0:
        pathlib.Path('outputs').mkdir(parents=True, exist_ok=True)
    calc = GPAW('licoo2_gs.gpw',
                nbands=-20,
                fixdensity=True,
                symmetry='off',
                kpts={
                    'path': licoo2.cell.bandpath().path,
                    'npoints': 120
                },
                convergence={'bands': -20},
                txt='outputs/licoo2_bands.txt')
    occasionally = 2

    class OccasionalWriter:
        def __init__(self):
            self.iter = 0

        def write(self):
            calc.write('filename.%03d.gpw' % self.iter)
            self.iter += occasionally

    calc.attach(OccasionalWriter().write, occasionally)

    calc.get_potential_energy()
    calc.write('licoo2_bands.gpw', mode='all')
    return calc


licoo2 = get_ase(system="Li")
write("lco.cif", licoo2)
parprint("Groundstate calculation started")
gs_calc(licoo2)

parprint("Done\n Band calculation started")
calc = band_calc()

parprint("Done\n Plotting bands")
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=False, emax=10.0, emin=-5.0)

with open('band_energies.pickle', 'wb') as handle:
    pickle.dump(calc.band_structure().todict(),
                handle,
                protocol=pickle.HIGHEST_PROTOCOL)
parprint("Done")
