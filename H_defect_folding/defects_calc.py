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

# ----------------- Defect calculation-----------------


def get_ase(n=3, m=1):
    parser = CifParser("small.cif")
    lco = parser.get_structures()[0]
    lco_sc = lco.copy()
    n = 3
    sc = np.array([[n, 0, 0], [0, n, 0], [0, 0, n]])
    lco_sc.make_supercell(sc)
    li = np.array(
        [i if j.species_string == "Li" else -1 for i, j in enumerate(lco_sc)])
    li = li[li >= 0]
    np.random.seed(4)
    index = np.random.choice(li.shape[0], m, replace=False)
    for i in index:
        lco_sc.replace(li[i], "H")
    percentage_li = lco_sc.composition.get_atomic_fraction(
        "Li") / lco_sc.composition.get_atomic_fraction("Co")
    # lco_sc.to("cif","lco_sc.cif")
    return ase.AseAtomsAdaptor().get_atoms(lco_sc), percentage_li


def gs_calc(licoo2):
    # check if the calcualtion has been already done
    if os.path.isfile('licoo2_gs.gpw') == False:
        parprint("restart file not found :(")
        if world.rank == 0:
            pathlib.Path('outputs').mkdir(parents=True, exist_ok=True)
        calc = GPAW(mode='lcao',
                    basis='dzp',
                    xc='PBE',
                    kpts=(3, 3, 3),
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
                nbands=-90,
                fixdensity=True,
                symmetry='off',
                kpts={
                    'path': licoo2.cell.bandpath().path,
                    'npoints': 60
                },
                convergence={'bands': -90},
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


# ----- GS calculation

num_defects = int(sys.argv[1])
licoo2, percentage = get_ase(m=num_defects)
write("lco_defect.cif", licoo2)
parprint("percentage of defect concentration = {:.3f}".format(percentage))
gs_calc(licoo2)

# ----- band calculation

parprint("Running bands calculation .....")
calc = band_calc()

# ----band structure plot
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=False, emax=10.0, emin=-5.0)

with open('band_energies_{}.pickle'.format(percentage), 'wb') as handle:
    pickle.dump(calc.band_structure().todict(),
                handle,
                protocol=pickle.HIGHEST_PROTOCOL)

# ------------------------------Band renormalization

a, calc = restart('licoo2_gs.gpw')
n = 3
M = np.array([[n, 0, 0], [0, n, 0], [0, 0, n]])

bp = a.cell.bandpath("BZG", npoints=90)
x, X, name = bp.get_linear_kpoint_axis()

Kpts = []
for k in bp.kpts:
    K = find_K_from_k(k, M)[0]
    Kpts.append(K)


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


for i, kpts in enumerate(chunks(bp.kpts, 25)):
    Kpts = []
    for k in kpts:
        K = find_K_from_k(k, M)[0]
        Kpts.append(K)
    parprint("running bands calculation")
    calc_bands = GPAW('licoo2_gs.gpw',
                      fixdensity=True,
                      kpts=Kpts,
                      symmetry='off',
                      nbands=-90,
                      convergence={'bands': -90})
    calc_bands.get_potential_energy()
    calc_bands.write('bands_3x3_defect.gpw', 'all')

    parprint("\n Unfolding......")
    parprint("\n length=" + str(len(kpts)))

    unfold = Unfold(name='3x3_defect',
                    calc='bands_3x3_defect.gpw',
                    M=M,
                    spinorbit=False)

    unfold.spectral_function(kpts=kpts, x=x, X=X, points_name=[])
    if world.rank == 0:
        os.rename("sf_3x3_defect.pckl", "sf_3x3_defect_{}.pckl".format(i))
        os.rename("weights_3x3_defect.pckl",
                  "weights_3x3_defect_{}.pckl".format(i))
    parprint("done ", i, "\n")
