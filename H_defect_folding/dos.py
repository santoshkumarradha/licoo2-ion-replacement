from gpaw import GPAW, restart
import pickle
from ase.parallel import parprint, world
_, calc = restart('bands_3x3_defect.gpw')

parprint("Running DOS calculation .....")
e, dos = calc.get_dos(spin=0, npts=10001, width=5e-3)
e_f = calc.get_fermi_level()

parprint("\n Done, Saving")
if world.rank == 0:
    with open('dos.pickle', 'wb') as handle:
        pickle.dump({
            "ef": e_f,
            "e": e,
            "dos": dos
        },
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL)