from pymatgen.core import Structure as p
from pymatgen.io import pwscf
from matplotlib import pyplot as plt
import numpy as np
import subprocess
def make_input_r3m(species="K",vol=1.05):
    s=p.from_file("nacoo2.cif")
    s.replace_species({"Na":species})
    s.scale_lattice(s.volume*vol)
    species_pesudo={}
    species_pesudo["H"]="H.pbe-kjpaw_psl.1.0.0.UPF"
    species_pesudo["Li"]="Li.pbe-s-kjpaw_psl.1.0.0.UPF"
    species_pesudo["Na"]="Na.pbe-spn-kjpaw_psl.1.0.0.UPF"
    species_pesudo["K"]="K.pbe-spn-kjpaw_psl.1.0.0.UPF"
    species_pesudo["Cs"]="Cs.pbe-spn-kjpaw_psl.1.0.0.UPF"
    
    pseudo={f"{species}":species_pesudo[species],
            "Co":"Co.pbe-n-kjpaw_psl.1.0.0.UPF",
            "O":"O.pbe-n-kjpaw_psl.1.0.0.UPF"}
    control={'pseudo_dir' : '/home/srr70/QE/pseudo/pslibrary.1.0.0/pbe/PSEUDOPOTENTIALS/',
            'verbosity' :'high',
            'prefix':'lco',
            'calculation': 'vc-relax',
            'restart_mode':'from_scratch',
            'nstep':100,
            'wf_collect':True,
            'outdir': './',
            'tprnfor':True,
            'tstress':True}
    system={"ecutwfc":80.0 ,
            "ecutrho" :240.0,
            "occupations":'smearing',
            "smearing":'mp',
            "degauss":0.02,}
    electrons={"diagonalization":'david'}
    pw=pwscf.PWInput(s,pseudo=pseudo,
                    control=control,
                    kpoints_grid=(8,8,2),
                    electrons=electrons,
                    system=system)
    pw.write_file("scf.in")
    

def make_input_tetra(species="K",vol=1.05):
    s=p.from_file("tetra.cif")
    s.replace_species({"K":species})
    s.scale_lattice(s.volume*vol)
    species_pesudo={}
    species_pesudo["H"]="H.pbe-kjpaw_psl.1.0.0.UPF"
    species_pesudo["Li"]="Li.pbe-s-kjpaw_psl.1.0.0.UPF"
    species_pesudo["Na"]="Na.pbe-spn-kjpaw_psl.1.0.0.UPF"
    species_pesudo["K"]="K.pbe-spn-kjpaw_psl.1.0.0.UPF"
    species_pesudo["Cs"]="Cs.pbe-spn-kjpaw_psl.1.0.0.UPF"
    
    pseudo={f"{species}":species_pesudo[species],
            "Co":"Co.pbe-n-kjpaw_psl.1.0.0.UPF",
            "O":"O.pbe-n-kjpaw_psl.1.0.0.UPF"}
    control={'pseudo_dir' : '/home/srr70/QE/pseudo/pslibrary.1.0.0/pbe/PSEUDOPOTENTIALS/',
            'verbosity' :'high',
            'prefix':'lco',
            'calculation': 'vc-relax',
            'restart_mode':'from_scratch',
            'nstep':100,
            'wf_collect':True,
            'outdir': './',
            'tprnfor':True,
            'tstress':True}
    system={"ecutwfc":80.0 ,
            "ecutrho" :240.0,
            "occupations":'smearing',
            "smearing":'mp',
            "degauss":0.02,}
    electrons={"diagonalization":'david'}
    pw=pwscf.PWInput(s,pseudo=pseudo,
                    control=control,
                    kpoints_grid=(8,8,3),
                    electrons=electrons,
                    system=system)
    pw.write_file("scf.in")



def calculate(species):
    cmd = f"mkdir {species}"
    subprocess.call(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    cmd=f"cp scf.in {species}"
    subprocess.call(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    cmd = "mpirun -np 32 /home/srr70/install/qe_6.6/q-e-qe-6.6/bin/pw.x -in scf.in > scf.out"
    subprocess.call(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,cwd=f"{species}")


# for r3m
# system=[["H",.8],["Li",.93],["Na",1],["K",1,1],["Cs",1.5]] 

#for tetra
system=[["H",.8],["Li",.88],["Na",.93],["K",1],["Cs",1.1]]

for i in system:
    species=i[0]
    vol=i[1]
    print(f"calculating for {species}")
    #make_input_r3m(species,vol)
    make_input_tetra(species,vol)
    calculate(species)