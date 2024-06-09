# Details

## Servers

* <http://swissparam.ch>
* <https://cgenff.silcsbio.com>

## Commands

```shell
perl sort_mol2_bonds.pl ./LIG.mol2 ./LIG.mol2

# or

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/autodocking:2024.2.11 sort_bonds ./LIG.mol2 ./LIG.mol2
```

Generate stream file from <https://cgenff.silcsbio.com> server.

```shell
python ./cgenff_charmm2gmx.py LIG ./LIG.mol2 ./LIG.str ./charmm36-jul2022.ff

# or

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/autodocking:2024.2.11 cgenff LIG ./LIG.mol2 ./LIG.str ./charmm36-jul2022.ff
```

Copy 1t4w.pdb, ligand files and charmm36-jul2022.ff force field in the folder.

In linux we can compact the command

```shell
export GMX="docker run --gpus all -it --rm -v ${PWD}:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx"
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx pdb2gmx -f 2gs2_repaired.pdb -o rec.gro -ter -ignh

# or on linux

$GMX pdb2gmx -f 2gs2_repaired.pdb -o rec.gro -ter -ignh
```

Select charmm36-jul2022 force field and its recommended water model.

Now generate ligand gro file

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx editconf -f lig_ini.pdb -o ligand.gro

# or on linux

$GMX editconf -f lig_ini.pdb -o ligand.gro
```

Make complex.gro file combining rec.gro and ligand.gro. Calculate the atom number in complex.gro.

Update topol.top file to add ligand parameters and topology. And then ligand in molecules.

```text
; Include ligand parameters
#include "lig.prm"
```

```text
; Include ligand topology
#include "lig.itp"
```

```text
LIG                 1
```

Create box for simulation.

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0

# or on linux

$GMX editconf -f complex.gro -o newbox.gro -bt dodecahedron -d 1.0
```

Solvate the system.

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

# or on linux

$GMX solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

# or on linux

$GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral

# or on linux

$GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname SOD -nname CLA -neutral
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

# or on linux

$GMX grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx mdrun -v -deffnm em -nb gpu

# or on linux

$GMX mdrun -v -deffnm em -nb gpu
```

Create index

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx make_ndx -f ligand.gro -o index_lig.ndx

# or on linux

$GMX make_ndx -f ligand.gro -o index_lig.ndx
```

```text
- Selection: 0 & ! a H*
- Then for quit: q

Here,
    0 means selecting 0 (System) from list
    ! a H* means excluding all H atoms

** Possible selection:
nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx genrestr -f ligand.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000

# or on linux

$GMX genrestr -f ligand.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000
```

```text
- Selection: Group     3 (   System_&_!H*) has    30 elements
```

Add posre_lig.itp file in topology

```text
#ifdef POSRES
#include "posre_lig.itp"
#endif
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx make_ndx -f em.gro -o index.ndx

# or on linux

$GMX make_ndx -f em.gro -o index.ndx
```

```text
- Selection: 1 | 13
- Selection: 15 | 14
- Then for quit: q

Here,
    1 Protein             :  n atoms
    13 LIG                :  n atoms
    14 CLA/SOD            :  n atoms
    15 Water              :  n atoms
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx mdrun -v -deffnm nvt -nb gpu

# or on linux

$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
$GMX mdrun -v -deffnm nvt -nb gpu
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx mdrun -v -deffnm npt -nb gpu

# or on linux

$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
$GMX mdrun -v -deffnm npt -nb gpu
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx mdrun -v -deffnm md -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu

# or on linux

$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
$GMX mdrun -v -deffnm md -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx mdrun -v -deffnm ie -rerun md.xtc -nb cpu

# or on linux

$GMX grompp -f ie.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o ie.tpr
$GMX mdrun -v -deffnm ie -rerun md.xtc -nb cpu
```

### Analysis

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx energy -f md.edr -o energy.xvg

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx energy -f nvt.edr -o temperature.xvg

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx energy -f npt.edr -o pressure.xvg

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx energy -f npt.edr -o density.xvg

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx energy -f ie.edr -o interaction_energy.xvg

# or on linux

$GMX energy -f md.edr -o energy.xvg
$GMX energy -f nvt.edr -o temperature.xvg
$GMX energy -f npt.edr -o pressure.xvg
$GMX energy -f npt.edr -o density.xvg
$GMX energy -f ie.edr -o interaction_energy.xvg
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx trjconv -s nvt.tpr -f nvt.xtc -o nvt_center.xtc -pbc mol -center -ur compact

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx trjconv -s npt.tpr -f npt.xtc -o npt_center.xtc -pbc mol -center -ur compact

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -pbc mol -center -ur compact

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx trjconv -s md.tpr -f md_center.xtc -o start.pdb -dump 0

# or on linux

$GMX trjconv -s nvt.tpr -f nvt.xtc -o nvt_center.xtc -pbc mol -center -ur compact

$GMX trjconv -s npt.tpr -f npt.xtc -o npt_center.xtc -pbc mol -center -ur compact

$GMX trjconv -s md.tpr -f md.xtc -o md_center.xtc -pbc mol -center -ur compact

$GMX trjconv -s md.tpr -f md_center.xtc -o start.pdb -dump 0
```

```text
- Selection: 1
- Selection: 0
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx hbond -s md.tpr -f md_center.xtc -num hb.xvg

# or on linux

$GMX hbond -s md.tpr -f md_center.xtc -num hb.xvg
$GMX hbond -s md.tpr -f md_center.xtc -n index.ndx -num hb.xvg
$GMX hbond -s md.tpr -f md_center.xtc -n index.ndx -num -dist hb.xvg
$GMX hbond -s md.tpr -f md_center.xtc -n index.ndx -ac hb_ang.xvg
```

```text
- Selection: 1
- Selection: 13
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx gyrate -s md.tpr -f md_center.xtc -o gyrate.xvg

# or on linux

$GMX gyrate -s md.tpr -f md_center.xtc -o gyrate.xvg
```

```text
- Selection: 1
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx sasa -s md.tpr -f md_center.xtc -n index.ndx -o sasa.xvg -odg odg.xvg -or or.xvg -tv tv.xgv -b 0 -tu ns
#  Group 18 "Protein_LIG" (7212 atoms)

# or on linux

$GMX sasa -s md.tpr -f md_center.xtc -n index.ndx -o sasa.xvg -odg odg.xvg -or or.xvg -tv tv.xgv -b 0 -tu ns

docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx sasa -s md.tpr -f md_center.xtc -n index.ndx -o sasa_p.xvg -odg odg_p.xvg -or or_p.xvg -tv tv_p.xgv -b 0 -tu ns
#  Group  1 "Protein" (7155 atoms)

# or on linux

$GMX sasa -s md.tpr -f md_center.xtc -n index.ndx -o sasa_p.xvg -odg odg_p.xvg -or or_p.xvg -tv tv_p.xgv -b 0 -tu ns
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx rms -s md.tpr -f md_center.xtc -o rmsd.xvg

# or on linux

$GMX rms -s md.tpr -f md_center.xtc -o rmsd.xvg
```

```text
- Selection: Group     4 (       Backbone)
- Selection: Group     4 (       Backbone)
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg

# or on linux

$GMX rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg
```

```text
- Selection: Group     4 (       Backbone)
```

B-factor analysis with RMSF

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg -oq bfac.pdb

# or on linux

$GMX rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg -oq bfac.pdb
```

```text
- Selection: Group     3 (        C-alpha)
```

RMSD of Complex

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx rms -s md.tpr -f md_center.xtc -o rmsd_protein.xvg

# or on linux

$GMX rms -s md.tpr -f md_center.xtc -o rmsd_protein.xvg
```

```text
- Selection: Group     4 (       Backbone)
- Selection: Group     4 (       Backbone)
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx rms -s md.tpr -f md_center.xtc -o rmsd_ligand.xvg

# or on linux

$GMX rms -s md.tpr -f md_center.xtc -o rmsd_ligand.xvg
```

```text
- Selection: Group    13 (            LIG)
- Selection: Group    13 (            LIG)
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx rms -s md.tpr -f md_center.xtc -n index.ndx -o rmsd_complex.xvg

# or on linux

$GMX rms -s md.tpr -f md_center.xtc -n index.ndx -o rmsd_complex.xvg
```

```text
- Selection: Group     1 (        Protein)
- Selection: Group    13 (            LIG)
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx covar -s md.tpr -f md_center.xtc -n index.ndx -o eigenval.xvg

$GMX covar -s md.tpr -f md_center.xtc -n index.ndx -o eigenval.xvg

- Selection: Group     3 (        C-alpha)
- Selection: Group    13 (            LIG)
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx anaeig -v eigenvec.trr -s md.tpr -f md_center.xtc -n index.ndx -comp eigcomp.xvg -rmsf eigrmsf.xvg -2d eig2dproj.xvg -first 1 -last 5

$GMX anaeig -v eigenvec.trr -s md.tpr -f md_center.xtc -n index.ndx -comp eigcomp.xvg -rmsf eigrmsf.xvg -2d eig2dproj.xvg -first 1 -last 5

- Selection: Group     3 (        C-alpha)
- Selection: Group    13 (            LIG)
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx anaeig -v eigenvec.trr -s md.tpr -f md_center.xtc -n index.ndx -filt filter1.pdb -first 1 -last 1 -skip 100

$GMX anaeig -v eigenvec.trr -s md.tpr -f md_center.xtc -n index.ndx -filt filter1.pdb -first 1 -last 1 -skip 100

- Selection: Group     3 (        C-alpha)
- Selection: Group    13 (            LIG)
```

```shell
docker run --gpus all -it --rm -v .:/host_pwd --workdir /host_pwd firesimulations/gromacs:2023.11.26.rdtscp gmx anaeig -v eigenvec.trr -s md.tpr -f md_center.xtc -n index.ndx -extr extreme1.pdb -first 1 -last 1 -nframes 30

$GMX anaeig -v eigenvec.trr -s md.tpr -f md_center.xtc -n index.ndx -extr extreme1.pdb -first 1 -last 1 -nframes 30

- Selection: Group     3 (        C-alpha)
- Selection: Group    13 (            LIG)
```

## Commands for Git Large File Storage

```shell
git lfs install
git lfs track "*.xtc"
git lfs push --all origin main
git lfs uninstall
```
