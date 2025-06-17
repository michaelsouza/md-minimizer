We are implementing the experiments described in the `README.md` file.
We tried to do it using the `in.springs` file, but we failed.
Now, we are trying to convert the `in.springs` file to a Python version that calls LAMMPS.
Our LAMMPS install has the following capabilities:

```bash
(md-minimizer) (base) michael@ohm:~/gitrepos/md-minimizer$ lmp -h 

Large-scale Atomic/Molecular Massively Parallel Simulator - 12 Jun 2025
Git info (release / patch_12Jun2025)

Usage example: lmp -var t 300 -echo screen -in in.alloy

List of command-line options supported by this LAMMPS executable:

-echo none/screen/log/both  : echoing of input script (-e)
-help                       : print this help message (-h)
-in none/filename           : read input from file or stdin (default) (-i)
-kokkos on/off ...          : turn KOKKOS mode on or off (-k)
-log none/filename          : where to send log output (-l)
-mdi '<mdi flags>'          : pass flags to the MolSSI Driver Interface
-mpicolor color             : which exe in a multi-exe mpirun cmd (-m)
-cite                       : select citation reminder style (-c)
-nocite                     : disable citation reminder (-nc)
-nonbuf                     : disable screen/logfile buffering (-nb)
-package style ...          : invoke package command (-pk)
-partition size1 size2 ...  : assign partition sizes (-p)
-plog basename              : basename for partition logs (-pl)
-pscreen basename           : basename for partition screens (-ps)
-restart2data rfile dfile ... : convert restart to data file (-r2data)
-restart2dump rfile dgroup dstyle dfile ... 
                            : convert restart to dump file (-r2dump)
-restart2info rfile         : print info about restart rfile (-r2info)
-reorder topology-specs     : processor reordering (-r)
-screen none/filename       : where to send screen output (-sc)
-skiprun                    : skip loops in run and minimize (-sr)
-suffix gpu/intel/kk/opt/omp: style suffix to apply (-sf)
-var varname value          : set index style variable (-v)

OS: Linux "Ubuntu 24.04.2 LTS" 5.15.167.4-microsoft-standard-WSL2 x86_64

Compiler: GNU C++ 13.3.0 with OpenMP 4.5
C++ standard: C++17
Embedded fmt library version: 10.2.0
Embedded JSON class version: 3.12.0

MPI v3.1: Open MPI v5.0.8, package: Open MPI conda@74b3e81c684d Distribution, ident: 5.0.8, repo rev: v5.0.8, May 30, 2025

Accelerator configuration:


FFT information:

FFT precision  = double
FFT engine  = mpiFFT
FFT library = KISS

Active compile time flags:

-DLAMMPS_GZIP
-DLAMMPS_PNG
-DLAMMPS_JPEG
-DLAMMPS_FFMPEG
-DLAMMPS_SMALLBIG
sizeof(smallint): 32-bit
sizeof(imageint): 32-bit
sizeof(tagint):   32-bit
sizeof(bigint):   64-bit

Available compression formats:

Extension: .gz     Command: gzip
Extension: .bz2    Command: bzip2
Extension: .zst    Command: zstd
Extension: .xz     Command: xz
Extension: .lzma   Command: xz
Extension: .lz4    Command: lz4


Installed packages:

BPM MC MOLECULE PYTHON 

List of individual style options included in this LAMMPS executable

* Atom styles:

angle           atomic          body            bond            bpm/sphere      
charge          ellipsoid       full            hybrid          line            
molecular       sphere          template        tri             

* Integrate styles:

respa           verlet          

* Minimize styles:

cg              fire/old        fire            hftn            quickmin        
sd              

* Pair styles:

born            bpm/spring      buck            buck/coul/cut   coul/cut        
coul/debye      coul/dsf        coul/wolf       meam/c          reax            
reax/c          mesont/tpm      dsmc            hbond/dreiding/lj               
hbond/dreiding/morse            hybrid          hybrid/omp      hybrid/molecular                
hybrid/molecular/omp            hybrid/overlay  hybrid/overlay/omp              
hybrid/scaled   hybrid/scaled/omp               lj/charmm/coul/charmm           
lj/charmm/coul/charmm/implicit  lj/charmmfsw/coul/charmmfsh     lj/cut          
lj/cut/coul/cut lj/cut/tip4p/cut                lj/expand       morse           
python          soft            table           tip4p/cut       yukawa          
zbl             zero            

* Bond styles:

bpm/rotational  bpm/spring      bpm/spring/plastic              fene            
fene/expand     gromos          harmonic        hybrid          morse           
quartic         table           zero            

* Angle styles:

charmm          cosine          cosine/squared  harmonic        hybrid          
table           zero            

* Dihedral styles:

charmm          charmmfsw       harmonic        hybrid          multi/harmonic  
opls            table           zero            

* Improper styles:

cvff            harmonic        hybrid          umbrella        zero            

* KSpace styles:

zero            

* Fix styles

adapt           addforce        atom/swap       ave/atom        ave/chunk       
ave/correlate   ave/grid        ave/histo       ave/histo/weight                
ave/time        aveforce        balance         bond/break      bond/create     
bond/create/angle               bond/swap       box/relax       
charge/regulation               cmap            deform          deposit         
ave/spatial     ave/spatial/sphere              lb/pc           
lb/rigid/pc/sphere              reax/c/bonds    reax/c/species  dt/reset        
efield          enforce2d       evaporate       external        gcmc            
gravity         halt            heat            indent          langevin        
lineforce       mol/swap        momentum        move            nph             
nph/sphere      npt             npt/sphere      nve             nve/bpm/sphere  
nve/limit       nve/noforce     nve/sphere      nvt             nvt/sllod       
nvt/sphere      pair            planeforce      press/berendsen press/langevin  
print           property/atom   python/invoke   python          python/move     
recenter        restrain        set             setforce        spring          
spring/chunk    spring/self     store/force     store/state     temp/berendsen  
temp/rescale    tfmc            thermal/conductivity            vector          
viscous         wall/harmonic   wall/lj1043     wall/lj126      wall/lj93       
wall/morse      wall/reflect    wall/region     wall/table      widom           

* Compute styles:

aggregate/atom  angle           angle/local     angmom/chunk    bond            
bond/local      centro/atom     centroid/stress/atom            chunk/atom      
chunk/spread/atom               cluster/atom    cna/atom        com             
com/chunk       coord/atom      count/type      mesont          dihedral        
dihedral/local  dipole          dipole/chunk    displace/atom   erotate/sphere  
erotate/sphere/atom             fragment/atom   global/atom     group/group     
gyration        gyration/chunk  heat/flux       improper        improper/local  
inertia/chunk   ke              ke/atom         msd             msd/chunk       
nbond/atom      omega/chunk     orientorder/atom                pair            
pair/local      pe              pe/atom         pressure        property/atom   
property/chunk  property/grid   property/local  rdf             reduce          
reduce/chunk    reduce/region   slice           stress/atom     temp            
temp/chunk      temp/com        temp/deform     temp/partial    temp/profile    
temp/ramp       temp/region     temp/sphere     torque/chunk    vacf            
vacf/chunk      vcm/chunk       

* Region styles:

block           cone            cylinder        ellipsoid       intersect       
plane           prism           sphere          union           

* Dump styles:

atom            cfg             custom          atom/mpiio      cfg/mpiio       
custom/mpiio    xyz/mpiio       grid            grid/vtk        image           
local           movie           xyz             

* Command styles

angle_write     balance         change_box      create_atoms    create_bonds    
create_box      delete_atoms    delete_bonds    box             kim_init        
kim_interactions                kim_param       kim_property    kim_query       
reset_ids       reset_atom_ids  reset_mol_ids   message         server          
dihedral_write  displace_atoms  info            minimize        read_data       
read_dump       read_restart    replicate       rerun           run             
set             velocity        write_coeff     write_data      write_dump      
write_restart
```