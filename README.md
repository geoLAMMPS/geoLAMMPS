# geoLAMMPS
This repository contains a fork of [LAMMPS](https://www.lammps.org/) which includes various additions, summarised below, to facilitate the simulation of soil mechanics problems with DEM. Development of this repository started in 2011 in the group of [Prof. Catherine O'Sullivan](https://profiles.imperial.ac.uk/cath.osullivan/) at Imperial College London. The main contributors to this repository since 2011 have been:
- [Dr Kevin Hanley](https://www.eng.ed.ac.uk/about/people/dr-kevin-j-hanley), now at the University of Edinburgh
- [Dr George Marketos](https://orcid.org/0000-0002-6441-0434), now at COWI UK, London
- [Dr Masahide Otsubo](https://orcid.org/0000-0001-6738-2160), now at Port and Airport Research Institute (PARI), Japan

Until 2024, this repository was private on Bitbucket, being accessible only by invitation from the repository administrators. 

## How to cite
If you publish anything which has used this geoLAMMPS code, we would appreciate it if you mention the code name and cite the most relevant paper(s) corresponding to whichever code features you have used that are not part of the main LAMMPS distribution. Most of the additions below have an associated paper which either describes the code addition or relies upon it for the research.

## Main additions compared to the main LAMMPS distribution

#### Advanced particle–particle and particle–wall contact models
These models are mostly not documented in the usual LAMMPS manner, but can be understood [from the examples provided](examples/geoLAMMPS) and the source code.

- *PairGranShmHistory* is our implementation of a 'simplified Hertz–Mindlin' contact model
- *PairGranHookeHistoryOldstyle* is our linear, Hookean contact model implementation
- *PairGranCMHistory* is the Cavarretta–Mindlin contact model for rough particles
  - The reference paper for *PairGranCMHistory* is Otsubo, M., O'Sullivan, C., Hanley, K. J. & Sim, W. W. (2017). [The influence of particle surface roughness on elastic stiffness and dynamic response](https://doi.org/10.1680/jgeot.16.P.050), Géotechnique, 67(5), 452–459.
- *PairGranCMDHistory* is the Cavarretta–Mindlin–Deresiewicz contact model for rough particles, accounting for partial slip effects
  - The same reference paper applies for the normal component of *PairGranCMDHistory* as for *PairGranCMHistory*.
- *PairGranHMDHistory* is the Hertz–Mindlin–Deresiewicz contact model, including partial slip
  - The reference publication for *PairGranHMDHistory* is Chapter 4 of the PhD thesis of Masahide Otsubo: [Particle scale analysis of soil stiffness and elastic wave propagation](https://doi.org/10.25560/44380), 2017, Imperial College London.

All of these contact models have the following common features:
- energy tracing (details below)
- the ability to enable advanced rolling and twisting/spin resistance models (details below)
- an [expanded range of per-contact ComputePairLocal entries](doc/USER/gran/Interpreting_ComputePairLocal_Output.doc)
- equivalent particle–wall contact models through *FixWallGranOldstyle* and *FixWallGranRegionOldstyle*
- access to expanded wall motion capabilities via added keywords in *FixWallGranOldstyle* and *FixWallGranRegionOldstyle*: *translate* to specify a wall velocity and *stresscontrol* to achieve a desired force/stress on a planar wall

#### Rolling and twisting/spin resistance models
Two rolling and twisting/spin resistance models are available:
1. An advanced rotational (rolling & twisting) resistance model with a physical basis derived by Dr Xin Huang (now at Tongji University, Shanghai)
   - This is [enabled with the *rolling* keyword and associated arguments in pair\_modify](doc/src/pair_modify.rst).
   - This rotational resistance model is only available when *PairGranShmHistory* or *PairGranHookeHistoryOldstyle* is selected.
   - The reference paper for this model is Huang, X., Hanley, K.J., O'Sullivan, C. & Kwok, C.-Y. (2017). [Implementation of rotational resistance models: a critical appraisal](https://doi.org/10.1016/j.partic.2016.08.007), Particuology, 34, 14–23.
2. A spin resistance model theoretically described by [Deresiewicz (1954)](https://doi.org/10.1115/1.4010818)
   - Accessible through the *D\_spin* keyword in pair\_modify (not documented)
   - This spin resistance model is only available when *PairGranShmHistory*, *PairGranCMHistory*, *PairGranCMDHistory* or *PairGranHMDHistory* is selected.
   - The reference publication is Chapter 4 of the PhD thesis of Masahide Otsubo: [Particle scale analysis of soil stiffness and elastic wave propagation](https://doi.org/10.25560/44380), 2017, Imperial College London.

#### Energy tracing
Using *ComputeEnergyGran*, various specified components of energy (kinetic energy, boundary work, strain energy at contacts, energy dissipated by damping and friction) may be computed and exported periodically for any simulation using a suitable granular pairstyle, including those five added contact models listed above.
- Per-contact energy terms can be exported using *DumpLocal* provided that the *trace\_energy* option of pair\_modify is enabled in the script.
- This is [fully documented](doc/src/compute_energy_gran.rst) and an [example has been provided](examples/geoLAMMPS/Energy_Tracing) to demonstrate the energy balance.
- The reference paper is Hanley, K.J., Huang, X. & O'Sullivan, C. (2018). [Energy dissipation in soil samples during drained triaxial shearing](https://doi.org/10.1680/jgeot.16.P.317), Géotechnique, 68(5), 421–433.

#### Stress/strain control of periodically-bounded samples
*FixMultistress* allows any reasonable permutation of strain and stress control to be applied to a sample bounded by periodic boundaries, whether in an orthogonal or a triclinic bounding box. It allows for the simulation of undrained triaxial tests by maintaining a constant volume, in addition to drained compression/extension triaxial tests, tests at constant mean effective stress (p'), constant deviator stress (q), cyclic loading, true triaxial tests at constant intermediate stress ratios, b, between 0 and 1, and more besides.
- This is [fully documented](doc/src/fix_multistress.rst) and [examples have been provided](examples/geoLAMMPS) to demonstrate this.
- Much of the functionality of *FixMultistress* [was incorporated into the main LAMMPS distribution](https://github.com/lammps/lammps/pull/4017) in 2024 as an extension of *FixDeform*. The prospective user should consider whether it would be more convenient to use that instead.
- For drained triaxial (constant minor principal effective stress), constant volume or constant p' stress paths, there are two relevant reference publications to choose from, both of which use all three of these conditions:
  - Huang, X., O'Sullivan, C., Hanley, K.J. & Kwok, C.-Y. (2014). [Discrete-element method analysis of the state parameter](https://doi.org/10.1680/geot.14.P.013), Géotechnique, 64(12), 954–965.
  - Huang, X., Hanley, K.J., O'Sullivan, C. & Kwok, C.-Y. (2014). [Exploring the influence of interparticle friction on critical state behaviour using DEM](https://doi.org/10.1002/nag.2259), International Journal for Numerical and Analytical Methods in Geomechanics, 38(12), 1276–1297.
- For simulations using *FixMultistress* to conduct true triaxial, constant intermediate stress ratio (b) simulations, the reference publication is Huang, X., Hanley, K.J., O'Sullivan, C., Kwok, C.-Y. & Wadee, M.A. (2014). [DEM analysis of the influence of the intermediate stress ratio on the critical-state behaviour of granular materials](https://doi.org/10.1007/s10035-014-0520-6), Granular Matter, 16(5), 641–655.
- A more comprehensive description of the equations being implemented for these conditions is provided in Chapter 2 of the PhD thesis of Xin Huang: [Exploring critical-state behaviour using DEM](https://doi.org/10.25560/25316), 2014, Imperial College London.
- For cyclic loading simulations using *FixMultistress*, the reference paper is Keishing, J., Huang, X. & Hanley, K.J. (2020). [Energy dissipation in soil samples during cyclic triaxial simulations](https://doi.org/10.1016/j.compgeo.2020.103481), Computers and Geotechnics, 121, 103481.

#### Particle crushing
*FixCrushing* allows particle crushing during a simulation. This is effected by reducing the diameter of a crushed particle once a failure criterion has been met, and inserting small particles into voids to conserve the total sample mass. This approach is computationally tractable and the implementation allows some randomness in the crushing responses.
- The version of *FixCrushing* in geoLAMMPS uses the two-parameter failure criterion for brittle materials proposed by [Christensen (2000)](https://doi.org/10.1177/108128650000500302).
- This is [fully documented](doc/src/fix_crushing.rst) and an [example has been provided](examples/geoLAMMPS/FixCrushing) to demonstrate this.
- The reference paper is Hanley, K.J., O'Sullivan, C. & Huang, X. (2015). [Particle-scale mechanics of sand crushing in compression and shearing using DEM](https://doi.org/10.1016/j.sandf.2015.09.011), Soils and Foundations, 55(5), 1100–1112.
- It is noted that a different failure criterion (based on the maximum force) was used for this paper, as described in Sec. 2.2, compared to the Christensen failure criterion provided in geoLAMMPS.

#### Other additions
A number of other smaller capabilities have been added to geoLAMMPS, most of which are documented and demonstrated by the [provided examples](examples/geoLAMMPS):
- Local damping [is provided by *FixDampLocal*](doc/src/fix_damp_local.rst). This is applied at the particle-scale and is convenient for removing energy from a system during, e.g., preparation of an isotropic sample for triaxial shearing.
- *FixMomentumGran* is used to [zero the angular momentum of each particle in a group on every timestep](doc/src/fix_momentum_gran.rst). This was originally developed to facilitate validation simulations using the approach of [Thornton (1979)](https://doi.org/10.1680/geot.1979.29.4.441), based on face-centred cubic arrays. An example of this [analytical validation has been provided](examples/geoLAMMPS/Thornton_FCC) with the source code.
- *timestep auto* is used to [compute a stable simulation timestep automatically](doc/src/timestep.rst) for the *sphere* atom\_style, SI units and a granular pairstyle.
- *FixFluidDrag* [provides a simple drag force model](doc/src/fix_fluiddrag.rst) for simulating the flow of water through a packing of particles. The *sphere* atom\_style and SI units are required.

## Licence
LAMMPS and this geoLAMMPS fork are released under a GNU Public License Version 2 (GPLv2).

> The code comes with no warranty of any kind. The contributors to this repository shall not be responsible for the use of your code or any information contained in the code, its documentation, or any other source referring to the code or its documentation.

## Compiling geoLAMMPS

geoLAMMPS can be built using either GNU Make or CMake in a similar manner to the standard LAMMPS distribution. The following are minimal examples to configure and compile geoLAMMPS in parallel. Many of the additions described above are part of the GRANULAR package.

A. Using GNU Make

```
cd <path-to-lammps>/src
make yes-GRANULAR 	# install the GRANULAR package in the source tree
make mpi          	# build the LAMMPS executable with MPI
```

B. Using CMake

```
cd <path-to-lammps>                         # change to the geoLAMMPS distribution directory
mkdir build; cd build                       # create and change to build directory
cmake -D BUILD_MPI=yes PKG_GRANULAR=yes ../cmake/ # include the GRANULAR package
cmake --build .                             # compilation (or type "make")
```