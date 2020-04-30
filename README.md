<img src="https://cloud.githubusercontent.com/assets/22405564/20934548/7163a14e-bbd3-11e6-84d3-e2e0ac073201.png" width="700">

## About

[RheoTool] is an open-source toolbox based on [OpenFOAM®] to simulate Generalized Newtonian Fluids (GNF) and viscoelastic fluids under pressure-driven and/or electrically-driven flows. In addition to continuum mechanics solvers, [RheoTool] also has a Brownian dynamics solver to simulate polymer molecules individually.  

The theory behind the single-phase flow solvers used in [RheoTool] can be found in [Pimenta F. and Alves M.A., 2017, J. Non-Newtonian Fluid Mech.](http://dx.doi.org/10.1016/j.jnnfm.2016.12.002) for pressure-driven flows, in [Pimenta F. and Alves M.A., 2018, arXiv:1802.02843](http://arxiv.org/abs/1802.02843) for electrically-driven flows and in [Pimenta F. and Alves M.A., 2020, arXiv:2004.14235](https://arxiv.org/abs/2004.14235) for non-isothermal flows with conjugate heat transfer. The Brownian dynamics solver was presented in [Pimenta F. et al., 2018, Biomicrofluidics](https://doi.org/10.1063/1.5037472).

Since version 4.0, [RheoTool] includes interfaces to the sparse matrix solvers of three external libraries: [Eigen](http://eigen.tuxfamily.org/), [Hypre](https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) and [PETSc](https://www.mcs.anl.gov/petsc/). Both coupled and segregated solution methods have been implemented, as described in [Pimenta F. and Alves M.A., 2019, Computer & Fluids](https://doi.org/10.1016/j.compfluid.2019.104279).

The library containing the viscoelastic fluid models has been developed based on the library used by _viscoelasticFluidFoam_ [(Favero et al., 2010, J. Non-Newtonian Fluid Mech.)](http://dx.doi.org/10.1016/j.jnnfm.2010.08.010), already present in [foam-extend 4.0]. 
 

## Features

* all the features are available for 2D/3D problems and generic grids;
* the code for continuum mechanics is fully-parallelized with MPI;
* all solvers support either moving or static meshes;
* non-isothermal flows, including conjugate heat transfer in multi-regions; 
* wide range of electrically-driven flow models;
* wide range of rheological models (viscoelastic, elastoviscoplastic and GNF);
* the log-conformation tensor approach is available for all viscoelastic models; 
* the transient flow solvers (_rheoFoam_ and _rheoEFoam_) are highly stable regarding pressure-stress-velocity coupling;
* the material functions of any rheological model can be obtained numerically (_rheoTestFoam_);
* wide range of (direct/iterative) sparse matrix solvers available through the interfaces to external libraries;
* segregated and coupled solvers;
* Brownian dynamics solver implementing bead-spring models; 
* a set of tutorials is included to illustrate the application of the solvers to different problems;
* the theory and the tutorials are described in a user-guide;
* a solver for two-phase flows is available (_rheoInterFoam_), where any GNF or viscoelastic model can be used for each phase (_under-development_);
* the tool is available for both [OpenFOAM®] and [foam-extend] versions.

## Installation

[RheoTool] can be either cloned using `git` via: `git clone https://github.com/fppimenta/rheoTool` or simply downloaded from the GitHub page at https://github.com/fppimenta/rheoTool.

The repository includes versions of [RheoTool] for: OpenFOAM® [v7.0](http://openfoam.org/version/7), OpenFOAM® [v6.0](http://openfoam.org/version/6) and foam-extend [v4.0](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0).

**Note: the [RheoTool] version for foam-extend is not updated since version 4.1.**

To install [RheoTool], please follow the instructions in Chapter 2 of the [user-guide](https://github.com/fppimenta/rheoTool/tree/master/doc).

## Docs

See the [user-guide](https://github.com/fppimenta/rheoTool/tree/master/doc).

## ChangeLog

Follow the version updates of rheoTool in the [ChangeLog](https://github.com/fppimenta/rheoTool/tree/master/ChangeLog). They can also be found in the [user-guide](https://github.com/fppimenta/rheoTool/tree/master/doc). 

## Third-Party 

[RheoTool] is using the following third-party packages:

* [Eigen](http://eigen.tuxfamily.org/).
* [Hypre](https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods).
* [PETSc](https://www.mcs.anl.gov/petsc/). 

## Screenshots

Here are some images from the tutorials included in [RheoTool]. 

<img src="https://user-images.githubusercontent.com/22405564/80705470-68921b80-8ade-11ea-96f6-3d3e7f4c68eb.png">

*This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks*

*OPENFOAM® is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com*

[RheoTool]:https://github.com/fppimenta/rheoTool
[OpenFOAM®]:http://openfoam.org/
[foam-extend]:https://github.com/Unofficial-Extend-Project-Mirror/Unofficial-Extend-Project-Mirror.github.io/wiki
[foam-extend 4.0]:https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0