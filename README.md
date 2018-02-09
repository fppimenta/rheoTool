<img src="https://cloud.githubusercontent.com/assets/22405564/20934548/7163a14e-bbd3-11e6-84d3-e2e0ac073201.png" width="700">

## About

[RheoTool] is an open-source toolbox based on [OpenFOAM®] to simulate Generalized Newtonian Fluids (GNF) and viscoelastic fluids under pressure-driven and/or electrically-driven flows. 

The theory behind the single-phase flow solvers used in [RheoTool] can be found in [Pimenta F. and Alves M.A., 2017, J. Non-Newtonian Fluid Mech.](http://dx.doi.org/10.1016/j.jnnfm.2016.12.002) for pressure-driven flows, and in [Pimenta F. and Alves M.A., 2018, arXiv:1802.02843](http://arxiv.org/abs/1802.02843) for electrically-driven flows.

The library containing the viscoelastic models has been developed based on the library used by _viscoelasticFluidFoam_ [(Favero et al., 2010, J. Non-Newtonian Fluid Mech.)](http://dx.doi.org/10.1016/j.jnnfm.2010.08.010), already present in [foam-extend 4.0]. 
 

## Features

* all the features are available for 2D/3D problems and generic grids;
* the code is fully-parallelized;
* wide range of electrically-driven flow models;
* wide range of viscoelastic and GNF models;
* the log-conformation tensor approach is available for all viscoelastic models; 
* the transient flow solvers (_rheoFoam_ and _rheoEFoam_) are highly stable regarding pressure-stress-velocity coupling;
* the material functions of any rheological model can be obtained numerically (_rheoTestFoam_);
* a set of tutorials is included to illustrate the application of the solvers to different problems;
* the theory and the tutorials are described in a user-guide;
* a solver for two-phase flows is available (_rheoInterFoam_), where any GNF or viscoelastic model can be used for each phase (_under-development_);
* the tool is available for both [OpenFOAM®] and [foam-extend] versions.

## Installation

[RheoTool] can be either cloned using `git` via: `git clone https://github.com/fppimenta/rheoTool` or simply downloaded from the GitHub page at https://github.com/fppimenta/rheoTool.

The repository includes versions of [RheoTool] for: OpenFOAM® [v2.2.2](http://openfoam.org/version/2-2-2), OpenFOAM® [v4.1](http://openfoam.org/version/4-1)/[v4.0](http://openfoam.org/version/4-0) and foam-extend [v4.0](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0). 

To install [RheoTool], please follow the instructions in Chapter 2 of the [user-guide](https://github.com/fppimenta/rheoTool/tree/master/doc).

## Docs

See the [user-guide](https://github.com/fppimenta/rheoTool/tree/master/doc).

## Third-Party 

[RheoTool] is using the following third-party packages:

* [Eigen](http://eigen.tuxfamily.org/). 

## Screenshots

Here are some images from the tutorials included in [RheoTool]. 

<img src="https://user-images.githubusercontent.com/22405564/36007099-b988edaa-0d38-11e8-8f9c-157e78f040a6.png">

[RheoTool]:https://github.com/fppimenta/rheoTool
[OpenFOAM®]:http://openfoam.org/
[foam-extend]:http://www.extend-project.de/
[foam-extend 4.0]:https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0