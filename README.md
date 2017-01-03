<img src="https://cloud.githubusercontent.com/assets/22405564/20934548/7163a14e-bbd3-11e6-84d3-e2e0ac073201.png" width="700">
 
## About

[RheoTool] is an open-source toolbox based on [OpenFOAM®] to simulate the flow of Generalized Newtonian Fluids (GNF) and viscoelastic fluids. The library containing the viscoelastic models has been developed based on the library used by _viscoelasticFluidFoam_ [(Favero et al., 2010, J. Non-Newtonian Fluid Mech.)](http://dx.doi.org/10.1016/j.jnnfm.2010.08.010), already present in [foam-extend 3.2]. 

The theory behind the single-phase flow solver used in [RheoTool] can be found in: [Pimenta F. and Alves M.A., 2017, J. Non-Newtonian Fluid Mech.](http://dx.doi.org/10.1016/j.jnnfm.2016.12.002).

## Features

* all the features are available for 2D/3D problems and generic grids;
* the code is fully-parallelized;
* wide range of viscoelastic and GNF models;
* the log-conformation tensor approach is available for all viscoelastic models; 
* the transient flow solver (_rheoFoam_) is highly stable regarding pressure-stress-velocity coupling;
* the material functions of any model can be obtained numerically (_rheoTestFoam_);
* a group of tutorials is included to exemplify the application of the solvers to different problems;
* the theory and the tutorials are described in a user-guide;
* a solver for two-phase flows is available (_rheoInterFoam_), where any GNF or viscoelastic model can be used for each phase (_in-development_);
* the tool is available for both [OpenFOAM®] and [foam-extend] versions.

## Installation

[RheoTool] can be either cloned using `git` via: `git clone https://github.com/fppimenta/rheoTool` or simply downloaded from the GitHub page at https://github.com/fppimenta/rheoTool.  

The repository includes versions of [RheoTool] for: OpenFOAM® [v2.2.2](http://openfoam.org/version/2-2-2), OpenFOAM® [v4.0](http://openfoam.org/version/4-0) and foam-extend [v3.2](https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-3.2/tree/nextRelease) (branch _nextRelease_; lastly tested on commit 6de4e26).  

To install [RheoTool], please follow the instructions in Chapter 2 of the [user-guide](https://github.com/fppimenta/rheoTool/tree/master/doc).

## Docs

See the [user-guide](https://github.com/fppimenta/rheoTool/tree/master/doc).

## Third-Party 

[RheoTool] is using the following third-party packages:

* [Eigen](http://eigen.tuxfamily.org/). 

## Screenshots

Here are some pictures from the tutorials included in [RheoTool]. 

<img src="https://cloud.githubusercontent.com/assets/22405564/20934556/778ed66a-bbd3-11e6-99a7-7d1ff75757fd.png">

[RheoTool]:https://github.com/fppimenta/rheoTool
[OpenFOAM®]:http://openfoam.org/
[foam-extend]:http://www.extend-project.de/
[foam-extend 3.2]:https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-3.2/tree/nextRelease

