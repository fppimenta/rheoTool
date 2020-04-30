/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 
Class
    Foam::sparseSolverBase

Description
    This class contains static methods and data which are used for inter-processor
    communication, like the adressing at processor or cyclic patches. It also 
    handles counters for specific solver types (hypre and petsc) which need a 
    registry of all their instances for initialization or/and destruction (ex
    to know when to initialize/finalize petsc or/and MPI). Both coupled and 
    segregated solvers derived from here.
   
\*---------------------------------------------------------------------------*/

#ifndef sparseSolverBase_H
#define sparseSolverBase_H

#include "fieldTypes.H"
#include "List.H"
#include "fvMesh.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Struct sharedDataS Declaration
\*---------------------------------------------------------------------------*/

 // Struct used to indentify inter-processor connections 
 namespace spSolver
 {
  struct sharedDataS
  {
    List< List<label > > fCn;      // faceCells of neighbours
    List< List<label > > fCo;      // faceCells of owner
    List< List<label > > procInfo; // Info about connections
    int ilower; // Lowest cell ID (row) in a given processor
    int iupper; // Highest cell ID (row) in a given processor
    word meshName; // Name of mesh to which it belongs
  };
 }
 
/*---------------------------------------------------------------------------*\
                       Class sparseSolverBase Declaration
\*---------------------------------------------------------------------------*/

// This class is restricted to non-template (Type-independent) data 
class sparseSolverBase
{ 
  public:
  
  // Type-independent inside sparseSolver. Deriveds aware of each others
  // for convenience in managing MPI. 
  static int counterHypre_;
  static int counterPetsc_;
     
  // List of sharedDataS
  static List<spSolver::sharedDataS> sharedData;  
  
  // Fills sharedData. The mesh in each process gets numbered in
  // such a way that if the mesh was reconstructed (global space),
  // last cell of process n would precede the first cell of process
  // n+1. This is to ensure that the matrix rows are contiguous across
  // processors. In practice, we just need to shift the local ID of the 
  // cells by the amount of cells in the preceding process.          
  static void setSharedData
  (
    const fvMesh&
  );
        
  // Fills sharedData if it is empty
  static int buildSharedDataOnDemand
  (
    const fvMesh&
  );  
  
  // Debug switch
  static bool debug_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}  // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
