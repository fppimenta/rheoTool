/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "eSparseLU.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace EigenIterDirSolvers 
{
    defineTypeNameAndDebug(SparseLU, 0);
    addToRunTimeSelectionTable(EigenIterDirSolver, SparseLU, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::EigenIterDirSolvers::SparseLU::SparseLU
(
    const dictionary& dict
)
:
EigenIterDirSolver(dict),
sparseLU_(NULL)
{
    sparseLU_.reset
    (
       new sparseLU
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::EigenIterDirSolvers::SparseLU::initialize  
(
   
)
{ 
   sparseLU_->setPivotThreshold(dict_.lookupOrDefault<scalar>("pivotThreshold", 1e-12));    
}


void Foam::EigenIterDirSolvers::SparseLU::setup  
(
  Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>& matrix,
  bool isSymmetric,
  bool reusePattern
)
{ 
 if (!reusePattern) 
   sparseLU_->analyzePattern(matrix);   
    
 sparseLU_->factorize(matrix);
 sparseLU_->compute(matrix);
}

void Foam::EigenIterDirSolvers::SparseLU::solve 
(
  const Eigen::VectorXd& source,
  Eigen::VectorXd& x
) 
{ 
  x = sparseLU_->solve(source);
  
  residual_ = 1e-16; // Should be O() machine precision upon convergence
  nIters_ = 1; // single iter for direct solvers  
}

// ************************************************************************* //
