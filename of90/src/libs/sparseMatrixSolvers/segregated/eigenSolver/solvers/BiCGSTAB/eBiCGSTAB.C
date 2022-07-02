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

#include "eBiCGSTAB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace EigenIterDirSolvers 
{
    defineTypeNameAndDebug(BiCGSTAB, 0);
    addToRunTimeSelectionTable(EigenIterDirSolver, BiCGSTAB, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EigenIterDirSolvers::BiCGSTAB::BiCGSTAB
(
    const dictionary& dict
)
:
EigenIterDirSolver(dict),
solverILU_(NULL),
solverDiag_(NULL),
precondType_(dict.subDict("preconditioner").lookup("preconditioner"))
{
    // Note: I know the pseudo-RTS for preconditioner is ugly, but since
    // there are only 2 'usable' preconditioners (actually only 1 is good - ILUT),
    // this highly simplifies the code, without significant overburden. 
    
    if (precondType_ == "ILUT")
    {
     solverILU_.reset
     (
        new ILUBiCGSTAB
     );
    }
    else if (precondType_ == "Diagonal")
    {
     solverDiag_.reset
     (
        new DiagBiCGSTAB
     );
    }
    else
    {
     FatalErrorIn
        (
            "Foam::EigenIterDirSolvers::BiCGSTAB::BiCGSTAB"
        )   << nl << "Unknown Eigen preconditioner type: " << precondType_ << nl 
            << "Available Eigen preconditioners for BiCGSTAB solver are: " << nl
            << " . ILUT" << nl << " . Diagonal "           
            << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EigenIterDirSolvers::BiCGSTAB::initialize  
(
   
)
{
   if (precondType_ == "ILUT")  
   {
     solverILU_->setMaxIterations(dict_.lookupOrDefault<int>("maxIter", 1000));
     solverILU_->setTolerance(dict_.lookupOrDefault<scalar>("tolerance", 1e-12)); 
  
     solverILU_->preconditioner().setDroptol(dict_.subDict("preconditioner").lookupOrDefault<scalar>("dropTol", 1e-12));  
     solverILU_->preconditioner().setFillfactor(dict_.subDict("preconditioner").lookupOrDefault<int>("fillFactor", 10));  
   }
   else
   {
     solverDiag_->setMaxIterations(dict_.lookupOrDefault<int>("maxIter", 1000));
     solverDiag_->setTolerance(dict_.lookupOrDefault<scalar>("tolerance", 1e-12));  
   }
}


void Foam::EigenIterDirSolvers::BiCGSTAB::setup  
(
  Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>& matrix,
  bool isSymmetric,
  bool reusePattern
)
{
  if (precondType_ == "ILUT")  
  {
    if (!reusePattern) 
      solverILU_->analyzePattern(matrix);   
    
    solverILU_->factorize(matrix);
    solverILU_->compute(matrix);
  }
  else
  {
    if (!reusePattern) 
      solverDiag_->analyzePattern(matrix);   
    
    solverDiag_->factorize(matrix);
    solverDiag_->compute(matrix);
  }
}

void Foam::EigenIterDirSolvers::BiCGSTAB::solve 
(
  const Eigen::VectorXd& source,
  Eigen::VectorXd& x
) 
{ 
  if (precondType_ == "ILUT")  
  {
    x = solverILU_->solve(source);
  
    residual_ = solverILU_->error();
    nIters_ = solverILU_->iterations(); 
  }
  else
  {
    x = solverDiag_->solve(source);
  
    residual_ = solverDiag_->error();
    nIters_ = solverDiag_->iterations(); 
  }
  
}

 
// ************************************************************************* //
