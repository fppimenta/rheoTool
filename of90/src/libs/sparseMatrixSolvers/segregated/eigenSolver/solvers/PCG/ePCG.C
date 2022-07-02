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

#include "ePCG.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace EigenIterDirSolvers 
{
    defineTypeNameAndDebug(ConjugateGradient, 0);
    addToRunTimeSelectionTable(EigenIterDirSolver, ConjugateGradient, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EigenIterDirSolvers::ConjugateGradient::ConjugateGradient
(
    const dictionary& dict
)
:
EigenIterDirSolver(dict),
solverICC_(NULL),
solverDiag_(NULL),
precondType_(dict.subDict("preconditioner").lookup("preconditioner"))
{
    // Note: I know the pseudo-RTS for preconditioner is ugly, but since
    // there are only 2 'usable' preconditioners (actually only 1 is good - ICC),
    // this highly simplifies the code, without significant overburden. 
    
    if (precondType_ == "ICC")
    {
     solverICC_.reset
     (
        new ICCPCG
     );
    }
    else if (precondType_ == "Diagonal")
    {
     solverDiag_.reset
     (
        new DiagPCG
     );
    }
    else
    {
     FatalErrorIn
        (
            "Foam::EigenIterDirSolvers::ConjugateGradient::ConjugateGradient"
        )   << nl << "Unknown Eigen preconditioner type: " << precondType_ << nl 
            << "Available Eigen preconditioners for PCG solver are: " << nl
            << " . ICC" << nl << " . Diagonal "           
            << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::EigenIterDirSolvers::ConjugateGradient::initialize  
(
   
)
{
   if (precondType_ == "ICC")  
   {
     solverICC_->setMaxIterations(dict_.lookupOrDefault<int>("maxIter", 1000));
     solverICC_->setTolerance(dict_.lookupOrDefault<scalar>("tolerance", 1e-12)); 
   }
   else
   {
     solverDiag_->setMaxIterations(dict_.lookupOrDefault<int>("maxIter", 1000));
     solverDiag_->setTolerance(dict_.lookupOrDefault<scalar>("tolerance", 1e-12));  
   }
}


void Foam::EigenIterDirSolvers::ConjugateGradient::setup  
(
  Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>& matrix,
  bool isSymmetric,
  bool reusePattern
)
{
  if (!isSymmetric)
  {
    FatalErrorIn
    (
      "Foam::EigenIterDirSolvers::ConjugateGradient::setup"
    )   << nl << "Solver Eigen:PCG cannot be used with assymetric matrices. "            
        << endl
        << exit(FatalError);
  }
  
  if (precondType_ == "ICC")  
  {
    if (!reusePattern) 
      solverICC_->analyzePattern(matrix);   
    
    solverICC_->factorize(matrix);
    solverICC_->compute(matrix);
  }
  else
  {
    if (!reusePattern) 
      solverDiag_->analyzePattern(matrix);   
    
    solverDiag_->factorize(matrix);
    solverDiag_->compute(matrix);
  }
}

void Foam::EigenIterDirSolvers::ConjugateGradient::solve 
(
  const Eigen::VectorXd& source,
  Eigen::VectorXd& x
) 
{ 
  if (precondType_ == "ICC")  
  {
    x = solverICC_->solve(source);
  
    residual_ = solverICC_->error();
    nIters_ = solverICC_->iterations(); 
  }
  else
  {
    x = solverDiag_->solve(source);
  
    residual_ = solverDiag_->error();
    nIters_ = solverDiag_->iterations(); 
  }
  
}

 
// ************************************************************************* //
