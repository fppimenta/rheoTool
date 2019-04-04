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

#include "hsGMRES.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJSolvers 
{
    defineTypeNameAndDebug(GMRES, 0);
    addToRunTimeSelectionTable(HypreIJSolver, GMRES, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJSolvers::GMRES::GMRES
(
    const dictionary& dict
)
:
HypreIJSolver(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJSolvers::GMRES::initialize 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
  
  HYPRE_GMRESSetMaxIter(solver, dict_.lookupOrDefault<int>("maxIter", 1000)); /* max iterations */
  HYPRE_GMRESSetTol(solver, dict_.lookupOrDefault<scalar>("relTol", 0)); /* conv. tolerance */
  HYPRE_GMRESSetAbsoluteTol(solver, dict_.lookupOrDefault<scalar>("tolerance", 1e-8)); /* conv. tolerance */
  HYPRE_GMRESSetKDim(solver, dict_.lookupOrDefault<int>("KrylovSpaceDim", 100)); /* size of Krylov space */
  HYPRE_GMRESSetPrintLevel(solver, dict_.lookupOrDefault<int>("printLevel", 0)); /* print solve info */
  HYPRE_GMRESSetLogging(solver, dict_.lookupOrDefault<int>("logging", 0)); /* needed to get run info later */
}


void Foam::HypreIJSolvers::GMRES::setup  
(
  HYPRE_Solver& solver,
  HYPRE_Solver& precond,
  HYPRE_ParCSRMatrix& parcsr_A,
  HYPRE_ParVector& par_b, 
  HYPRE_ParVector& par_x,
  HYPRE_PtrToSolverFcn precondSolvePtr,
  HYPRE_PtrToSolverFcn precondSetupPtr,
  bool isSymmetric
) const
{
  if (precondSolvePtr != 0)
    HYPRE_GMRESSetPrecond(solver, precondSolvePtr, precondSetupPtr, precond);
         
  HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b, par_x); 
}

void Foam::HypreIJSolvers::GMRES::solve 
(
  HYPRE_Solver& solver,
  HYPRE_ParCSRMatrix& parcsr_A,
  HYPRE_ParVector& par_b, 
  HYPRE_ParVector& par_x
) 
{
    
  HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b, par_x);

  // Get stats
  HYPRE_GMRESGetNumIterations(solver, &nIters_);
  HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &residual_);
   
}

 
void Foam::HypreIJSolvers::GMRES::destroy 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_ParCSRGMRESDestroy(solver);
}
 
// ************************************************************************* //
