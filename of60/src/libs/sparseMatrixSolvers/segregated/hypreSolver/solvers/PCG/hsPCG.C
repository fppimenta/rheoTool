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

#include "hsPCG.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJSolvers 
{
    defineTypeNameAndDebug(PCG, 0);
    addToRunTimeSelectionTable(HypreIJSolver, PCG, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJSolvers::PCG::PCG
(
    const dictionary& dict
)
:
HypreIJSolver(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJSolvers::PCG::initialize 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
  
  HYPRE_PCGSetMaxIter(solver, dict_.lookupOrDefault<int>("maxIter", 1000)); /* max iterations */
  HYPRE_PCGSetTol(solver, dict_.lookupOrDefault<scalar>("relTol", 0)); /* conv. tolerance */
  HYPRE_PCGSetAbsoluteTol(solver, dict_.lookupOrDefault<scalar>("tolerance", 1e-8)); /* conv. tolerance */
  HYPRE_PCGSetTwoNorm(solver, dict_.lookupOrDefault<bool>("useTwoNorm", true)); /* use two norm check */
  HYPRE_PCGSetRecomputeResidual(solver, dict_.lookupOrDefault<bool>("recomputeEndResidual", false)); /* recompute end residual */
  HYPRE_PCGSetPrintLevel(solver, dict_.lookupOrDefault<int>("printLevel", 0)); /* print solve info */
  HYPRE_PCGSetLogging(solver, dict_.lookupOrDefault<int>("logging", 0)); /* needed to get run info later */
}


void Foam::HypreIJSolvers::PCG::setup  
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
  if (!isSymmetric)
  {
    FatalErrorIn
    (
      "Foam::HypreIJSolvers::PCG::setup"
    )   << nl << "Solver Hypre:PCG cannot be used with assymetric matrices. "            
        << endl
        << exit(FatalError);
  }
  
  if (precondSolvePtr != 0)
    HYPRE_PCGSetPrecond(solver, precondSolvePtr, precondSetupPtr, precond);
         
  HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x); 
}

void Foam::HypreIJSolvers::PCG::solve 
(
  HYPRE_Solver& solver,
  HYPRE_ParCSRMatrix& parcsr_A,
  HYPRE_ParVector& par_b, 
  HYPRE_ParVector& par_x
)
{

  HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

  // Get stats
  HYPRE_PCGGetNumIterations(solver, &nIters_);
  HYPRE_PCGGetFinalRelativeResidualNorm(solver, &residual_); 
}


void Foam::HypreIJSolvers::PCG::destroy 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_ParCSRPCGDestroy(solver);
}
 
// ************************************************************************* //
