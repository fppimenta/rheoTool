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

#include "hsBiCGSTAB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJSolvers 
{
    defineTypeNameAndDebug(BiCGSTAB, 0);
    addToRunTimeSelectionTable(HypreIJSolver, BiCGSTAB, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJSolvers::BiCGSTAB::BiCGSTAB
(
    const dictionary& dict
)
:
HypreIJSolver(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJSolvers::BiCGSTAB::initialize 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);
  
  HYPRE_BiCGSTABSetMaxIter(solver, dict_.lookupOrDefault<int>("maxIter", 1000)); /* max iterations */
  HYPRE_BiCGSTABSetTol(solver, dict_.lookupOrDefault<scalar>("relTol", 0)); /* conv. tolerance */
  HYPRE_BiCGSTABSetAbsoluteTol(solver, dict_.lookupOrDefault<scalar>("tolerance", 1e-8)); /* conv. tolerance */
  HYPRE_BiCGSTABSetPrintLevel(solver, dict_.lookupOrDefault<int>("printLevel", 0)); /* print solve info */
  HYPRE_BiCGSTABSetLogging(solver, dict_.lookupOrDefault<int>("logging", 0)); /* needed to get run info later */
}


void Foam::HypreIJSolvers::BiCGSTAB::setup  
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
    HYPRE_BiCGSTABSetPrecond(solver, precondSolvePtr, precondSetupPtr, precond);
         
  HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x); 
}

void Foam::HypreIJSolvers::BiCGSTAB::solve 
(
  HYPRE_Solver& solver,
  HYPRE_ParCSRMatrix& parcsr_A,
  HYPRE_ParVector& par_b, 
  HYPRE_ParVector& par_x
) 
{
     
  HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x);

  // Get stats
  HYPRE_BiCGSTABGetNumIterations(solver, &nIters_);
  HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &residual_);
  
}

void Foam::HypreIJSolvers::BiCGSTAB::destroy 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_ParCSRBiCGSTABDestroy(solver);
}

// ************************************************************************* //
