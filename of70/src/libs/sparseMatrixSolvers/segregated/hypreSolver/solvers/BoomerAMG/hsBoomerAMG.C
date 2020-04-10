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

#include "hsBoomerAMG.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJSolvers 
{
    defineTypeNameAndDebug(BoomerAMG, 0);
    addToRunTimeSelectionTable(HypreIJSolver, BoomerAMG, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJSolvers::BoomerAMG::BoomerAMG
(
    const dictionary& dict
)
:
HypreIJSolver(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJSolvers::BoomerAMG::initialize 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_BoomerAMGCreate(&solver);
  
  HYPRE_BoomerAMGSetPrintLevel(solver, dict_.lookupOrDefault<int>("printLevel", 0));  
  HYPRE_BoomerAMGSetTol(solver, dict_.lookupOrDefault<scalar>("relTol", 0.)); 
  HYPRE_BoomerAMGSetConvergeType(solver, dict_.lookupOrDefault<int>("convergenceType", 0));
  HYPRE_BoomerAMGSetMaxIter(solver, dict_.lookupOrDefault<int>("maxIter", 1000));  
  HYPRE_BoomerAMGSetMinIter(solver, dict_.lookupOrDefault<int>("minIter", 0));  
  HYPRE_BoomerAMGSetRelaxType(solver, dict_.lookupOrDefault<int>("relaxationType", 6)); 
 
  HYPRE_BoomerAMGSetNumSweeps(solver, dict_.lookupOrDefault<int>("nSweeps", 1));
  HYPRE_BoomerAMGSetCycleType(solver, dict_.lookupOrDefault<int>("cycleType", 1));
  HYPRE_BoomerAMGSetVariant(solver, dict_.lookupOrDefault<int>("AMGVariant", 0)); 
  HYPRE_BoomerAMGSetRelaxOrder(solver, dict_.lookupOrDefault<int>("relaxOrder", 0)); 
  HYPRE_BoomerAMGSetMaxLevels(solver, dict_.lookupOrDefault<int>("maxLevels", 25)); 
  HYPRE_BoomerAMGSetCoarsenType(solver, dict_.lookupOrDefault<int>("coarsenType", 10));
  HYPRE_BoomerAMGSetRestriction(solver, dict_.lookupOrDefault<int>("restrictionType", 0));  
  
  HYPRE_BoomerAMGSetLogging(solver, dict_.lookupOrDefault<int>("logging", 0));  
  
  if (dict_.lookupOrDefault<bool>("recoversOldDefault", true))
    HYPRE_BoomerAMGSetOldDefault(solver);   
}


void Foam::HypreIJSolvers::BoomerAMG::setup  
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
  // No preconditioning 
  
  HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
}

void Foam::HypreIJSolvers::BoomerAMG::solve 
(
  HYPRE_Solver& solver,
  HYPRE_ParCSRMatrix& parcsr_A,
  HYPRE_ParVector& par_b, 
  HYPRE_ParVector& par_x
) 
{
     
  HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

  // Get stats
  HYPRE_BoomerAMGGetNumIterations(solver, &nIters_);
  HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &residual_);
 
}


void Foam::HypreIJSolvers::BoomerAMG::destroy 
(
  HYPRE_Solver& solver
) const
{
  HYPRE_BoomerAMGDestroy(solver);
}

// ************************************************************************* //
