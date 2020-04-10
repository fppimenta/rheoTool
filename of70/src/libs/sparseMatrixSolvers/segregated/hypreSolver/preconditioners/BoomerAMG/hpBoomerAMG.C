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

#include "hpBoomerAMG.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJPreconditioners 
{
    defineTypeNameAndDebug(BoomerAMG, 0);
    addToRunTimeSelectionTable(HypreIJPreconditioner, BoomerAMG, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJPreconditioners::BoomerAMG::BoomerAMG
(
    const dictionary& dict
)
:
HypreIJPreconditioner(dict)
{
   
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJPreconditioners::BoomerAMG::initialize
(
  HYPRE_Solver& precond
) const
{
   HYPRE_BoomerAMGCreate(&precond);
  
   HYPRE_BoomerAMGSetPrintLevel(precond, dict_.lookupOrDefault<int>("printLevel", 0));  
   HYPRE_BoomerAMGSetTol(precond, dict_.lookupOrDefault<scalar>("relTol", 0.)); 
   HYPRE_BoomerAMGSetConvergeType(precond, dict_.lookupOrDefault<int>("convergenceType", 0));
   HYPRE_BoomerAMGSetMaxIter(precond, dict_.lookupOrDefault<int>("maxIter", 10));  
   HYPRE_BoomerAMGSetMinIter(precond, dict_.lookupOrDefault<int>("minIter", 0));  
   HYPRE_BoomerAMGSetRelaxType(precond, dict_.lookupOrDefault<int>("relaxationType", 6)); 
  
   HYPRE_BoomerAMGSetNumSweeps(precond, dict_.lookupOrDefault<int>("nSweeps", 1));
   HYPRE_BoomerAMGSetCycleType(precond, dict_.lookupOrDefault<int>("cycleType", 1));
   HYPRE_BoomerAMGSetVariant(precond, dict_.lookupOrDefault<int>("AMGVariant", 0)); 
   HYPRE_BoomerAMGSetRelaxOrder(precond, dict_.lookupOrDefault<int>("relaxOrder", 0)); 
   HYPRE_BoomerAMGSetMaxLevels(precond, dict_.lookupOrDefault<int>("maxLevels", 25)); 
   HYPRE_BoomerAMGSetCoarsenType(precond, dict_.lookupOrDefault<int>("coarsenType", 10));
   HYPRE_BoomerAMGSetRestriction(precond, dict_.lookupOrDefault<int>("restrictionType", 0));
   
   HYPRE_BoomerAMGSetLogging(precond, dict_.lookupOrDefault<int>("logging", 0));    
   
   if (dict_.lookupOrDefault<bool>("recoversOldDefault", true))
     HYPRE_BoomerAMGSetOldDefault(precond);   
}

extern "C" HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::BoomerAMG::setupPtr() const
{
   return (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup;
}

 
extern "C" HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::BoomerAMG::solvePtr() const
{
   return (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve;
}


void Foam::HypreIJPreconditioners::BoomerAMG::destroy
(
  HYPRE_Solver& solver
) const
{
  HYPRE_BoomerAMGDestroy(solver);
}


// ************************************************************************* //
