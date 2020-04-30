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

#include "ParaSails.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJPreconditioners 
{
    defineTypeNameAndDebug(ParaSails, 0);
    addToRunTimeSelectionTable(HypreIJPreconditioner, ParaSails, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJPreconditioners::ParaSails::ParaSails
(
    const dictionary& dict
)
:
HypreIJPreconditioner(dict)
{
   if (!Pstream::parRun())
   {
     FatalErrorIn
        (
            "Foam::HypreIJPreconditioners::ParaSails::ParaSails()"
        )   << nl << "Hypre preconditioner ParaSails cannot be run in serial. "
            << "Either perform a parallel run or choose other preconditioner." << endl
            << exit(FatalError);
   }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJPreconditioners::ParaSails::initialize
(
  HYPRE_Solver& precond
) const
{
   HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);
   HYPRE_ParaSailsSetSym(precond, dict_.lookupOrDefault<int>("printLevel", 0));
   HYPRE_ParaSailsSetLoadbal(precond, dict_.lookupOrDefault<scalar>("loadBalance", 0));
   
}

extern "C" HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::ParaSails::setupPtr() const
{
   return (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup; 
}

 
extern "C" HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::ParaSails::solvePtr() const
{
   return (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve;
}


void Foam::HypreIJPreconditioners::ParaSails::destroy
(
  HYPRE_Solver& solver
) const
{
  HYPRE_ParaSailsDestroy(solver);
}


// ************************************************************************* //
