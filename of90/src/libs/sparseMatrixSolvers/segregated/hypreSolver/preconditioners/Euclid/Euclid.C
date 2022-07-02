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

#include "Euclid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJPreconditioners 
{
    defineTypeNameAndDebug(Euclid, 0);
    addToRunTimeSelectionTable(HypreIJPreconditioner, Euclid, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJPreconditioners::Euclid::Euclid
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
            "Foam::HypreIJPreconditioners::Euclid::Euclid()"
        )   << nl << "Hypre preconditioner Euclid cannot be run in serial. "
            << "Either perform a parallel run or choose other preconditioner." << endl
            << exit(FatalError);
   }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJPreconditioners::Euclid::initialize
(
  HYPRE_Solver& precond
) const
{
   HYPRE_EuclidCreate(MPI_COMM_WORLD, &precond);
   
   HYPRE_EuclidSetLevel(precond, dict_.lookupOrDefault<int>("ILUklevel", 2));
   HYPRE_EuclidSetBJ(precond, dict_.lookupOrDefault<bool>("enableJacobiILU", false));
   HYPRE_EuclidSetSparseA(precond, dict_.lookupOrDefault<scalar>("dropTol", 0));
   HYPRE_EuclidSetRowScale(precond, dict_.lookupOrDefault<bool>("enableRowScaling", false));
}

extern "C" HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::Euclid::setupPtr() const
{
   return (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup;
}

 
extern "C" HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::Euclid::solvePtr() const
{
   return (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve;
}


void Foam::HypreIJPreconditioners::Euclid::destroy
(
  HYPRE_Solver& solver
) const
{
  HYPRE_EuclidDestroy(solver);
}


// ************************************************************************* //
