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

#include "none.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJPreconditioners 
{
    defineTypeNameAndDebug(none, 0);
    addToRunTimeSelectionTable(HypreIJPreconditioner, none, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HypreIJPreconditioners::none::none
(
    const dictionary& dict
)
:
HypreIJPreconditioner(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HypreIJPreconditioners::none::initialize
(
  HYPRE_Solver& precond
) const
{
   
}

HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::none::setupPtr() const
{
   return nullptr;
}

 
HYPRE_PtrToSolverFcn Foam::HypreIJPreconditioners::none::solvePtr() const
{
   return nullptr;
}


void Foam::HypreIJPreconditioners::none::destroy
(
  HYPRE_Solver& solver
) const
{
 
}


// ************************************************************************* //
