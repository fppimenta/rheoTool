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

#include "HypreIJSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HypreIJSolver, 0);
    defineRunTimeSelectionTable(HypreIJSolver, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HypreIJSolver::HypreIJSolver
(
    const dictionary& dict
)
:
dict_(dict),
nIters_(0),
residual_(1.)
{
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HypreIJSolver::HypreIJSolver::printInfo 
(
  word tName,
  word pcName,
  scalar initResidual,
  scalar finalResidual
) const
{
  if (pcName == " ")
  {
    
    Info << "Hypre:" << solverName() << ": Solving for "<< tName  << ", Initial residual = " << initResidual 
         << " Final residual = " << finalResidual << ", No Iterations " << nIterations() << endl; 
  }
  else
  {
    
    Info << "Hypre:" << solverName() << ":" << pcName << ": Solving for "<< tName  << ", Initial residual = " << initResidual  
         << " Final residual = " << finalResidual << ", No Iterations " << nIterations() << endl; 
  }
}
 

} //End namespace

// ************************************************************************* //
