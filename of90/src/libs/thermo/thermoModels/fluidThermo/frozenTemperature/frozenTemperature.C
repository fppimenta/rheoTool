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

#include "frozenTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace fluidThermoModels
 {
   defineTypeNameAndDebug(frozenTemperature, 0);
   addToRunTimeSelectionTable(fluidThermoModel, frozenTemperature, dictionary);
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermoModels::frozenTemperature::frozenTemperature
(
  const word& name,
  const fvMesh& mesh
)
:
  fluidThermoModel(name, mesh)
{
  
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::fluidThermoModels::frozenTemperature::correct
(
  const volVectorField& U,
  const surfaceScalarField& phi,
  const volSymmTensorField& tau,
  autoPtr<coupledSolver>& cpsT,
  fvModels& fvModel,
  fvConstraints& fvConstraint,
  int nNonOrtoC
)
{
   // Do nothing: T is constant over time
}
 

 

// ************************************************************************* //
