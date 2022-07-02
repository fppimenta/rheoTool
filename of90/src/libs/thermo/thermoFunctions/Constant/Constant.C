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

#include "Constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace thermoFunctions 
 {
   defineTypeNameAndDebug(Constant, 0);
   addToRunTimeSelectionTable(thermoFunction, Constant, dictionary);
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoFunctions::Constant::Constant
(
  const word& name,
  const fvMesh& mesh,
  const dictionary& dict
)
:
  thermoFunction(name, mesh)
{
  // T might not exist. Keep nullptr but do not dereference.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermoFunctions::Constant::multiply
(
  volScalarField& vf,
  const dimensionedScalar& ds
) const
{
  // Do nothing
}

void Foam::thermoFunctions::Constant::multiply
(
  volScalarField& vf
) const
{
  // Do nothing
}

Foam::tmp<Foam::volScalarField> Foam::thermoFunctions::Constant::createField
(
  const dimensionedScalar& ds
) const
{
  tmp<volScalarField> tvf
  (
   new volScalarField 
   (
     IOobject
     (
       "tmpField",
       mesh().time().timeName(),
       mesh(),
       IOobject::NO_READ,
       IOobject::NO_WRITE,
       false
     ),
     mesh(),
     ds,
     extrapolatedCalculatedFvPatchField<scalar>::typeName
   )
  );
    
  return tvf; 
}

// ************************************************************************* //
