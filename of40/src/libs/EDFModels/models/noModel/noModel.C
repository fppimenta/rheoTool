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

#include "noModel.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noModel, 0);
    addToRunTimeSelectionTable(EDFEquation, noModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noModel::noModel
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    EDFEquation(name, phi),
    zeroFe_
    (
        IOobject
        (
            "zeroForce",
            phi.time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi.mesh(),
        dimensionedVector
        (
                "zero",
                dimVelocity*dimDensity/dimTime,
                pTraits<vector>::zero
        ),
        extrapolatedCalculatedFvPatchField<vector>::typeName
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::noModel::Fe() const
{
    return
    (
         zeroFe_ 
    );     
}

void Foam::noModel::correct()
{

       // Do nothing
}

// ************************************************************************* //
