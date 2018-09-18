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

#include "EDFEquation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(EDFEquation, 0);
    defineRunTimeSelectionTable(EDFEquation, dictionary);
}

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

const Foam::dimensionedScalar
Foam::EDFEquation::epsilonK_
( 
        "vacuumPermittivity", 
        dimensionSet(-1, -3, 4, 0, 0, 2, 0),
        8.8541878176e-12 
);

const Foam::dimensionedScalar
Foam::EDFEquation::AK_
( 
        "AvogradoNumber", 
        dimensionSet(0, 0, 0, 0, -1, 0, 0),
        6.022140857e+23
);

const Foam::dimensionedScalar
Foam::EDFEquation::eK_
( 
        "elementaryCharge", 
        dimensionSet( 0, 0, 1, 0, 0, 1, 0 ),
        1.6021766208e-19
);

const Foam::dimensionedScalar
Foam::EDFEquation::kbK_
( 
        "BoltzmannConstant", 
        dimensionSet( 1, 2, -2, -1, 0, 0, 0 ),
        1.38064852e-23
);

const Foam::dimensionedScalar
Foam::EDFEquation::FK_
( 
       "FaradayConstant", 
       (Foam::EDFEquation::eK_ * Foam::EDFEquation::AK_).dimensions(),
       (Foam::EDFEquation::eK_ * Foam::EDFEquation::AK_).value() 
); 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquation::EDFEquation
(
    const word& name,
    const surfaceScalarField& phi
)
:
    name_(name),
    phi_(phi)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * //

bool Foam::EDFEquation::checkForPhiE
(
 const word& name,
 const surfaceScalarField& phi
)
{

 IOobject phiEheader
 (
       "phiE" + name,
       phi.time().timeName(),
       phi.mesh(),
       IOobject::MUST_READ
 );
 
 if (phiEheader.typeHeaderOk<volScalarField>(true))
 {         
   Info<< "Field phiE found. Electric potential split: psi and phiE. \n" << endl;
   return true;         
 }
 else
 {
   Info<< "Field phiE not found. Single electric potential: psi. \n" << endl;
   return false;   
 }
     
}


// ************************************************************************* //
