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

#include "HerschelBulkley.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HerschelBulkley, 0);
    addToRunTimeSelectionTable(constitutiveEq, HerschelBulkley, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HerschelBulkley::HerschelBulkley
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi),
    rho_(dict.lookup("rho")),
    tau0_(dict.lookup("tau0")),
    eta0_(dict.lookup("eta0")),
    k_(dict.lookup("k")),
    n_(dict.lookup("n")),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
                "0",
                dimensionSet(1, -1, -2, 0, 0, 0, 0),
                pTraits<symmTensor>::zero  
        ),
        extrapolatedCalculatedFvPatchField<symmTensor>::typeName
    ),
    eta_
    (
        IOobject
        (
            "eta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        strainRate()*dimensionedScalar("zeroU", dimensionSet(1, -1, 0, 0, 0, 0, 0), 0) //Just to ensure dimensions and BCs
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::HerschelBulkley::correct()
{
    volScalarField strainRate_(strainRate());
    
    dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    eta_ = 
        min
        (
             eta0_,
            (tau0_ + k_*rtone*pow(tone*strainRate_, n_))
           /(max(strainRate_, dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
        );
}


// ************************************************************************* //
