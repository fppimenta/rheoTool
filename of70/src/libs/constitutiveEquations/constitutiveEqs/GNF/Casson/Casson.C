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

#include "Casson.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(Casson, 0);
    addToRunTimeSelectionTable(constitutiveEq, Casson, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::Casson::Casson
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi),
    rho_(dict.lookup("rho")),
    reg_(dict.lookup("PapanastasiouRegularization")),
    tau0_(dict.lookup("tau0")),
    etaInf_(dict.lookup("etaInf")),
    etaMin_(reg_ == false ? dict.lookup("etaMin") : dimensionedScalar("1",dimPressure*dimTime,1.)),
    etaMax_(reg_ == false ? dict.lookup("etaMax") : dimensionedScalar("1",dimPressure*dimTime,1.)),
    m_(reg_ == true ? dict.lookup("m") : dimensionedScalar("1",dimTime,1.) ),
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
    ),
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))
{
  // Initialize eta_
  correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::Casson::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
  
 if (reg_)
  {
    volScalarField
    strainRate_ =
    max
    (
        strainRate(gradU),
        dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
    );
    
    eta_ = sqr( sqrt(tau0_/strainRate_)*( 1. - exp(-sqrt(m_*strainRate_)) ) + sqrt(etaInf_) );  
  }
 else 
  {
    eta_ = 
    max(
        etaMin_,
        min
        (
            etaMax_,
            sqr
            (
                sqrt
                (
                   tau0_
                   /max
                    (
                        strainRate(),
                        dimensionedScalar("VSMALL", dimless/dimTime, VSMALL)
                    )
                ) + sqrt(etaInf_)
            )
         )
       ); 
  }
  
  // Update viscosity for non-isothermal flow
  thermoEtaPtr_->multiply(eta_);
         
}


// ************************************************************************* //
