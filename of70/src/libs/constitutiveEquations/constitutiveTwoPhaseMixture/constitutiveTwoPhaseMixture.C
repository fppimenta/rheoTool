/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "constitutiveTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constitutiveTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveTwoPhaseMixture::constitutiveTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "constitutiveProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(U.mesh(), *this),
    phase1_
    (
        constitutiveEq::New(
                             word("." + phase1Name_),
                             U,
                             phi,
                             subDict(phase1Name_ == "1" ? "phase1": phase1Name_).subDict("parameters")
                           )
    ), 
    phase2_
    (
        constitutiveEq::New(
                             word("." + phase2Name_),
                             U,
                             phi,
                             subDict(phase2Name_ == "2" ? "phase2": phase2Name_).subDict("parameters")
                           )
    ),
    rho1_(phase1_->rho()),
    rho2_(phase2_->rho()),
    U_(U),
    phi_(phi),
    tauMF_
    (
        IOobject
        (
            "tauMF",
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
        calculatedFvPatchScalarField::typeName
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix>
Foam::constitutiveTwoPhaseMixture::divTauMF(volVectorField& U)
{

   const volScalarField bAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    
   const volScalarField bAlpha2(1.-bAlpha1);
   
   tauMF_ = phase1_->tau()*bAlpha1 + phase2_->tau()*bAlpha2;
     
   if (!(phase1_->isGNF() && phase2_->isGNF()) && U.time().outputTime())
    {
        tauMF_.write();     
    }
    
   return (
             phase1_->divTauS(U, bAlpha1) 
           + phase2_->divTauS(U, bAlpha2)
           + fvc::div(tauMF_,"div(Sum(tau))")          
         );
         
}

Foam::tmp<Foam::volSymmTensorField>
Foam::constitutiveTwoPhaseMixture::tauTotalMF() const
{

   const volScalarField bAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );
    
   const volScalarField bAlpha2(1.-bAlpha1);
   
   return phase1_->tauTotal()*bAlpha1 + phase2_->tauTotal()*bAlpha2;
}

bool Foam::constitutiveTwoPhaseMixture::read()
{

    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
        
}

// ************************************************************************* //
