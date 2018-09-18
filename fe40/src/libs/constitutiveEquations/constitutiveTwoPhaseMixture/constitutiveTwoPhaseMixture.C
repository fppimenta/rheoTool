/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "constitutiveTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

word constitutiveTwoPhaseMixture::getPhaseName(const word& key) const
{
    if (isDict(key))
    {
        return key;
    }
    else
    {
        return word(lookup(key)); 
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveTwoPhaseMixture::constitutiveTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word& alpha1Name
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
    phase1Name_(getPhaseName("phase1")),
    phase2Name_(getPhaseName("phase2")),
    phase1_
    (
        constitutiveEq::New(
                             word("1"),
                             U,
                             phi,
                             subDict("phase1").subDict("parameters")
                           )
    ),
    phase2_
    (
        constitutiveEq::New(
                             word("2"),
                             U,
                             phi,
                             subDict("phase2").subDict("parameters")
                           )
    ),
    rho1_(phase1_->rho()), 
    rho2_(phase2_->rho()),
    U_(U),
    phi_(phi),
    alpha1_(U_.db().lookupObject<const volScalarField> (alpha1Name)),
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
        calculatedFvPatchField<symmTensor>::typeName
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
    
   const volScalarField bAlpha2(1.0-bAlpha1);
   
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
