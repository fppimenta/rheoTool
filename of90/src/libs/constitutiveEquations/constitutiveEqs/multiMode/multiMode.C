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

#include "multiMode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(multiMode, 0);
    addToRunTimeSelectionTable(constitutiveEq, multiMode, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::multiMode::multiMode
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            symmTensor::zero
        )
    ),
    models_(),
    etaS_("0", dimPressure*dimTime, 0.)
{
    PtrList<entry> modelEntries(dict.lookup("models"));
    models_.setSize(modelEntries.size());

    forAll (models_, modelI)
    {
        models_.set
        (
            modelI,
            constitutiveEq::New
            (
                word(name + modelEntries[modelI].keyword()),
                U,
                phi,
                modelEntries[modelI].dict()
            )
        );
        
        dimensionedScalar etaI(modelEntries[modelI].dict().lookup("etaS"));
        etaS_.value() += etaI.value();	
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::constitutiveEqs::multiMode::hasThermo
(
) const
{
    constitutiveEqProtectedIntr cPI;
 
    bool hasTherm(true);
    forAll (models_, i)
        hasTherm = (hasTherm && cPI.hasThermo(models_(i)));

    return hasTherm;
}

Foam::tmp<Foam::volScalarField>
Foam::constitutiveEqs::multiMode::etaSThermo
() const
{
    constitutiveEqProtectedIntr cPI;

    tmp<volScalarField> etaSTh(cPI.etaSThermo(models_(0)));

    for (label i = 1; i < models_.size(); i++)
    {
        etaSTh.ref() += cPI.etaSThermo(models_(i));
    }

    return etaSTh;
}

Foam::tmp<Foam::volScalarField>
Foam::constitutiveEqs::multiMode::etaPThermo
() const
{
    constitutiveEqProtectedIntr cPI;

    tmp<volScalarField> etaPTh(cPI.etaPThermo(models_(0)));

    for (label i = 1; i < models_.size(); i++)
    {
        etaPTh.ref() += cPI.etaPThermo(models_(i));
    }

    return etaPTh;
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::constitutiveEqs::multiMode::divTau
(
  const volVectorField& U
) const
{
    tmp<fvVectorMatrix> divMatrix = models_[0].divTau(U);

    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix.ref() += models_[i].divTau(U);
    }

    return divMatrix;
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::constitutiveEqs::multiMode::divTauS
(
  const volVectorField& U, 
  const volScalarField& alpha
) const
{
    tmp<fvVectorMatrix> divMatrix = models_[0].divTauS(U, alpha);

    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix.ref() += models_[i].divTauS(U, alpha);
    }

    return divMatrix;
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::constitutiveEqs::multiMode::divTauThermo
(
  const volVectorField& U
) const
{
    if (!this->hasThermo())
     NotImplemented;
    
    tmp<fvVectorMatrix> divMatrix = models_[0].divTauThermo(U);

    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix.ref() += models_[i].divTauThermo(U);
    }

    return divMatrix;
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::constitutiveEqs::multiMode::divTauSThermo
(
  const volVectorField& U, 
  const volScalarField& alpha
) const
{
    if (!this->hasThermo())
     NotImplemented;
     
    tmp<fvVectorMatrix> divMatrix = models_[0].divTauSThermo(U, alpha);

    for (label i = 1; i < models_.size(); i++)
    {
        divMatrix.ref() += models_[i].divTauSThermo(U, alpha);
    }

    return divMatrix;
}


Foam::tmp<Foam::volSymmTensorField> Foam::constitutiveEqs::multiMode::tau() const
{
    tau_ *= 0;

    for (label i = 0; i < models_.size(); i++)
    {
        tau_ += models_[i].tau();
    }

    return tau_;
}

const Foam::dimensionedScalar Foam::constitutiveEqs::multiMode::rho() const
{
    // It is unlikely (wrong) to have different densities,
    // but average to be sure.
    
    dimensionedScalar rho_( models_[0].rho() );

    int cnt(1);
    
    for (label i = 1; i < models_.size(); i++)
    {
        rho_ += models_[i].rho();
        cnt++;
    }

    return rho_/cnt;
}


void Foam::constitutiveEqs::multiMode::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
    forAll (models_, i)
    {
        Info<< "Model mode "  << i+1 << endl;
        models_[i].correct(alpha,gradU);
    }

    tau();
}

void Foam::constitutiveEqs::multiMode::divTauImplCoupled() const
{
  forAll (models_, i)
    models_[i].divTauImplCoupled();    
}


// ************************************************************************* //
