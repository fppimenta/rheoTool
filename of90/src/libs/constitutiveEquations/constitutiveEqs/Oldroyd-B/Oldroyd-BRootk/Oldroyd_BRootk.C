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

#include "Oldroyd_BRootk.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(Oldroyd_BRootk, 0);
    addToRunTimeSelectionTable(constitutiveEq, Oldroyd_BRootk, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::Oldroyd_BRootk::Oldroyd_BRootk
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    theta_
    (
        IOobject
        (
            "theta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    eigVals_
    (
        IOobject
        (
            "eigVals" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
         extrapolatedCalculatedFvPatchField<tensor>::typeName
    ),
    eigVecs_
    (
        IOobject
        (
            "eigVecs" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
         extrapolatedCalculatedFvPatchField<tensor>::typeName
    ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    lambda_(dict.lookup("lambda")),
    k_(dict.lookup("k")),
    thermoLambdaPtr_(thermoFunction::New("thermoLambda", U.mesh(), dict)),  
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))  
{
 checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::Oldroyd_BRootk::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
    // Update temperature-dependent properties
    volScalarField lambda(thermoLambdaPtr_->createField(lambda_));
    volScalarField etaP(thermoEtaPtr_->createField(etaP_));
 
    // Decompose grad(U).T()
    #include "boilerLog.H"

    scalar k(k_.value());
  
    volTensorField opert
    (
      IOobject
      (
        "opert",
        U().time().timeName(),
        U().mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
      ),
      U().mesh(),
      dimensionedTensor
      (
       "zero",
       dimless,
       pTraits<tensor>::zero
      ),
      extrapolatedCalculatedFvPatchField<tensor>::typeName
    );

    forAll(opert, i)
    {
        opert[i].xx() = Foam::pow(eigVals_[i].xx(), 1./k-1.);
        opert[i].yy() = Foam::pow(eigVals_[i].yy(), 1./k-1.);
        opert[i].zz() = Foam::pow(eigVals_[i].zz(), 1./k-1.);

        opert[i] = (eigVecs_[i] & opert[i] & eigVecs_[i].T()) - theta_[i];
    }

    fvSymmTensorMatrix thetaEqn
    (
         fvm::ddt(theta_)
       + fvm::div(phi(), theta_)
       ==
       symm
       (  
         (omega&theta_)
       - (theta_&omega)
       + (2.0/k_) * (B&theta_)
       + (1.0/(lambda*k)) * opert         
       ) 
    );

   
    thetaEqn.relax();
    thetaEqn.solve();
  
    // Diagonalization of theta

    calcEig(theta_, eigVals_, eigVecs_);
    forAll(theta_, cellI)
    {
       eigVals_[cellI].xx()=Foam::pow(Foam::log(eigVals_[cellI].xx()), k);
       eigVals_[cellI].yy()=Foam::pow(Foam::log(eigVals_[cellI].yy()), k);
       eigVals_[cellI].zz()=Foam::pow(Foam::log(eigVals_[cellI].zz()), k);
    }

    // Convert from theta to tau

    tau_ = (etaP/lambda) * symm( (eigVecs_ & eigVals_ & eigVecs_.T()) - Itensor);

    tau_.correctBoundaryConditions();

}


// ************************************************************************* //
