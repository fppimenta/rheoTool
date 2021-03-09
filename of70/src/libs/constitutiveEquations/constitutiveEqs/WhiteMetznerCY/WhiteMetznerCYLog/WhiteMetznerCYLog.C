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

#include "WhiteMetznerCYLog.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(WhiteMetznerCYLog, 0);
    addToRunTimeSelectionTable(constitutiveEq, WhiteMetznerCYLog, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::WhiteMetznerCYLog::WhiteMetznerCYLog
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
    m_(dict.lookup("m")),
    n_(dict.lookup("n")),
    K_(dict.lookup("K")),
    L_(dict.lookup("L")),
    a_(dict.lookup("a")),
    b_(dict.lookup("b")),
    thermoLambdaPtr_(thermoFunction::New("thermoLambda", U.mesh(), dict)),  
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))  
{
 checkForStab(dict);
 
 // Check if parameters allow the use of the log version
 
 if ( m_.value()!=n_.value() || K_.value()!=L_.value() || a_.value()!=b_.value())
  {
     FatalErrorIn("Foam::constitutiveEqs::WhiteMetznerCYLog::WhiteMetznerCYLog\n")
            << "The Log version of the WhiteMetznerCY model can only be used if:\n"
            << "\n   m=n   and   K=L   and   a=b\n"
            << "\n Check if this is the case, otherwise use the stress version (WhiteMetznerCY).\n" 
            << abort(FatalError);
  }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::WhiteMetznerCYLog::correct()
{
    // Decompose grad(U).T()

    volTensorField L(fvc::grad(U()));

    dimensionedScalar c1( "zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.);
    volTensorField   B(c1 * eigVecs_); 
    volTensorField   omega(B);
    volTensorField   M(eigVecs_.T() & L.T() & eigVecs_);

    decomposeGradU(M, eigVals_, eigVecs_, omega, B);

  
    // Solve the constitutive Eq in theta = log(c)
 
    dimensionedTensor Itensor
    ( 
        "Identity", 
        dimensionSet(0, 0, 0, 0, 0, 0, 0), 
        tensor::I 
    );

    // Effective viscosity and relaxation time
    volScalarField etaP(etaP_*
        Foam::pow(1 + Foam::pow(K_* sqrt(2.0)*mag(symm(L)),a_), (n_- 1)/a_));

    volScalarField lambda(lambda_*
        Foam::pow(1 + Foam::pow( L_* sqrt(2.0)*mag(symm(L)),b_), (m_- 1)/b_));

    // Update temperature-dependent properties
    thermoLambdaPtr_->multiply(lambda);
    thermoEtaPtr_->multiply(etaP);
    
    fvSymmTensorMatrix thetaEqn
    (
         fvm::ddt(theta_)
       + fvm::div(phi(), theta_)
       ==
       symm
       (  
         (omega&theta_)
       - (theta_&omega)
       + 2.0 * B
       + (1.0/lambda) 
        * (
             eigVecs_ &
               (
                  inv(eigVals_)  
                - Itensor
               )
           & eigVecs_.T() 
          )
       ) 
       
    );

   
    thetaEqn.relax();
    thetaEqn.solve();
  
    // Diagonalization of theta

    calcEig(theta_, eigVals_, eigVecs_);
 
    // Convert from theta to tau

    tau_ = (etaP/lambda) * symm( (eigVecs_ & eigVals_ & eigVecs_.T()) - Itensor);

    tau_.correctBoundaryConditions();

}


// ************************************************************************* //
