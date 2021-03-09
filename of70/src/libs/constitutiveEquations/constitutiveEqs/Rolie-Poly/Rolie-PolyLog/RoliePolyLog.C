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

#include "RoliePolyLog.H"
#include "addToRunTimeSelectionTable.H"
 
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(RoliePolyLog, 0);
    addToRunTimeSelectionTable(constitutiveEq, RoliePolyLog, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::RoliePolyLog::RoliePolyLog
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
    lambdaR_(dict.lookup("lambdaR")),
    lambdaD_(dict.lookup("lambdaD")),
    beta_(dict.lookup("beta")),
    delta_(dict.lookup("delta")),
    chiMax_(dict.lookup("chiMax")),
    thermoLambdaRPtr_(thermoFunction::New("thermoLambdaR", U.mesh(), dict)),  
    thermoLambdaDPtr_(thermoFunction::New("thermoLambdaD", U.mesh(), dict)),  
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))  
{
 checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::RoliePolyLog::correct()
{
    // Update temperature-dependent properties
    volScalarField lambdaR(thermoLambdaRPtr_->createField(lambdaR_));
    volScalarField lambdaD(thermoLambdaDPtr_->createField(lambdaD_));
    volScalarField etaP(thermoEtaPtr_->createField(etaP_));
 
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
    
    volSymmTensorField A_(symm(eigVecs_ & eigVals_ & eigVecs_.T()));

    volScalarField trA(tr(A_));
    
    volScalarField M1
     ( 
       2.*(1.-Foam::sqrt(3./trA))/lambdaR
     );
     
    if (chiMax_.value() > 1.)
     {
    
       M1 *= 
            (
               (  3.-(trA/3.)/sqr(chiMax_) ) * (1.-1./sqr(chiMax_))
             / ( (1.-(trA/3.)/sqr(chiMax_) ) * (3.-1./sqr(chiMax_)  ) ) 
            ); 
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
       + 2.0 * B
       - (1.0/lambdaD) * 
         ( 
           (
              eigVecs_ &
                (
                  inv(eigVals_)  
                )
            & eigVecs_.T() 
           ) &
           ((A_-Itensor) + M1 * lambdaD * ( A_ + beta_*pow(trA/3., delta_)*(A_-Itensor) ))
          )
       ) 
       
    );
    

    thetaEqn.relax();
    thetaEqn.solve();
  
    // Diagonalization of theta

    calcEig(theta_, eigVals_, eigVecs_);

    // Convert from theta to tau
    
    A_ = symm(eigVecs_ & eigVals_ & eigVecs_.T());

    tau_ = (etaP/lambdaD) * ( A_ - symm(Itensor));
    
    if (chiMax_.value() > 1.)
     {
         
       trA = tr(A_); //Update
        
       tau_ *= 
            (
               (  3.-(trA/3.)/sqr(chiMax_) ) * (1.-1./sqr(chiMax_))
             / ( (1.-(trA/3.)/sqr(chiMax_) ) * (3.-1./sqr(chiMax_)  ) ) 
            );
     }

    tau_.correctBoundaryConditions();

}


// ************************************************************************* //
