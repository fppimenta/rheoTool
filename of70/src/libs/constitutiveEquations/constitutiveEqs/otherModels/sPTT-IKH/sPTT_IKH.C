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

#include "sPTT_IKH.H" 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace constitutiveEqs 
 {
    defineTypeNameAndDebug(sPTT_IKH, 0);
    addToRunTimeSelectionTable(constitutiveEq, sPTT_IKH, dictionary);
 }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::constitutiveEqs::sPTT_IKH::sPTT_IKH
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
    A_
    (
        IOobject
        (
            "A" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ), 
    lambda_
    (
        IOobject
        (
            "lambda" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),   
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    eta0_(dict.lookup("eta0")),
    etaStab_(dict.lookup("etaStab")),
    G_(dict.lookup("G")),
    k0_(dict.lookup("k0")),
    eps_(dict.lookup("epsilon")),
    C0_(dict.lookup("C0")),
    q_(dict.lookup("q")),
    k1_(dict.lookup("k1")),
    k2_(dict.lookup("k2")),
    k3_(dict.lookup("k3")),
    n1_(dict.lookup("n1")),
    n2_(dict.lookup("n2")),
    n3_(dict.lookup("n3")),
    m1_(dict.lookup("m1")),
    m2_(dict.lookup("m2")),
    m3_(dict.lookup("m3"))
{
    // Stabilization
    checkForStab(dict);    
}
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::sPTT_IKH::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
    dimensionedScalar zeroD("0",dimless/eta0_.dimensions(),0.);
    
    dimensionedScalar small("S",dimPressure*eta0_.dimensions(),1e-16);
    
    dimensionedScalar oneDim("1",dimTime,1);

//- Solve for tau

    // Velocity gradient tensor
    volTensorField L( gradU == nullptr ? fvc::grad(U())() : *gradU );
    
    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);
   
    // Inter vars
    volSymmTensorField tauEff(tau_-C0_*Foam::pow(lambda_,m3_)*A_);
   
    volScalarField sigmaEff(Foam::mag(dev(tauEff))/Foam::sqrt(2.));
   
    volScalarField preFac
    (
      (1. + eps_*tr(tauEff)/G_)
    * Foam::max(zeroD, (sigmaEff - k0_*Foam::pow(lambda_, m2_))/(2.*sigmaEff*eta0_*Foam::pow(lambda_,m1_)+small) )  
    );
   
    fvSymmTensorMatrix tauEqn
    (
         fvm::ddt(tau_)
       + fvm::div(phi(), tau_) 
       - twoSymm(C)
     ==
         G_*twoD 
       - fvm::Sp(2*G_*preFac, tau_)
       + 2.*G_*preFac*C0_*Foam::pow(lambda_,m3_)*A_
    );
 
    tauEqn.relax();
    tauEqn.solve();
    
    volScalarField magDvp = Foam::mag(preFac*(tau_-C0_*Foam::pow(lambda_,m3_)*A_))/Foam::sqrt(2.);
 
//- Solve for A
    volTensorField CA = A_ & L;
    
    fvSymmTensorMatrix AEqn
    (
         fvm::ddt(A_)
       + fvm::div(phi(), A_) 
       - twoSymm(CA)
     ==
         2.*preFac*tau_
       - fvm::Sp(2.*preFac*C0_*Foam::pow(lambda_,m3_), A_)
       - fvm::Sp(q_*Foam::sqrt(2.)*magDvp*(1.+eps_*tr(A_)), A_)          
    );
 
    AEqn.relax();
    AEqn.solve();
 

//- Solve for lambda
    fvScalarMatrix lambdaEqn 
    (
       fvm::ddt(lambda_)
     ==
       k1_+k2_*Foam::pow(2.*magDvp*oneDim, n1_)
     - fvm::Sp(k1_+k2_*Foam::pow(2.*magDvp*oneDim, n1_), lambda_)
     - k3_*Foam::pow(2.*magDvp*oneDim, n2_)*Foam::pow(lambda_,n3_)
    );
 
    lambdaEqn.relax();
    lambdaEqn.solve();
}


// ************************************************************************* //
