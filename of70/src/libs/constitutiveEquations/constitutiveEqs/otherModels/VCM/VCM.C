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

#include "VCM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(VCM, 0);
    addToRunTimeSelectionTable(constitutiveEq, VCM, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::VCM::VCM
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
    nA_
    (
        IOobject
        (
            "nA" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    nB_
    (
        IOobject
        (
            "nB" + name,
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
    B_
    (
        IOobject
        (
            "B" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),    
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookupOrDefault<dimensionedScalar>("etaP",dimensionedScalar("0", etaS_.dimensions(), 0.))),
    lambdaA_(dict.lookup("lambdaA")),
    DA_(dict.lookup("DA")),
    DB_(dict.lookup("DB")),
    chi_(dict.lookup("chi")),
    cAEq_(dict.lookup("cAEq")),
    cBEq_(dict.lookup("cBEq")),
    G0_(dict.lookup("G0")),
    eps_(dict.lookup("epsilon"))
{
    checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::VCM::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
 // Velocity gradient tensor and gammaDot
    volTensorField L( gradU == nullptr ? fvc::grad(U())() : *gradU );   
    volSymmTensorField gDot = twoSymm(L);
  
 // Breakage rate
    volScalarField cA = cAEq_ + (1./3.)*chi_*(gDot&&(A_/nA_));
     
 // Solve for nA
  
    fvScalarMatrix nAEqn
     (
         fvm::ddt(nA_) 
       + fvm::div(phi(), nA_)
      == 
         fvm::laplacian(2.*DA_,nA_)
       + (0.5/lambdaA_)*cBEq_*nB_*nB_
       - fvm::SuSp(cA/lambdaA_, nA_)  
     //  - fvc::div(fvc::div(A_*DA_, "divLinear"), "divLinear")
     );
  
    nAEqn.relax();   
    nAEqn.solve();

 // Solve for nB  
    
    fvScalarMatrix nBEqn
     (
         fvm::ddt(nB_) 
       + fvm::div(phi(), nB_)
      == 
         fvm::laplacian(2*DB_,nB_)
       - fvm::Sp(cBEq_*nB_/lambdaA_, nB_)
       + (2./lambdaA_)*cA*nA_ 
    //   - fvc::div(fvc::div(2.*B_*DB_, "divLinear"), "divLinear")
     );
  
    nBEqn.relax();   
    nBEqn.solve();


 // Update breakage rate
    cA = cAEq_ + (1./3.)*chi_*(gDot&&(A_/nA_));

 // Solve for A
    
    // Convected derivate term
    volTensorField CA = A_ & L;
    
    fvSymmTensorMatrix AEqn
    (
         fvm::ddt(A_)
       + fvm::div(phi(), A_)
       - twoSymm(CA)
       + fvm::Sp(1./lambdaA_, A_)
       - (nA_/lambdaA_) * symm(tensor::I)
       - fvm::laplacian(DA_, A_)
     ==
         cBEq_*nB_*B_/lambdaA_
       - fvm::SuSp(cA/lambdaA_, A_)  
    );
 
    AEqn.relax();
    AEqn.solve();
    
 // Solve for B
  
    // Convected derivate term
    volTensorField CB = B_ & L;
    
    // Compute lambdaB from given parameters
    dimensionedScalar lambdaB = lambdaA_*eps_;
    
    fvSymmTensorMatrix BEqn
    (
         fvm::ddt(B_)
       + fvm::div(phi(), B_) 
       - twoSymm(CB)     
       + fvm::Sp(1./lambdaB, B_)
       - 0.5*(nB_/lambdaB) * symm(tensor::I)
       - fvm::laplacian(DB_, B_)
     ==
       - fvm::Sp(2.*cBEq_*nB_/lambdaA_, B_)
       + 2.*cA*A_/lambdaA_  
    );
 
    BEqn.relax();
    BEqn.solve();
    
    tau_ = G0_*( A_ + 2.*B_ - (nA_ + nB_)* symm(tensor::I) );
    tau_.correctBoundaryConditions(); 
}


// ************************************************************************* //
