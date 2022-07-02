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

#include "RRM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(RRM, 0);
    addToRunTimeSelectionTable(constitutiveEq, RRM, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::RRM::RRM
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
    S_
    (
        IOobject
        (
            "S" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    Lstar_
    (
        IOobject
        (
            "Lstar" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    lambdaS_(dict.lookup("lambdaS")),
    G0_(dict.lookup("G0")),
    Dr0_(dict.lookup("Dr0")),
    m_(dict.lookup("m")),
    alpha_(dict.lookup("alpha")),
    beta_(dict.lookup("beta")),
    k_(dict.lookup("k"))
{
  checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
void Foam::constitutiveEqs::RRM::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
    // Evolve Lstar    
    volScalarField Rsby = lambdaS_ / (1. - Foam::sqr( Lstar_/(alpha_ + beta_/ ( strainRate(gradU)/Dr0_ + 1e-16 ) ) ) ); 
  
    fvScalarMatrix LstarEqn
    (
        fvm::ddt(Lstar_)
     ==
        Rsby*Dr0_
      - fvm::Sp(Rsby*Dr0_, Lstar_) 
      + k_*Dr0_*Foam::sqrt(1.5)*mag(dev(S_))  
         
    );
     
    LstarEqn.relax();
    LstarEqn.solve();
    
    // Dr
    volScalarField Dr = Dr0_*Foam::pow(Lstar_, -3) * (1. + Foam::log(Lstar_)/m_); 

    // Velocity gradient tensor  
    volTensorField K( gradU == nullptr ? fvc::grad(U())() : *gradU );
    
    // Convected derivate term
    volTensorField C = S_ & K; 
 
    // Rate of strain
    volSymmTensorField D = symm(K);
 
    // Orientation tensor transport equation
    fvSymmTensorMatrix SEqn
    (
         fvm::ddt(S_)
     ==
         2.*Dr*symm(tensor::I)
       - fvm::Sp(6.*Dr, S_)
       + twoSymm(C)
       - 2.*KDDFourthMoment(S_, D) 
         
    );
 
    SEqn.relax();
    SEqn.solve();
    
    // Compute tau from S (update KDDFourthMoment; may use old to save time)
    tau_ = (G0_/Lstar_)*( 3.*dev(S_) + (0.5/Dr)*KDDFourthMoment(S_, D) );  
    
    tau_.correctBoundaryConditions();
}


// ************************************************************************* //
