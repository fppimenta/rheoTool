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

#include "BMP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(BMP, 0);
    addToRunTimeSelectionTable(constitutiveEq, BMP, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::BMP::BMP
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
    Phi_
    (
        IOobject
        (
            "Phi" + name,
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
    lambda_(dict.lookup("lambda")),
    G0_(dict.lookup("G0")),
    k_(dict.lookup("k")),
    Phi0_(dict.lookup("Phi0")),
    PhiInf_(dict.lookup("PhiInf"))
{
 checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::BMP::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());
    
    // Fluidity transport equation  
    fvScalarMatrix PhiEqn
    (
         fvm::ddt(Phi_)
       + fvm::div(phi(), Phi_) 
     ==
       - fvm::Sp(1./lambda_, Phi_)
       + Phi0_/lambda_
       + k_ * (PhiInf_ - Phi_) * (tau_ && symm(L))  
    );
    
    PhiEqn.relax();
    PhiEqn.solve();

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
         fvm::ddt(tau_)
       + fvm::div(phi(), tau_) 
     ==
        G0_*twoD
      + twoSymm(C)
      - fvm::Sp(G0_*Phi_, tau_)
    );
 
    tauEqn.relax();
    tauEqn.solve();
 
}


// ************************************************************************* //
