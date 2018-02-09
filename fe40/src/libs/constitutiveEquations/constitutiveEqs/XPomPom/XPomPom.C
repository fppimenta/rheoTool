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

#include "XPomPom.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(XPomPom, 0);
    addToRunTimeSelectionTable(constitutiveEq, XPomPom, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XPomPom::XPomPom
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
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    lambdaS_(dict.lookup("lambdaS")),
    lambdaB_(dict.lookup("lambdaB")),
    alpha_(dict.lookup("alpha")),
    q_(dict.lookup("q"))
{
 checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::XPomPom::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);
    
    dimensionedTensor Itensor
    ( 
        "Identity", 
        dimensionSet(0, 0, 0, 0, 0, 0, 0), 
        tensor::I 
    );
    
    volScalarField lambda( Foam::sqrt( 1. + tr(tau_)/(3.*etaP_/lambdaB_) ) );
    
    volScalarField f
    ( 
        2.*(lambdaB_/lambdaS_)*Foam::exp( (2./q_)*(lambda-1.) ) * (1. - 1./lambda)
     + (1./(lambda*lambda)) * ( 1. - (alpha_/3.) * tr(tau_&tau_) / Foam::sqr(etaP_/lambdaB_) )
    );

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
         fvm::ddt(tau_)
       + fvm::div(phi(), tau_) 
     ==
        etaP_/lambdaB_*twoD
      + twoSymm(C)
      - fvm::Sp(f/lambdaB_, tau_)
      - symm
        (
           (alpha_ / etaP_) * (tau_&tau_)
         + (etaP_/(lambdaB_*lambdaB_)) * (f - 1.) * Itensor
        )
    );
 
    tauEqn.relax();
    tauEqn.solve();
 
}


// ************************************************************************* //
