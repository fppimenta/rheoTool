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

#include "RoliePoly.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RoliePoly, 0);
    addToRunTimeSelectionTable(constitutiveEq, RoliePoly, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RoliePoly::RoliePoly
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
    lambdaR_(dict.lookup("lambdaR")),
    lambdaD_(dict.lookup("lambdaD")),
    beta_(dict.lookup("beta")),
    delta_(dict.lookup("delta")),
    chiMax_(dict.lookup("chiMax")),
    solveInTau_(dict.lookupOrDefault<Switch>("solveInTau", false)),
    A_
    (
        IOobject
        (
            "A" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            solveInTau_ == true ? (IOobject::NO_WRITE) : (IOobject::AUTO_WRITE)
        ),
        U.mesh(),
	dimensionedSymmTensor("I", dimless, symmTensor::I),
	tau_.boundaryField().types()
    )
{
 checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::RoliePoly::correct()
{

   // Velocity gradient tensor
    volTensorField L = fvc::grad(U());

if (!solveInTau_)
 {
   
    dimensionedSymmTensor Ist
    ( 
        "Identity", 
        A_.dimensions(), 
        symmTensor::I
    );
    
    // Convected derivate term
    volTensorField C = A_ & L;
    
    volScalarField trA(tr(A_));
    
    volScalarField M1
     ( 
       2.*(1.-Foam::sqrt(3./trA))/lambdaR_
     );
     
    if (chiMax_.value() > 1.)
     {
    
       M1 *= 
            (
               (  3.-(trA/3.)/sqr(chiMax_) ) * (1.-1./sqr(chiMax_))
             / ( (1.-(trA/3.)/sqr(chiMax_) ) * (3.-1./sqr(chiMax_)  ) ) 
            ); 
     }
     
    // Stress transport equation
    fvSymmTensorMatrix AEqn
    (
        fvm::ddt(A_)
      + fvm::div(phi(), A_)
     ==
        twoSymm(C)
      - fvm::Sp(1./lambdaD_, A_)
      + (1./lambdaD_)*Ist  
      - M1 * ( A_ + beta_*pow(trA/3., delta_)*(A_-Ist) )    
    );

    AEqn.relax();
    AEqn.solve();
    
    tau_ = (etaP_/lambdaD_) * (A_ - Ist);
    
    if (chiMax_.value() > 1.)
     {
     
       trA = tr(A_); //Update
        
       tau_ *= 
            (
               (  3.-(trA/3.)/sqr(chiMax_) ) * (1.-1./sqr(chiMax_))
             / ( (1.-(trA/3.)/sqr(chiMax_) ) * (3.-1./sqr(chiMax_)  ) ) 
            );
     }
    
    tau_.boundaryField().updateCoeffs();
    tau_.correctBoundaryConditions();
 }
else
 {
      // TODO: To be implemented 
      
       FatalErrorIn
        (
            "Foam::RoliePoly::correct()"
        )   << "Solving in tau not implemented!"
            << abort(FatalError);
 } 
 
}


// ************************************************************************* //
