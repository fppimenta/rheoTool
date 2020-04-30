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

#include "coupledSolver.H" 
#include "blockOperators.H" 
#include "Oldroyd_B.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(Oldroyd_B, 0);
    addToRunTimeSelectionTable(constitutiveEq, Oldroyd_B, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::Oldroyd_B::Oldroyd_B
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
    lambda_(dict.lookup("lambda")),
    thermoLambdaPtr_(thermoFunction::New("thermoLambda", U.mesh(), dict)),  
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))  
{
 checkForStab(dict);
 
 checkIfCoupledSolver(U.mesh().solutionDict(), tau_); 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
void Foam::constitutiveEqs::Oldroyd_B::correct()
{
 // Update temperature-dependent properties
 volScalarField lambda = thermoLambdaPtr_->createField(lambda_);
 volScalarField etaP = thermoEtaPtr_->createField(etaP_);

 // Velocity gradient tensor
 volTensorField L = fvc::grad(U());

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
    twoSymm(C)
  - fvm::Sp(1/lambda, tau_)
 );
 
 tauEqn.relax();
    
 if (!solveCoupled_)
 {
   solve(tauEqn == etaP/lambda*twoD);
 }
 else
 {   
   // Get the solver
   coupledSolver& cps = U().time().lookupObjectRef<coupledSolver>( word("Uptau." + U().mesh().name()) ); 
      
   // Insert tauEqn
   cps.insertEquation
   (
     tau_.name(),
     tau_.name(),
     tauEqn
   );

   // Insert term (gradU + gradU.T)
   cps.insertEquation
   (
     tau_.name(),
     U().name(),
     fvmb::twoSymmGrad(-etaP/lambda, U())
   );   
 } 
}

// ************************************************************************* //
