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
#include "XPomPom.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(XPomPom, 0);
    addToRunTimeSelectionTable(constitutiveEq, XPomPom, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::XPomPom::XPomPom
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
 q_(dict.lookup("q")),
 n_(dict.lookup("n")),
 thermoLambdaBPtr_(thermoFunction::New("thermoLambdaB", U.mesh(), dict)), 
 thermoLambdaSPtr_(thermoFunction::New("thermoLambdaS", U.mesh(), dict)),  
 thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))  
{
 checkForStab(dict);
 
 checkIfCoupledSolver(U.mesh().solutionDict(), tau_); 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::XPomPom::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
  // Update temperature-dependent properties
  volScalarField lambdaB = thermoLambdaBPtr_->createField(lambdaB_);
  volScalarField lambdaS = thermoLambdaSPtr_->createField(lambdaS_);
  volScalarField etaP = thermoEtaPtr_->createField(etaP_);

  // Velocity gradient tensor
  volTensorField L( gradU == nullptr ? fvc::grad(U())() : *gradU );

  // Convected derivate term
  volTensorField C = tau_ & L;

  // Twice the rate of deformation tensor
  volSymmTensorField twoD = twoSymm(L);
    
  dimensionedTensor Itensor("Identity", dimless, tensor::I);
   
  volScalarField lambda( Foam::sqrt( 1. + tr(tau_)/(3.*etaP/lambdaB) ) );
    
  volScalarField f
  ( 
   n_.value() == 0
   ?
        2.*(lambdaB/lambdaS)*Foam::exp( (2./q_)*(lambda-1.) ) * (1. - 1./lambda)
     + (1./(lambda*lambda)) * ( 1. - (alpha_/3.) * tr(tau_&tau_) / Foam::sqr(etaP/lambdaB) )
   :
        2.*(lambdaB/lambdaS)*Foam::exp( (2./q_)*(lambda-1.) ) * (1. - 1./Foam::pow(lambda, n_+1))
     + (1./(lambda*lambda)) * ( 1. - (alpha_/3.) * tr(tau_&tau_) / Foam::sqr(etaP/lambdaB) )
  );
  

  // Stress transport equation
  fvSymmTensorMatrix tauEqn
  (
       fvm::ddt(tau_)
     + fvm::div(phi(), tau_) 
   ==
      twoSymm(C)
    - fvm::Sp(f/lambdaB, tau_)
    - symm
      (
         (alpha_ / etaP) * (tau_&tau_)
       + (etaP/(lambdaB*lambdaB)) * (f - 1.) * Itensor
      )
  );
 
  tauEqn.relax();
  
  if (!solveCoupled_)
  {
    solve(tauEqn == etaP/lambdaB*twoD);
  }
  else
  {   
    // Get the solver
    coupledSolver& cps = U().time().lookupObjectRef<coupledSolver>(word("Uptau."+U().mesh().name())); 
    
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
     fvmb::twoSymmGrad(-etaP/lambdaB, U())
   );   
 } 
 
}


// ************************************************************************* //
