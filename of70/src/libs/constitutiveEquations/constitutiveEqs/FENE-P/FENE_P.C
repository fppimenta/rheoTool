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
#include "FENE_P.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(FENE_P, 0);
    addToRunTimeSelectionTable(constitutiveEq, FENE_P, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::FENE_P::FENE_P
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
    L2_(dict.lookup("L2")),
    lambda_(dict.lookup("lambda")),
    solveInTau_(dict.lookupOrDefault<Switch>("solveInTau", false)),
    modifiedForm_(dict.lookupOrDefault<Switch>("modifiedForm", false)),
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
	dimensionedSymmTensor("I", L2_.dimensions(), symmTensor::I),
	tau_.boundaryField().types()
    ),
    varf_
    (
        IOobject
        (
            "varf_" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1./ ( (L2_ + tr(tau_)*lambda_/etaP_)/(L2_ - 3.) )
    ),
    thermoLambdaPtr_(thermoFunction::New("thermoLambda", U.mesh(), dict)),  
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict))  
{
 checkForStab(dict);
 
 checkIfCoupledSolver(U.mesh().solutionDict(), tau_); 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::FENE_P::correct()
{

// Update temperature-dependent properties
volScalarField lambda = thermoLambdaPtr_->createField(lambda_);
volScalarField etaP = thermoEtaPtr_->createField(etaP_);

dimensionedSymmTensor Ist
( 
 "Identity", 
 varf_.dimensions(), 
 symmTensor::I
);
    
// Velocity gradient tensor
volTensorField L = fvc::grad(U());

if (!solveInTau_)
{
    // Convected derivate term
    volTensorField C = A_ & L;
    
    // Update varf
    dimensionedScalar a = L2_/(L2_-3); 
    varf_ = (1./(1.-tr(A_)/L2_));

    // Stress transport equation
    fvSymmTensorMatrix AEqn
    (
        fvm::ddt(A_)
      + fvm::div(phi(), A_)
     ==
        twoSymm(C)
     -  fvm::Sp(varf_/lambda, A_)
     +  (a/lambda) * Ist
    );

    AEqn.relax();
    AEqn.solve();
    
    varf_ = (1./(1.-tr(A_)/L2_));
    
    tau_ = (etaP/lambda) * (varf_*A_ - a*Ist);
    
    tau_.correctBoundaryConditions(); 
}
else
{
    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);
    
    // Update f
    dimensionedScalar a = L2_/(L2_-3); 
    varf_ = 1./ ( (L2_ + tr(tau_)*lambda/(a*etaP))/(L2_ - 3.) );

    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
     ==
        twoSymm(C)
      - fvm::Sp(1/(lambda*varf_), tau_)
    );
    
    if(!modifiedForm_)
    {
       tauEqn += (1/varf_)*(tau_ + (a*etaP/lambda)*Ist)*( fvc::ddt(varf_) + fvc::div(phi(), varf_, "div(phi,1|f)") ); 
    }

    tauEqn.relax();
    
    if (!solveCoupled_)
    {
      solve(tauEqn == (etaP*a/lambda)*twoD);
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
        fvmb::twoSymmGrad(-(etaP*a/lambda), U())
      );         
    }  
}
 
}


// ************************************************************************* //
