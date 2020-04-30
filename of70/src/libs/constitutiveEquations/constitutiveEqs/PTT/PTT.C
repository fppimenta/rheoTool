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
#include "PTT.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace constitutiveEqs 
 {
    defineTypeNameAndDebug(PTT, 0);
    addToRunTimeSelectionTable(constitutiveEq, PTT, dictionary);
    
    template<>
    const char* NamedEnum
    <
      PTT::PTTFunctions,
      3
    >::names[] =
    {
      "linear",
      "exponential",
      "generalized"
    };
  
    const NamedEnum<PTT::PTTFunctions, 3> PTT::PTTFunctionNames_;
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::PTT::PTT
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
    epsilon_(dict.lookup("epsilon")),
    zeta_(dict.lookup("zeta")),
    lambda_(dict.lookup("lambda")),
    MLrtol_(dict.lookupOrDefault<scalar>("rTolMittagLeffler", 1e-12)),
    MLmaxIter_(dict.lookupOrDefault<int>("maxIterMittagLeffler", 200)),
    thermoLambdaPtr_(thermoFunction::New("thermoLambda", U.mesh(), dict)),  
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict)),
    PTTFunction_(PTTFunctionNames_.read(dict.lookup("destructionFunctionType")))   
{
 checkForStab(dict);
 
 checkIfCoupledSolver(U.mesh().solutionDict(), tau_); 
 
 if (PTTFunction_ == pfGen)
 {
    alpha_ = dict.lookup("alpha");
    beta_ = dict.lookup("beta");
  
    // Check alpha and beta consistency
    if (alpha_.value()<=0 || beta_.value()<=0)
    {
      FatalErrorInFunction << "Both alpha and beta should be positive values "
      << "for the Mittag-Leffler function to converge." << endl << abort(FatalError);
    }

    // First value of field is gamma(beta)
    gammaFunValues_.append(tgamma(beta_.value()));
    
    // ...remaining values are for the denominator of Eab
    int k(0);   
    while (k < MLmaxIter_ && gammaFunValues_.last()<1e+100)  // last condition is to avoid overflow
    {
      gammaFunValues_.append(tgamma(alpha_.value()*k+beta_.value()));     
      k++;
    }
    
    // It might happen that the previous while exits in the condition of overflow (last one), 
    // in which case we need to redefine MLmaxIter_ to the max iterations performed
    // inside the while, otherwise later Eab would use undefined values of gammaFunValues_. 
    MLmaxIter_ = k;
 }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::PTT::correct()
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
      - 0.5*zeta_*(symm(tau_ & twoD) + symm(twoD & tau_))
    );
    
    
    switch (PTTFunction_) 
    {    
      case pfLinear :      
        tauEqn += fvm::Sp(epsilon_/etaP*tr(tau_) + 1./lambda, tau_);
        break;
      
      case pfExpt : 
        tauEqn += fvm::Sp( (1./lambda)*Foam::exp(epsilon_*lambda/etaP*tr(tau_)), tau_);
        break;
      
      case pfGen :
        scalar gammaBeta(gammaFunValues_[0]);
        volScalarField z = epsilon_*lambda*tr(tau_)/etaP;
        volScalarField Eab(z);        
        forAll(z, i)
        {
         scalar zi = z[i];
         scalar sum(0.0);
         scalar sumOld(0.0);
         scalar error(1.0);
         int k(0);

         while (k < MLmaxIter_ && error > MLrtol_)  
         {
            scalar Eabk = Foam::pow(zi, k)/gammaFunValues_[k+1];
            sumOld = sum;
            sum += Eabk;
            error = Foam::mag((sumOld-sum)/(sumOld+1e-12));            
            k++;
         }
         
         // Verification for warning
         if (k == MLmaxIter_)    
         {
           WarningInFunction
            << "Computation of the Mittag-Leffler function does not converged." << nl
            << "Iterations: " << k << ", relative error: " << error << "." << nl
            << "Try increasing maxIterMittagLeffler parameter in constitutiveProperties." << nl << endl;
         }
         
         Eab[i] = sum;
        } 
         
        tauEqn += fvm::Sp( (1./lambda)*gammaBeta*Eab, tau_); 
        break;   
    }
    
    tauEqn.relax();
    
    
    if (!solveCoupled_)
    {
      solve(tauEqn == etaP/lambda*twoD);
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
        fvmb::twoSymmGrad(-etaP/lambda, U())
      );   
    }
 
}


// ************************************************************************* //
