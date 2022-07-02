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

#include "Saramito.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(Saramito, 0);
    addToRunTimeSelectionTable(constitutiveEq, Saramito, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::Saramito::Saramito
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
    tau0_(dict.lookup("tau0")),
    n_(dict.lookup("n")),
    k_( n_.value() == 1 ? etaP_*1. : dict.lookup("k")),
    dims_(dict.lookupOrDefault<vector>("dims", U.mesh().solutionD() )),
    nDims(scalar(dims_.x()+dims_.y()+dims_.z())),
    Itensor("Identity", dimless, symm(tensor::I)),
    writeII_(dict.lookupOrDefault<Switch>("writeSecondInvariantTauDev", false)),
    funcPTT(0)    
{
    // Stabilization 
    checkForStab(dict);
 
    // Adjust Itensor
    Itensor.value().xx() = dims_.x();
    Itensor.value().yy() = dims_.y();
    Itensor.value().zz() = dims_.z();
  
    // Check PTT function
    if (n_.value() == 1.)
    {
      word pttfun(dict.lookup("PTTfunction"));
  
      if (pttfun == "none")
      {
        funcPTT = 0;
        Info << "PTT function: none.\n";
      }
      else if (pttfun == "linear")
      {
        funcPTT = 1;
        Info << "PTT function: linear.\n";
      }
      else if (pttfun == "exponential")
      {
        funcPTT = 2;
        Info << "PTT function: exponential.\n";
      }
      else
      {
        FatalErrorIn("Foam::constitutiveEqs::Saramito::Saramito\n")
            << "\nThe PTT function specified does not exist.\n"
            << "\nAvailable PTT functions are:\n"
            << "\n. none" <<"\n. linear" << "\n. exponential" 
            << abort(FatalError);
      }   
   
   }
   else
   {
      // if n != 1 PTT regularization is not allowed
      funcPTT = 0;
      Info << "PTT function: none.\n";
   }
} 

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::Saramito::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
    // Velocity gradient tensor
    volTensorField L( gradU == nullptr ? fvc::grad(U())() : *gradU );

    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);
      
    // 2nd invariant of deviatoric stress  
    volScalarField tauDMag = Foam::mag(tau_-Itensor*tr(tau_)/nDims)/Foam::sqrt(2.);
    if (writeII_ && tau_.time().outputTime())
     {
       tauDMag.rename("II_tauD");
       tauDMag.write();
     }
    
    // Pre-factor 
    dimensionedScalar ZScalar("Zero", dimless/(dimPressure*dimTime), 0.);
        
    dimensionedScalar small("small", dimPressure*dimPressure*dimTime, 1e-16);
    
    dimensionedScalar oneU("1", Foam::pow(dimPressure, 1.-n_.value()), 1.);
    
    dimensionedScalar oneK("1", Foam::pow(dimless/(dimPressure*dimTime), 1.-1./n_.value()), 1.);
    
    // If n = 1 there is no need to waste CPU time in two pow() operations
    volScalarField fac
    ( 
      n_.value() == 1.
      ?
      Foam::max( ZScalar, (tauDMag - tau0_)/(k_*tauDMag + small) )
      :
      Foam::pow
      (
        Foam::max(ZScalar, (tauDMag - tau0_)/(k_*Foam::pow(tauDMag, n_.value())*oneU + small) ),
        1./n_.value()
      )*oneK    
    ); 
    
    // Stress transport equation
    fvSymmTensorMatrix tauEqn
    (
        fvm::ddt(tau_)
      + fvm::div(phi(), tau_)
      - twoSymm(C)
     ==
        etaP_/lambda_*twoD   
    );
     
    // Add remaining terms of Gordon derivative if zeta is not zero
    if (zeta_.value() != 0.)
    {
      tauEqn += 0.5*zeta_*(symm(tau_ & twoD) + symm(twoD & tau_));
    }
    
    // No PTT term
    if (funcPTT==0)
    {
       tauEqn += fvm::Sp(fac*etaP_/lambda_, tau_); 
    }
    // Add linear PTT term
    else if (funcPTT==1)
    {
       tauEqn += fvm::Sp(fac*(epsilon_*tr(tau_) + etaP_/lambda_), tau_); 
    }
    // Add exp PTT term
    else if (funcPTT==2)
    {
       tauEqn += fvm::Sp(fac*( (etaP_/lambda_)*Foam::exp(epsilon_*lambda_/etaP_*tr(tau_)) ), tau_); 
    }

    tauEqn.relax();
    tauEqn.solve();
 
}


// ************************************************************************* //
