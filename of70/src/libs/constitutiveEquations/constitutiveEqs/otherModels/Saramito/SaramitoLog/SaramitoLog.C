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

#include "SaramitoLog.H"
#include "addToRunTimeSelectionTable.H"
 
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(SaramitoLog, 0);
    addToRunTimeSelectionTable(constitutiveEq, SaramitoLog, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::SaramitoLog::SaramitoLog
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
    theta_
    (
        IOobject
        (
            "theta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    eigVals_
    (
        IOobject
        (
            "eigVals" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
        extrapolatedCalculatedFvPatchField<tensor>::typeName
    ),
    eigVecs_
    (
        IOobject
        (
            "eigVecs" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
        extrapolatedCalculatedFvPatchField<tensor>::typeName
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
    dims_(dict.lookup("dims")),
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

void Foam::constitutiveEqs::SaramitoLog::correct
(
  const volScalarField* alpha,
  const volTensorField* gradU
)
{
   // Decompose grad(U).T()

    volTensorField L( gradU == nullptr ? fvc::grad(U())() : *gradU );
 
    dimensionedScalar c1( "zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.);
    volTensorField   B = c1 * eigVecs_; 
    volTensorField   omega = B;
    volTensorField   M = (eigVecs_.T() & ( L.T() - zeta_*symm(L) ) & eigVecs_);

    decomposeGradU(M, eigVals_, eigVecs_, omega, B);
    
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
 
    // Solve the constitutive Eq in theta = log(c)  
    fvSymmTensorMatrix thetaEqn
    (
         fvm::ddt(theta_)
       + fvm::div(phi(), theta_)
       ==
       symm
       (  
         (omega&theta_)
       - (theta_&omega)
       + 2.0 * B
       )
    );
    
    // No PTT term
    if (funcPTT==0)
    {
       thetaEqn -= 
       symm
       (
          (fac*etaP_/lambda_)
        * (
             eigVecs_ &
               (
                  inv(eigVals_)  
                - dimensionedTensor("I",dimless,tensor::I)
               )
           & eigVecs_.T() 
          )
       ) ;
    }
    // Add linear PTT term
    else if (funcPTT==1)
    {
       thetaEqn -= 
       symm
       (
          (fac*etaP_/lambda_) * ( 1. + (epsilon_/(1-zeta_))*( tr( (eigVecs_ & eigVals_ & eigVecs_.T()) ) - 3.) )
        * (
             eigVecs_ &
               (
                  inv(eigVals_)  
                - dimensionedTensor("I",dimless,tensor::I)
               )
           & eigVecs_.T() 
          )
       ) ;
    }
    // Add exp PTT term
    else if (funcPTT==2)
    {
       thetaEqn -=  
       symm 
       (
          (fac*etaP_/lambda_) * Foam::exp( (epsilon_/(1-zeta_))*( tr( (eigVecs_ & eigVals_ & eigVecs_.T()) ) - 3.) )
        * (
             eigVecs_ &
               (
                  inv(eigVals_)  
                - tensor::I
               )
           & eigVecs_.T() 
          )
       ) ;
    }
   
    thetaEqn.relax();
    thetaEqn.solve();
  
    // Diagonalization of theta

    calcEig(theta_, eigVals_, eigVecs_);

    // Convert from theta to tau

    tau_ = (etaP_/(lambda_*(1-zeta_))) * symm( (eigVecs_ & eigVals_ & eigVecs_.T()) - tensor::I);

    tau_.correctBoundaryConditions();

}


// ************************************************************************* //
