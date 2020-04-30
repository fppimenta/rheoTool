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

#include "ML_IKH.H" 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(ML_IKH, 0);
    addToRunTimeSelectionTable(constitutiveEq, ML_IKH, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::constitutiveEqs::ML_IKH::ML_IKH::lambda::lambda
(
   const word& name,
   const volVectorField& U,
   const surfaceScalarField& phi,
   const dictionary& dict
)
:  
    lambdaI_
    (
        IOobject
        (
            "lambda" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("1",dimless,1),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    )
{}

Foam::constitutiveEqs::ML_IKH::ML_IKH
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
    A_
    (
        IOobject
        (
            "A" + name,
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
    C_(dict.lookup("C")),  
    D_(dict.lookup("D")),  
    a_(dict.lookup("a")),
    n_(dict.lookup("n")),
    b_(dict.lookup("b")),
    ky_(dict.lookup("ky")),
    kh_(dict.lookup("kh")),
    k1_(dict.lookup("k1")),
    k2_(dict.lookup("k2")),
    k3_(dict.lookup("k3")),
    q_(dict.lookup("q")),
    lambdaE_(dict.lookup("lambdaE")),
    isSC_(dict.lookup("isStressControlled")), 
    dims_(dict.lookupOrDefault<vector>("dims", U.mesh().solutionD() )),
    nDims(scalar(dims_.x()+dims_.y()+dims_.z())),
    Itensor("Identity", dimless, symm(tensor::I)),   
    lambdas_()
{
    // Stabilization
    checkForStab(dict);
    
    // Adjust Itensor
    Itensor.value().xx() = dims_.x();
    Itensor.value().yy() = dims_.y();
    Itensor.value().zz() = dims_.z();
    
    // Build lambdas
    int nLambdas_ = C_.size();
    lambdas_.setSize(nLambdas_);

    std::stringstream int2str;
    forAll (lambdas_, i)
    {  
        int2str.str(std::string());
        int2str << i; 
      
        lambdas_.set
        (
            i,
            new lambda
            (
                name + int2str.str(),
                U,
                phi,
                dict
            )
        );  
    }  
}
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::ML_IKH::correct()
{
//- Aux numbers     
    dimensionedScalar ZScalar("0",dimless/dimTime,0.);
     
    dimensionedScalar small("S",dimPressure,1e-16);
    
    dimensionedScalar small1("small",dimensionSet(2, -2, -3, 0, 0, 0, 0),1e-16);
    
    dimensionedScalar onePa("1",dimPressure,1);
    
//- Assemble lambda
    volScalarField lambda = lambdas_[0].lambdaI() * C_[0];
    for (int i = 1; i<C_.size(); i++)   
      lambda += lambdas_[i].lambdaI() * C_[i];   

//- Solve for tau
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());
    
    // Convected derivate term
    volTensorField C = tau_ & L;

    // Twice the rate of deformation tensor
    volSymmTensorField twoD = twoSymm(L);
    
    // Inter vars 
    volSymmTensorField kback = kh_*symm(A_ + 2.*(A_&A_));
    volSymmTensorField tauEff = tau_ - kback;
    volScalarField sigmaBar = Foam::mag(tau_-Itensor*tr(tau_)/nDims)/Foam::sqrt(2.);
 
    volSymmTensorField Dp = tauEff*(sigmaBar-lambda*ky_)/(2.*lambda*etaP_*sigmaBar+small1);
   
    forAll(Dp, i)
    {
      if ( sigmaBar[i] < (lambda[i]*ky_.value()) )
       {
         Dp[i] *= 0.;
       }
    }
    
    volScalarField preFac =
    Foam::max
    ( 
      ZScalar,
     (1./(lambda*lambdaE_))*(sigmaBar - lambda*ky_)/(sigmaBar+small)
    );
    
    fvSymmTensorMatrix tauEqn
    (
         fvm::ddt(tau_)
       + fvm::div(phi(), tau_) 
       - twoSymm(C)
     ==
         (etaP_/lambdaE_)*twoD      
       - fvm::Sp(preFac, tau_)
       + preFac * kback
    );
 
    tauEqn.relax();
    tauEqn.solve();
    
//- Solve for A
    volTensorField W = (L - L.T())/2.; 
  
    fvSymmTensorMatrix AEqn
    (
         fvm::ddt(A_)
       + fvm::div(phi(), A_) 
       + symm( (W&A_) - (A_&W) )
     ==
         symm( (Dp & A_) + (A_ & Dp) )
       + Dp
       - fvm::Sp(q_*Foam::sqrt(1./2.)*Foam::mag(Dp), A_)       
    );
 
    AEqn.relax();
    AEqn.solve();
 
//- Solve for lambdaI
    volScalarField Phi(isSC_==true ? Foam::max(small*0., sigmaBar - lambda*ky_) : Foam::sqrt(2.)*Foam::mag(Dp));
    
    forAll(lambdas_, i)
    {
          volScalarField& lambdaI = lambdas_[i].lambdaI();
          
          fvScalarMatrix lambdaIEqn
          (
               fvm::ddt(lambdaI)
             ==
               D_[i]*
               (
                 - k1_*Foam::pow(Phi/onePa, a_)*Foam::pow(lambdaI, n_)
                 + k2_*Foam::pow(Phi/onePa, b_)*(1.-lambdaI)
                 + k3_                    
               )
             - fvm::Sp(k3_*D_[i], lambdaI)
          );
          
          lambdaIEqn.solve(phi().mesh().solverDict("lambdaI"));
    }
    
}


// ************************************************************************* //
