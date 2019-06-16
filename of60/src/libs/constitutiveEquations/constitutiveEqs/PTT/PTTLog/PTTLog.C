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

#include "PTTLog.H"
#include "addToRunTimeSelectionTable.H"
 
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace constitutiveEqs 
 {
    defineTypeNameAndDebug(PTTLog, 0);
    addToRunTimeSelectionTable(constitutiveEq, PTTLog, dictionary);
    
    template<>
    const char* NamedEnum
    <
      PTTLog::PTTFunctions,
      3
    >::names[] =
    {
      "linear",
      "exponential",
      "generalized"
    };
  
    const NamedEnum<PTTLog::PTTFunctions, 3> PTTLog::PTTFunctionNames_;
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::PTTLog::PTTLog
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
    MLrtol_(dict.lookupOrDefault<scalar>("rTolMittagLeffler", 1e-12)),
    MLmaxIter_(dict.lookupOrDefault<int>("maxIterMittagLeffler", 200)),
    PTTFunction_(PTTFunctionNames_.read(dict.lookup("destructionFunctionType")))
{
 checkForStab(dict);
 
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

void Foam::constitutiveEqs::PTTLog::correct()
{
    // Decompose grad(U).T()

    volTensorField L = fvc::grad(U());

    dimensionedScalar c1( "zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.);
    volTensorField   B = c1 * eigVecs_; 
    volTensorField   omega = B;
    volTensorField   M = (eigVecs_.T() & ( L.T() - zeta_*symm(L) ) & eigVecs_);

    decomposeGradU(M, eigVals_, eigVecs_, omega, B);
    
    dimensionedTensor Itensor
    ( 
        "Identity", 
        dimensionSet(0, 0, 0, 0, 0, 0, 0), 
        tensor::I 
    );
        
    // Select function    
    volTensorField extFun = (1.0/lambda_) * ( eigVecs_ & (inv(eigVals_) - Itensor) & eigVecs_.T() ); 
      
    switch (PTTFunction_) 
    {    
      case pfLinear :      
        extFun *= ( 1. + (epsilon_/(1-zeta_))*( tr( (eigVecs_ & eigVals_ & eigVecs_.T()) ) - 3.) );
        break;
      
      case pfExpt : 
        extFun *= Foam::exp( (epsilon_/(1-zeta_))*( tr( (eigVecs_ & eigVals_ & eigVecs_.T()) ) - 3.) );
        break;
      
      case pfGen :
        scalar gammaBeta(gammaFunValues_[0]);
        volScalarField z = (epsilon_/(1-zeta_))*(tr((eigVecs_ & eigVals_ & eigVecs_.T())) - 3.);
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
            << "Try increasing maxIterMittagLeffler parameter in constitutiveProperties." << endl;
         }
         
         Eab[i] = sum;
        } 
         
        extFun *= gammaBeta*Eab;
        break;   
    }

  
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
        + extFun
       )       
    );
   
    thetaEqn.relax();
    thetaEqn.solve();
  
    // Diagonalization of theta

    calcEig(theta_, eigVals_, eigVecs_);

    // Convert from theta to tau

    tau_ = (etaP_/(lambda_*(1-zeta_))) * symm( (eigVecs_ & eigVals_ & eigVecs_.T()) - Itensor);

    tau_.correctBoundaryConditions();

}


// ************************************************************************* //
