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

#include "constitutiveEq.H"
#include "coupledSolver.H" 
#include "blockOperators.H" 
#include <Eigen/Dense> // For eigen decomposition
#include "jacobi.H"    // Only required for jacobi decomposition

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(constitutiveEq, 0);
  defineRunTimeSelectionTable(constitutiveEq, dictionary);
  
  template<>
  const char* NamedEnum
  <
    constitutiveEq::stabOptions,
    3
  >::names[] =
  {
    "none",
    "BSD",
    "coupling"
  };
  
  const NamedEnum<constitutiveEq::stabOptions, 3> constitutiveEq::stabOptionNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constitutiveEq::constitutiveEq
(
  const word& name,
  const volVectorField& U,
  const surfaceScalarField& phi
)
:
  name_(name),
  U_(U),
  phi_(phi),
  solveCoupled_(false)
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<fvVectorMatrix> constitutiveEq::divTau
(   
   const volVectorField& U  
) const
{
   
 if (isGNF())
 {   
   return
   (
      fvm::laplacian( eta()/rho(), U, "laplacian(eta,U)")
    + (fvc::grad(U) & fvc::grad(eta()/rho()))
   ); 
 }
 else
 {  
   switch (stabOption_) 
   {    
     case soNone : // none       
     return
     (
       solveCoupled_
       ?
          fvm::laplacian(etaS()/rho(), U, "laplacian(eta,U)")
       :
          fvc::div(tau()/rho(), "div(tau)")
        + fvm::laplacian(etaS()/rho(), U, "laplacian(eta,U)")
     );
   
     case soBSD : // BSD   
     return
     (
       solveCoupled_
       ?
        - fvc::laplacian(etaP()/rho(), U, "laplacian(eta,U)")
        + fvm::laplacian( (etaP()+ etaS())/rho(), U, "laplacian(eta,U)")
       :
          fvc::div(tau()/rho(), "div(tau)")
        - fvc::laplacian(etaP()/rho(), U, "laplacian(eta,U)")
        + fvm::laplacian( (etaP()+ etaS())/rho(), U, "laplacian(eta,U)")
     );
  
     case soCoupling : // coupling       
     return
     (
       solveCoupled_
       ?
        - fvc::div((etaP()/rho())*fvc::grad(U),"div(grad(U))")
        + fvm::laplacian( (etaP() + etaS())/rho(), U, "laplacian(eta,U)")
       
       :
          fvc::div(tau()/rho(), "div(tau)")
        - fvc::div((etaP()/rho())*fvc::grad(U),"div(grad(U))")
        + fvm::laplacian( (etaP() + etaS())/rho(), U, "laplacian(eta,U)")
     );
     
     default: // This will never happen
     return tmp<fvVectorMatrix>(); 
   }      
 }    
}

tmp<fvVectorMatrix> constitutiveEq::divTauS
(   
   const volVectorField& U,  
   const volScalarField& alpha
) const
{       
 if (isGNF())
 {       
   return
   (
      fvm::laplacian( eta()*alpha, U, "laplacian(eta,U)")
    + fvc::div(eta()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
   ); 
 }
 else
 { 
   switch (stabOption_) 
   {    
     case soNone : // none    
     return
     (        
        fvm::laplacian(etaS()*alpha, U, "laplacian(eta,U)")
      + fvc::div(etaS()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
     );
   
     case soBSD : // BSD 
     return
     (
       - fvc::laplacian(etaP()*alpha, U, "laplacian(eta,U)")
       + fvm::laplacian( (etaP()+ etaS())*alpha, U, "laplacian(eta,U)")
       + fvc::div(etaS()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
     );
      
     case soCoupling : // coupling 	      
     return
     (        
       - fvc::div( (etaP()*alpha) * fvc::grad(U), "div(grad(U))")
       + fvm::laplacian( (etaP() + etaS())*alpha, U, "laplacian(eta,U)")
       + fvc::div(etaS()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
     ); 
     
     default: // This will never happen
     return tmp<fvVectorMatrix>();  
   }      
 }      
}

tmp<fvVectorMatrix> constitutiveEq::divTauThermo
(   
   const volVectorField& U  
) const
{
   
 if (isGNF())
 {   
   return
   (
      fvm::laplacian( eta()/rho(), U, "laplacian(eta,U)")
    + fvc::div((eta()/rho())*dev2(T(fvc::grad(U))), "div(eta*dev2(T(gradU)))")
   ); 
 }
 else
 {  
   switch (stabOption_) 
   {    
     case soNone : // none       
     return
     ( 
       solveCoupled_
       ?
          fvm::laplacian(etaSThermo()/rho(), U, "laplacian(eta,U)")
        + fvc::div((etaSThermo()/rho())*dev2(T(fvc::grad(U))), "div(eta*dev2(T(gradU)))")
       :
          fvc::div(tau()/rho(), "div(tau)")
        + fvm::laplacian(etaSThermo()/rho(), U, "laplacian(eta,U)")
        + fvc::div((etaSThermo()/rho())*dev2(T(fvc::grad(U))), "div(eta*dev2(T(gradU)))")
     );
   
     case soBSD : // BSD   
     return
     (
       solveCoupled_
       ?
        - fvc::laplacian(etaPThermo()/rho(), U, "laplacian(eta,U)")
        + fvm::laplacian( (etaPThermo()+ etaSThermo())/rho(), U, "laplacian(eta,U)")
        + fvc::div((etaSThermo()/rho())*dev2(T(fvc::grad(U))), "div(eta*dev2(T(gradU)))")
       :
          fvc::div(tau()/rho(), "div(tau)")
        - fvc::laplacian(etaPThermo()/rho(), U, "laplacian(eta,U)")
        + fvm::laplacian( (etaPThermo()+ etaSThermo())/rho(), U, "laplacian(eta,U)")
        + fvc::div((etaSThermo()/rho())*dev2(T(fvc::grad(U))), "div(eta*dev2(T(gradU)))")
     );
  
     case soCoupling : // coupling             
     return
     (
       solveCoupled_
       ?
        - fvc::div((etaPThermo()/rho())*fvc::grad(U), "div(grad(U))")
        + fvm::laplacian( (etaPThermo() + etaSThermo())/rho(), U, "laplacian(eta,U)")
        + fvc::div((etaSThermo()/rho())*dev2(T(fvc::grad(U))), "div(eta*dev2(T(gradU)))")
       :
          fvc::div(tau()/rho(), "div(tau)")
        - fvc::div((etaPThermo()/rho())*fvc::grad(U), "div(grad(U))")
        + fvm::laplacian( (etaPThermo() + etaSThermo())/rho(), U, "laplacian(eta,U)")
        + fvc::div((etaSThermo()/rho())*dev2(T(fvc::grad(U))), "div(eta*dev2(T(gradU)))")
       
     );
     
     default: // This will never happen
     return tmp<fvVectorMatrix>(); 
   }      
 }    
}

tmp<fvVectorMatrix> constitutiveEq::divTauSThermo
(   
   const volVectorField& U,  
   const volScalarField& alpha
) const
{       
 if (isGNF())
 {       
   return
   (
      fvm::laplacian(eta()*alpha, U, "laplacian(eta,U)")
    + fvc::div(eta()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
   ); 
 }
 else
 { 
   switch (stabOption_) 
   {    
     case soNone : // none    
     return
     (        
        fvm::laplacian(etaSThermo()*alpha, U, "laplacian(eta,U)")
      + fvc::div(etaSThermo()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
     );
   
     case soBSD : // BSD 
     return
     (
       - fvc::laplacian(etaPThermo()*alpha, U, "laplacian(eta,U)")
       + fvm::laplacian( (etaPThermo()+ etaSThermo())*alpha, U, "laplacian(eta,U)")
       + fvc::div(etaSThermo()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
     );
      
     case soCoupling : // coupling 	      
     return
     (        
       - fvc::div( (etaPThermo()*alpha) * fvc::grad(U), "div(grad(U))")
       + fvm::laplacian( (etaPThermo() + etaSThermo())*alpha, U, "laplacian(eta,U)")
       + fvc::div(etaSThermo()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
     ); 
     
     default: // This will never happen
     return tmp<fvVectorMatrix>();  
   }      
 }      
}

void constitutiveEq::divTauImplCoupled
(
) const
{
 if (solveCoupled_ && !isGNF())
 {
  coupledSolver& cps = const_cast<coupledSolver&>
  (
    U().time().lookupObject<coupledSolver>(word("Uptau."+U().mesh().name()))
  );
   
  cps.insertEquation
  (
    U().name(),
    tau()().name(),
    fvmb::div(-1./rho().value(), tau())
  );
 }
}

void constitutiveEq::decomposeGradU
(
  const volTensorField& M,
  const volTensorField& eigVals, 
  const volTensorField& eigVecs,
  volTensorField& omega, 
  volTensorField& B
)
{

 forAll(M, cellI)
 {
   const tensor& eigValsR = eigVals[cellI];
   const tensor& MR = M[cellI];
   tensor& omegaR = omega[cellI];
        
   B[cellI].xx()=MR.xx();
   B[cellI].yy()=MR.yy(); 
   B[cellI].zz()=MR.zz(); 
         
   omegaR.xy() = ( eigValsR.yy()*MR.xy() + eigValsR.xx()*MR.yx() ) / ( eigValsR.yy() - eigValsR.xx() + 1e-16);
   omegaR.xz() = ( eigValsR.zz()*MR.xz() + eigValsR.xx()*MR.zx() ) / ( eigValsR.zz() - eigValsR.xx() + 1e-16);
   omegaR.yz() = ( eigValsR.zz()*MR.yz() + eigValsR.yy()*MR.zy() ) / ( eigValsR.zz() - eigValsR.yy() + 1e-16);

   omegaR.yx() = -omegaR.xy();
   omegaR.zx() = -omegaR.xz();
   omegaR.zy() = -omegaR.yz(); 
 }
 
 omega = ( eigVecs & omega & eigVecs.T() );

 B = ( eigVecs & B & eigVecs.T() );
 
}

void constitutiveEq::calcEig
(
  const volSymmTensorField& theta,
  volTensorField& vals,
  volTensorField& vecs
)
{
 
 // Eigen decomposition using a QR algorithm of Eigen library 
 Eigen::Matrix3d theta_eig(Eigen::Matrix3d::Zero(3,3));
 Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigSol;
 
 forAll(theta, cellI)
 {
   // Transfer theta from OF to Eigen        
   const symmTensor& thetaR = theta[cellI];  
       
   theta_eig(0,0)=thetaR.xx();
   theta_eig(1,1)=thetaR.yy();
   theta_eig(2,2)=thetaR.zz();

   theta_eig(0,1)=thetaR.xy();
   theta_eig(1,0)=thetaR.xy();

   theta_eig(0,2)=thetaR.xz();
   theta_eig(2,0)=thetaR.xz();

   theta_eig(1,2)=thetaR.yz();
   theta_eig(2,1)=thetaR.yz();
    
   // Compute eigenvalues/vectors in Eigen

   eigSol.compute(theta_eig);
   Eigen::Vector3d eival = eigSol.eigenvalues();
   Eigen::Matrix3d eivect = eigSol.eigenvectors();

   // Transfer eigenvalues/vectors from Eigen to OF 
   tensor& vecsR = vecs[cellI];  

   vecsR.xx()=eivect(0,0);
   vecsR.yx()=eivect(1,0);
   vecsR.zx()=eivect(2,0);

   vecsR.xy()=eivect(0,1);
   vecsR.yy()=eivect(1,1);      
   vecsR.zy()=eivect(2,1);      

   vecsR.xz()=eivect(0,2);  
   vecsR.yz()=eivect(1,2);
   vecsR.zz()=eivect(2,2);

   vals[cellI] *= 0.;
   vals[cellI].xx()=Foam::exp(eival(0));
   vals[cellI].yy()=Foam::exp(eival(1));
   vals[cellI].zz()=Foam::exp(eival(2));
  
 }

/*
 // Eigen decomposition using the iterative jacobi algorithm 

  forAll(theta, cellI)
    {
       int N=3;
       int NROT=0;
       jacobi(theta[cellI], N, vals[cellI], vecs[cellI], NROT);

    }

*/

}

void constitutiveEq::checkForStab
(
 const dictionary& dict
)
{
  stabOption_ = stabOptionNames_.read
  (
     dict.lookup("stabilization")
  ); 
}

void constitutiveEq::checkIfCoupledSolver
(
  const dictionary& dict,
  volSymmTensorField& tau
)
{

 const dictionary* cSDict = dict.subDictPtr("coupledSolvers");
 
 if (cSDict != NULL)
 {
   const dictionary* cSDictUpTau = cSDict->subDictPtr("Uptau");
   if (cSDictUpTau != NULL)
   {
    solveCoupled_ = 
    (
        readBool(cSDictUpTau->lookup("solveCoupledTau"))
     && readBool(cSDictUpTau->lookup("solveCoupledUp"))
    );
   }
 } 
 
 if (solveCoupled_)
 {
   coupledSolver& cps = U_.time().lookupObjectRef<coupledSolver>(word("Uptau."+tau.mesh().name()));  
   cps.insertField(tau);    
 }
 
}

tmp<volSymmTensorField> constitutiveEq::tauTotal() const
{
 volTensorField L(fvc::grad(U_));
   
 // Last term is only meaningfull for two-phase flows
 if (isGNF())
 {      
   return eta()*(symm(L+L.T()) - (2./3.)*tr(L)*symmTensor::I); 
 }
 else
 {
   return tau() + etaSThermo()*(symm(L+L.T()) - (2./3.)*tr(L)*symmTensor::I); 
 }
}

} //End namespace

// ************************************************************************* //
