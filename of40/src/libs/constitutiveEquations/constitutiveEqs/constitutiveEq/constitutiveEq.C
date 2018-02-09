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
#include <Eigen/Dense> // For eigen decomposition
#include "jacobi.H"    // Only required for jacobi decomposition

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constitutiveEq, 0);
    defineRunTimeSelectionTable(constitutiveEq, dictionary);



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
    stabMeth_(0)
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
           if (stabMeth_ == 0) // none
	     {
	     
	      return
	       (
		  fvc::div(tau()/rho(), "div(tau)")
		+ fvm::laplacian(etaS()/rho(), U, "laplacian(eta,U)")
	       );

	     }
	    else if (stabMeth_ == 1) // BSD 
	     {
	      
	      return
	       (
		  fvc::div(tau()/rho(), "div(tau)")
		- fvc::laplacian(etaP()/rho(), U, "laplacian(eta,U)")
		+ fvm::laplacian( (etaP()+ etaS())/rho(), U, "laplacian(eta,U)")
	       );

	     }   
	    else // coupling
	     {
	      
	      return
	       (
		  fvc::div(tau()/rho(), "div(tau)")
		- (etaP()/rho()) * fvc::div(fvc::grad(U))
		+ fvm::laplacian( (etaP() + etaS())/rho(), U, "laplacian(eta,U)")
	       );

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
           if (stabMeth_ == 0) // none
	     {
	     
	      return
	       (        
		  fvm::laplacian(etaS()*alpha, U, "laplacian(eta,U)")
		+ fvc::div(etaS()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
	       );

	     }
	   else if (stabMeth_ == 1) // BSD 
	     {
	      
	      return
	       (
		- fvc::laplacian(etaP()*alpha, U, "laplacian(eta,U)")
		+ fvm::laplacian( (etaP()+ etaS())*alpha, U, "laplacian(eta,U)")
		+ fvc::div(etaS()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
	       );

	     } 
	   else // coupling
	     {
	      
	      return
	       (        
		- fvc::div( (etaP()*alpha) * fvc::grad(U), "div(grad(U))")
		+ fvm::laplacian( (etaP() + etaS())*alpha, U, "laplacian(eta,U)")
		+ fvc::div(etaS()*alpha*dev2(T(fvc::grad(U))), "div(eta*alpha*dev2(T(gradU)))")
	       );

	     }       
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
        
         B[cellI].xx()=M[cellI].xx();
         B[cellI].yy()=M[cellI].yy(); 
         B[cellI].zz()=M[cellI].zz();
 
         
         omega[cellI].xy() = ( eigVals[cellI].yy()*M[cellI].xy() + eigVals[cellI].xx()*M[cellI].yx() ) / ( eigVals[cellI].yy() - eigVals[cellI].xx() + 1e-16);
         omega[cellI].xz() = ( eigVals[cellI].zz()*M[cellI].xz() + eigVals[cellI].xx()*M[cellI].zx() ) / ( eigVals[cellI].zz() - eigVals[cellI].xx() + 1e-16);
         omega[cellI].yz() = ( eigVals[cellI].zz()*M[cellI].yz() + eigVals[cellI].yy()*M[cellI].zy() ) / ( eigVals[cellI].zz() - eigVals[cellI].yy() + 1e-16);

         omega[cellI].yx() = -omega[cellI].xy();
         omega[cellI].zx() = -omega[cellI].xz();
         omega[cellI].zy() = -omega[cellI].yz(); 
      }
 
      omega = ( eigVecs & omega & eigVecs.T() );
   
      B = ( eigVecs & B & eigVecs.T() );


}

void constitutiveEq::calcEig(const volSymmTensorField& theta, volTensorField& vals, volTensorField& vecs)
{
 
 // Eigen decomposition using a QR algorithm of Eigen library 

    Eigen::Matrix3d theta_eig(Eigen::Matrix3d::Zero(3,3));
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigSol;
 
    forAll(theta, cellI)
     {
       // Transfer theta from OF to Eigen
    
       theta_eig(0,0)=theta[cellI].xx();
       theta_eig(1,1)=theta[cellI].yy();
       theta_eig(2,2)=theta[cellI].zz();

       theta_eig(0,1)=theta[cellI].xy();
       theta_eig(1,0)=theta[cellI].xy();

       theta_eig(0,2)=theta[cellI].xz();
       theta_eig(2,0)=theta[cellI].xz();

       theta_eig(1,2)=theta[cellI].yz();
       theta_eig(2,1)=theta[cellI].yz();
    
       // Compute eigenvalues/vectors in Eigen

       eigSol.compute(theta_eig);
       Eigen::Vector3d eival = eigSol.eigenvalues();
       Eigen::Matrix3d eivect = eigSol.eigenvectors();

       // Transfer eigenvalues/vectors from Eigen to OF 

       vecs[cellI].xx()=eivect(0,0);
       vecs[cellI].yx()=eivect(1,0);
       vecs[cellI].zx()=eivect(2,0);

       vecs[cellI].xy()=eivect(0,1);
       vecs[cellI].yy()=eivect(1,1);      
       vecs[cellI].zy()=eivect(2,1);      

       vecs[cellI].xz()=eivect(0,2);  
       vecs[cellI].yz()=eivect(1,2);
       vecs[cellI].zz()=eivect(2,2);


       vals[cellI]=tensor::zero;
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

tmp<volScalarField> constitutiveEq::strainRate()
{
   return ( sqrt(2.0)*mag( symm(fvc::grad(U_)) ) );
}

void constitutiveEq::checkForStab(const dictionary& dict)
{
  // Allow coupling by default
  word stab_(dict.lookupOrDefault<word>("stabilization", "coupling"));
  
  if (stab_ == "none")
   {
     stabMeth_ = 0;
     Info << "Selected stabilization method: none.\n";
   }
  else if (stab_ == "BSD")
   {
     stabMeth_ = 1;
     Info << "Selected stabilization method: BSD.\n";
   }
  else if (stab_ == "coupling")
   {
     stabMeth_ = 2;
     Info << "Selected stabilization method: coupling.\n";
   }
  else
   {
       FatalErrorIn("constitutiveEq::checkForStab(const dictionary& dict)\n")
            << "\nThe stabilizatin method specified does not exist.\n"
            << "\nAvailable methods are:\n"
            << "\n. none" <<"\n. BSD" << "\n. coupling" 
            << abort(FatalError);
   }   
}

} //End namespace

// ************************************************************************* //
