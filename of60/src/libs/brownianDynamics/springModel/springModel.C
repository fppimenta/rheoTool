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
 
#include <Eigen/Dense>  
#include "springModel.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(springModel, 0);
    defineRunTimeSelectionTable(springModel, dictFS);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

springModel::springModel
(
    const dictionary& dict,
    const volVectorField& U,
    sPCloudInterface& sPCI
)
:
U_(U),
sPCref_(sPCI),
maxIter_(dict.subDict("springModelProperties").lookupOrDefault<int>("maxIter", 10)),
relTol_(dict.subDict("springModelProperties").lookupOrDefault<scalar>("relTol", 1e-6)),
tresholdF_(dict.subDict("springModelProperties").lookupOrDefault<scalar>("tresholdF", .95)),
solverType_(dict.subDict("springModelProperties").lookupOrDefault<word>("solver", "TDMA")),
timeSch_(dict.subDict("springModelProperties").lookupOrDefault<word>("timeScheme", "semiImplicit")),
linkM_(sPCI.linkM()),       
nMolc_(sPCI.nMolc()),
mx_(sPCI.mx()),       
mU_(sPCI.mU()),              
mAct_(sPCI.mAct()),      
mIds_(sPCI.mIds()),
mD_(sPCI.mD()),
mSigma_(sPCI.mSigma()),
mSpr_(sPCI.mSpr()),
Nks_(sPCI.Nks()),
D_(sPCI.D()),
Ls_(sPCI.Ls()),
isTethered_(sPCI.isTethered()),
isHI_(sPCI.isHI())
{
   
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::springModel::fSpring()
{
  // What follows to the tracking is mU_, but we also update mx_ to know
  // if the springs would be overstretched
  
  scalar dt(U().mesh().time().deltaTValue());

  forAll(mU_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
          mU_[mi] = fSpringI(mx_[mi], mi, true, true);
          mx_[mi] += dt*mU_[mi];          
       }
   }

}

// TDMA method solving for Ax+B=0, where A is tridiagonal
// and diagonal dominant. On return, B is X. B is scalarField
// and not a matrix to avoid conversion between types on return
void springModel::TDMA(scalarSquareMatrix& A, scalarField& B) 
{
  int n(B.size());  
  
  A[0][1] /= A[0][0];
  B[0] /= -A[0][0];  

  for (int i = 1; i < n -1; i++) 
  {
     A[i][i+1] /= A[i][i] - A[i][i-1]*A[i-1][i];
     
     B[i] = (-B[i] - A[i][i-1]*B[i-1]) / (A[i][i] - A[i][i-1]*A[i-1][i]);      
  }

  B[n-1] = (-B[n-1] - A[n-1][n-2]*B[n-2]) / (A[n-1][n-1] - A[n-1][n-2]*A[n-2][n-1]);
  

  for (int i = n-2; i>=0; i--)
  {
    B[i] -= A[i][i+1]*B[i+1];
  }
}

void springModel::deleteMolecule(label mi)
{

  sPCref_.deleteMolecule(mi);
  
  // Uncomment/comment following sections to control the behavior
  // after a spring is overstreched: warning or error.

/*
  FatalErrorIn("springModel::deleteMolecule()")
  << "\nSprings are overstretched.\n"
  << "Try to reduce the lagrangian time-step and/or to increase the number of cycles in the implicit loop."
  << exit(FatalError);
*/

  WarningIn("springModel::deleteMolecule()")
  << "\nSprings in molecule " << mi << " are overstretched.\n"
  << " The molecule will be deleted." 
  << "Try to reduce the lagrangian time-step and/or to increase the number of cycles in the implicit loop." 
  << endl;
														
}

void springModel::checkSpringsLength
(
  const PtrList<Field<vector > >& mxStar,
  const PtrList<Field<vector > >& mx0
) 
{
 // Explicit scheme just checks if the spring is bounded.
 // Hookean model is unbounded and has its own checkSpringsLength().
 if (timeSch_ == "explicit")
  {
    forAll(mx_, mi)
     {
      if (mAct_[mi][0] != -1)
       {  
        label gI(mIds_[mi][0][2]); // All beads belong to the same group, thus check for the first bead and save              
        
        forAll(mSpr_[mi], bi)
 	 {
           label& b0 = mSpr_[mi][bi][0];
           label& b1 = mSpr_[mi][bi][1];  
                
           if ( ( mag(mx_[mi][b0] - mx_[mi][b1])/Ls_[gI] ) > 1. )
            {             
              deleteMolecule(mi);
            }
          }
       }
     } 
  }
 else if (timeSch_ == "semiImplicit")  
  {
  
   forAll(mx_, mi)
   {
       if (mAct_[mi][0] != -1)
       {            
          label gI(mIds_[mi][0][2]); // All beads belong to the same group, thus check for the first bead and save
         
           forAll(mSpr_[mi], bi)
 	   {
             label& b0 = mSpr_[mi][bi][0];
             label& b1 = mSpr_[mi][bi][1];
              
             // If one single spring in a molecule is overstretched, then the molecule needs implicit evaluation
             if ( ( mag(mx_[mi][b0] - mx_[mi][b1])/Ls_[gI] ) > tresholdF_ )
             {             
                implicitfSpring(mi, mxStar[mi]/Ls_[gI], mx0[mi]/Ls_[gI]);
            
                // After the implicit correction check if the springs are really bounded and 
                // remove the molecule if not
                forAll(mSpr_[mi], bii)
 	        {
                  label& b00 = mSpr_[mi][bii][0];
                  label& b11 = mSpr_[mi][bii][1];
                
                  if ( ( mag(mx_[mi][b00] - mx_[mi][b11])/Ls_[gI] ) > 1. )
                   {             
                      deleteMolecule(mi);
                   }
                }
               
                break;             
             }                          
 	   }
       } 	   
   } 	    	              
  }
 else
  {
   FatalErrorIn("springModel::springModel()")
    << "\nUnknow time integration scheme:" << timeSch_ 
    << "\nValid time integration schemes are: \n .explicit \n .semiImplicit\n"
    << exit(FatalError);
  } 
    
}

void springModel::implicitfSpring
(
  label mi, 
  const Field<Foam::vector>& xStar, 
  const Field<Foam::vector>& x0
)
{
  scalar dt(U().mesh().time().deltaTValue());
   
  label n(xStar.size());
  
  label gI(mIds_[mi][0][2]);
  scalar Ls(Ls_[gI]);
   
  // Create vector of positions from the initial guess
  vectorField xn(x0);     
    
  // J- Jacobian; f- the function
  // If only accounting with spring force and no HI, J is a tridiagonal symmetric matrix.
  // To keep generality make J full.
  List<scalarSquareMatrix> J(3, scalarSquareMatrix(n,n,Zero)); 
  List<scalarField> f(3, scalarField(n,0.));  
  
  scalar error(GREAT);
  int nIter(0);
  while (error>relTol_ && nIter<maxIter_)
   { 
     vectorField  xStarPlusHI = xStar;
     if (isHI_)
       xStarPlusHI += fSpringI(xn, mi, false, false)*dt/Ls;
       
     // Build J and f
     for (int cmpI=0; cmpI<3 ; cmpI++)
      {
       
       // J[cmpI] and f[cmpI] come on return
       jacobianSIM
       (
         mi,
         cmpI,
         xn,
         J[cmpI]
       );
        
       fSIM
       (
         mi,
         cmpI,
         xStarPlusHI.component(cmpI),  // Update the off-diag components of spring force due to HI
         xn, 
         f[cmpI]
       );
      
       // Correct if tethered (needs improvement)
       if (isTethered_)
        {
          // dx/dt = 0 for the fixed bead
          f[cmpI][0] = 0.;
          
          // follows from previous
          J[cmpI][0][0] = 1.;
          for (int j=1; j<J[cmpI].m(); j++)
           {
             J[cmpI][0][j] = 0.;
           }
        }  
             
      }
  
     // Compute new xn 
     vectorField xOld(xn); // Save x to compute relative error
     # include "solveSystem.H"
      
     // Compute the error: max variation in the dimLess spring vector
     error = max(mag(xOld - xn));
 
     nIter++;
     
     if (debug)
      {
          Pout<< "sPCloudInterface::implicitfSpring()" << nl
              << "Entering implicit correction loop for molecule: " << mi << nl
              << "Error afer iteration " << nIter << " is: " << error << endl;
      }

   }
   
  // Warning if exit with rel error below a critical treshold. Change to warn
  // if exit after maxIters are exceeded? 
  if (error>.01 && relTol_<.01)
   {
      WarningIn("springModel::implicitfSpring()")
      << "\nNewton-Raphson process only converged up to a relative tolerance "
      << "of " << error << " for molecule " << mi << "."
      << endl;
   }
  
  // Set fields according to the newly computed positions        
  // mU is the only field realy needed for tracking purposes, used in trackToFace()
  mU_[mi] = (xn-xStar)*Ls/dt;
  
  // mx is set to pass the re-check for overstretch after call to implicitfSpring()
  // and also to be used, if needed, as initial guess in the next time-step, since
  // it is conservative 
  mx_[mi] = xn*Ls;
   
  // Correct if tethered
  if (isTethered_)
   {
     mU_[mi][0] *= 0.;
     mx_[mi][0] = x0[0]*Ls;
   }    
}

// * * * * * * * * * * * * * *    Selector    * * * * * * * * * * * * * * //

autoPtr<springModel> springModel::New
(
    const dictionary& dict,
    const volVectorField& U,
    sPCloudInterface& sPCI
)
{
    word typeName = dict.subDict("springModelProperties").lookup("springModel");

    Info<< "Selecting spring model: " << typeName << endl;

    dictFSConstructorTable::iterator cstrIter =
        dictFSConstructorTablePtr_->find(typeName);

    if (cstrIter == dictFSConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "springModel::New(const dictionary& dict, const volVectorField&)"
        )   << "Unknown springModel type " << typeName
            << endl << endl
            << "Valid springModel types are :" << endl
            << dictFSConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<springModel>(cstrIter()(dict, U, sPCI));
}


} //End namespace

// ************************************************************************* //
