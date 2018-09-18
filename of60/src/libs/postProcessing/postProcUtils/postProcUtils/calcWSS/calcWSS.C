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

#include "fvCFD.H"
 
#include "calcWSS.H"
#include "addToRunTimeSelectionTable.H"
#include "constitutiveModel.H"
#include "immiscibleConstitutiveTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ppUtils
{
    defineTypeNameAndDebug(calcWSS, 0);
    addToRunTimeSelectionTable(ppUtil, calcWSS, dictFS);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ppUtils::calcWSS::calcWSS
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U),
isTwoPhaseFlow_(dict.lookup("isTwoPhaseFlow")),
incPoly_(dict.lookup("includePolymericStresses")),
WSSmag_
(
 IOobject
  (
    "WSSmag",
    U.time().timeName(),
    U.mesh(),
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
 U.mesh(),
 dimensionedScalar
 (
  "zero",
  dimensionSet( 1, -1, -2, 0, 0, 0, 0),
  0.
 )
)  
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ppUtils::calcWSS::update()
{

if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//

     if (WSSmag_.time().outputTime()) // Only spends time when it is outputed
      {
        // Ensure zero WSS to start        
        forAll(WSSmag_.boundaryField(), patchI)
	{
	   WSSmag_.boundaryFieldRef()[patchI] *= 0.;
	}
	
	// Total extra-stress
        GeometricField<symmTensor, fvPatchField, volMesh> tau
        (
         IOobject
          (
           "tauTotal",
           U().mesh(),
           IOobject::NO_READ,
           IOobject::NO_WRITE,
           false
          ),
         U().mesh(),
         dimensionedSymmTensor("0",dimPressure,symmTensor::zero) 
        );
	 
	 
        if (isTwoPhaseFlow_) 
        {
          immiscibleConstitutiveTwoPhaseMixture& constEq_ = const_cast<immiscibleConstitutiveTwoPhaseMixture&>
          (
            U().mesh().lookupObject<immiscibleConstitutiveTwoPhaseMixture>("constitutiveProperties")
          );
          
          tau = constEq_.tauTotalMF(); 
          
          if (!incPoly_)
	  {
	    tau -= constEq_.tauMF(); 
	  } 
        } 
        else
        {
          constitutiveModel& constEq_ = const_cast<constitutiveModel&>
          (
            U().mesh().lookupObject<constitutiveModel>("constitutiveProperties")
          );

          tau = constEq_.tauTotal(); 
          
          if (!incPoly_)
	  {
	    tau -= constEq_.tau(); 
	  } 
        } 
      
        forAll(WSSmag_.boundaryField(), patchI)
	 {	          
           vectorField n(U().mesh().Sf().boundaryField()[patchI]/U().mesh().magSf().boundaryField()[patchI]);
 
           vectorField tracV(n & tau.boundaryField()[patchI]);

           vectorField nTracV(n * ( n & tracV));

           vectorField tTracV(tracV-nTracV);

           WSSmag_.boundaryFieldRef()[patchI] = mag(tTracV);
        }
       
      } // if outTime 

//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 


// ************************************************************************* //
