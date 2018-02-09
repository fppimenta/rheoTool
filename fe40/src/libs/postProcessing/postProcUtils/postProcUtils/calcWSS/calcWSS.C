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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calcWSS, 0);
    addToRunTimeSelectionTable(ppUtil, calcWSS, dictFS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcWSS::calcWSS
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U),
isVE_(dict.lookup("isViscoelastic")),
incSolv_(dict.lookup("includeSolventStresses")),
incPoly_(dict.lookup("includePolymericStresses")),
etaSWW_("etaSWW", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.),
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
{

// Compute/read solvent viscosity in case of viscoelastic fluid
// considering the possibility of being a multimode model

if (isVE_ && incSolv_)
 {
 
   IOdictionary constitutiveProperties
    (
        IOobject
        (
            "constitutiveProperties",
            U.time().constant(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );
    
   word CM = constitutiveProperties.subDict("parameters").lookup("type");
   
   if(CM=="multiMode")
     {
         PtrList<entry> modelEntries(constitutiveProperties.subDict("parameters").lookup("models"));

         scalar etaCum(0.0);

            forAll (modelEntries, modelI)
            {
               dimensionedScalar etaI(modelEntries[modelI].dict().lookup("etaS"));
               etaCum += etaI.value();	      
            }
         etaSWW_.value() = etaCum;
     }
    else
     {
        dimensionedScalar etaStmp(constitutiveProperties.subDict("parameters").lookup("etaS"));                 
        etaSWW_.value() = etaStmp.value();               
     }     
 }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcWSS::update()
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
		WSSmag_.boundaryField()[patchI] *= 0.;
	 }
        
        // Case GNF
         
        if (!isVE_)
         {   
           const volScalarField& eta_ = U().mesh().lookupObject<volScalarField>("eta");
           
           volTensorField L = fvc::grad(U());
           volTensorField tau_ = eta_ * ( L + L.T() ) ;

           forAll(WSSmag_.boundaryField(), patchI)
	    {
	      vectorField n(U().mesh().Sf().boundaryField()[patchI]/U().mesh().magSf().boundaryField()[patchI]);
	      
	      WSSmag_.boundaryField()[patchI] = mag( n & tau_.boundaryField()[patchI] );
	    }
	 
	 }   
       
       // Case VE
       	 
	else 
	 {
	   //- Polymeric contribution (needs to remove normal stresses in a general case)
	   
	   if (incPoly_)
	    {
	      forAll(WSSmag_.boundaryField(), patchI)
	       {
	          const volSymmTensorField& tau_ = U().mesh().lookupObject<volSymmTensorField>("tau");
           
                  vectorField n(U().mesh().Sf().boundaryField()[patchI]/U().mesh().magSf().boundaryField()[patchI]);
 
                  vectorField tracV(n & tau_.boundaryField()[patchI]);

                  vectorField nTracV(n * ( n & tracV));

                  vectorField tTracV(tracV-nTracV);

                  WSSmag_.boundaryField()[patchI] = mag(tTracV);
	       }
	    }
	   
	   //- Solvent contribution
	   
	   if (incSolv_)
	    {
	   
              volTensorField L = fvc::grad(U());
              volTensorField tau_ = etaSWW_ * ( L + L.T() ) ;

              forAll(WSSmag_.boundaryField(), patchI)
	       {
	         vectorField n(U().mesh().Sf().boundaryField()[patchI]/U().mesh().magSf().boundaryField()[patchI]);
	      
	         WSSmag_.boundaryField()[patchI] += mag( n & tau_.boundaryField()[patchI] );
	       }
	    }	 
	 }
	 
       
      } // if outTime 

//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 


// ************************************************************************* //
