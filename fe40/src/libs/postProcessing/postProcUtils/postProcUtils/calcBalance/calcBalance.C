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
 
#include "calcBalance.H"
#include "addToRunTimeSelectionTable.H"
#include "EDFEquation.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calcBalance, 0);
    addToRunTimeSelectionTable(ppUtil, calcBalance, dictFS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcBalance::calcBalance
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U),
D_(),
z_(),
specL_(),
T_()
{
 
   //- Read fluid properties      
    const dictionary& elecDict = U.mesh().lookupObject<IOdictionary>("electricProperties");
    
    dimensionedScalar  dT_(elecDict.subDict("parameters").lookup("T")); T_ = dT_.value();
        
    PtrList<entry> specEntries_(elecDict.subDict("parameters").lookup("species"));
    
    D_.clear(); z_.clear(); specL_.clear();
    D_.setSize(specEntries_.size(), 0.); z_.setSize(specEntries_.size(), 1); specL_.setSize(specEntries_.size(), " ");
    
    forAll (specEntries_, specI)
    {  
    
       dimensionedScalar dD_(specEntries_[specI].dict().lookup("D")); D_[specI]=dD_.value();
       dimensionedScalar dz_(specEntries_[specI].dict().lookup("z")); z_[specI]=dz_.value();
       specL_[specI] = specEntries_[specI].keyword();
    }
   
   //- Create output streams 
    createFile();
   
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
void Foam::calcBalance::createFile()
{
  // Note: ppDir_ is the directory for the output of postPorcessing data.
  // It is defined in the base class and mkdir is called there also.
  
  //- Create file for fluxes   
    if (outS_.empty())
     {     
      // Create name of dir
        if (Pstream::master()) 
        {
          // Open the file
          outS_.reset(new OFstream(ppDir_/"Balance.txt"));    
           
          // Header
          outS_() << "Time" << tab;
           
          forAll(specL_, si)
            {                  
              outS_() << word("Average_concentration_of_" + specL_[si] + "_mol.m-3") << tab 
                      << word("Average_Flux_of_" + specL_[si] + "_mol.m-2.s-1") << tab;            
            }           
                    
          outS_() << endl; 
        }      
     }
}

void Foam::calcBalance::update()
{

if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//


    scalar eK_(Foam::EDFEquation::eK_.value());
    scalar kbK_(Foam::EDFEquation::kbK_.value());
    
    const volScalarField& psi = U().mesh().lookupObject<volScalarField>("psi");
    const volScalarField& phiE = U().mesh().lookupObject<volScalarField>("phiE");
    const surfaceScalarField& phi = U().mesh().lookupObject<surfaceScalarField>("phi");

    // Write time
    if (Pstream::master()) 
     {       	     
	outS_() << U().mesh().time().value() << tab; 
     }
     
    forAll(specL_, si)
     {  
      
      const volScalarField& specI = U().mesh().lookupObject<volScalarField>(specL_[si]);
      scalar zi(z_[si]);
      scalar Di(D_[si]);
      
     // Compute average concentration in the bulk
     
      scalar med(0.); 
      scalar v(0.);
      forAll(psi, ci)
       {
        med += (specI[ci])*U().mesh().V()[ci];
        v += U().mesh().V()[ci];
       } 
     
      reduce(med, sumOp<scalar>()); 
      reduce(v, sumOp<scalar>()); 
    
    // Compute average flux on boundaries
    
     scalar fluxBC(0.);
     scalar ttArea(0.);
     
     forAll(psi.boundaryField(), pI)
      {
      
       if (!psi.boundaryField()[pI].coupled())
        {          
          const fvPatch& pt = U().boundaryField()[pI].patch();
           
           {
             scalarField J = 
             specI.boundaryField()[pI] * phi.boundaryField()[pI]
	    -Di*
	     (
	        specI.boundaryField()[pI].snGrad()
	      + (eK_*zi/(kbK_*T_)) * (psi.boundaryField()[pI].snGrad() + phiE.boundaryField()[pI].snGrad()) * specI.boundaryField()[pI]
	     )*pt.magSf();
	
	     fluxBC += gSum(J);
	     ttArea += gSum(pt.magSf());
           }
	                 
        } // if patch not coupled     
              
     } // forAll patches
     
     // Write              
      if (Pstream::master()) 
       {             	     
	  outS_() << med/v << tab << fluxBC/ttArea << tab; 		                      
       }
	       
   } // forAll speci  
        
   if (Pstream::master()) 
     {       	     
	outS_() << endl; 
     }  
      
           
//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 
 

// ************************************************************************* //
