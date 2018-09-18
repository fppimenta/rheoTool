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

#include "calcFfl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ppUtils
{
    defineTypeNameAndDebug(calcFfl, 0);
    addToRunTimeSelectionTable(ppUtil, calcFfl, dictFS);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ppUtils::calcFfl::calcFfl
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U),
shaft_(U.mesh().boundaryMesh().findPatchID("shaft")),
piston_(U.mesh().boundaryMesh().findPatchID("piston"))
{
  //- Create output streams 
    createFile();
   
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
void Foam::ppUtils::calcFfl::createFile()
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
          outS_.reset(new OFstream(ppDir_/"Ffl.txt"));    
        }      
     }
}

void Foam::ppUtils::calcFfl::update()
{

if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//
 
           const volScalarField& p = U().mesh().lookupObject<volScalarField>("p");
           const dictionary& constDict = U().mesh().lookupObject<IOdictionary>("constitutiveProperties");
           dimensionedScalar rho_(constDict.subDict("parameters").lookup("rho"));
           const volScalarField& eta_ = U().mesh().lookupObject<volScalarField>("eta");
  
           // Compute Ffl
 
           volTensorField L = fvc::grad(U());

	   volSymmTensorField F = symm( L + L.T() ) * eta_ - p * symmTensor::I * rho_;

           vector Fpatch =  gSum( ( -U().mesh().boundaryMesh()[piston_].faceAreas() ) & F.boundaryField()[piston_] )
                          + gSum( ( -U().mesh().boundaryMesh()[shaft_].faceAreas() ) & F.boundaryField()[shaft_] );
           
           // Note: the angle of our wedge is 5º, but we want to report the force acting on all the surface. 
           // We do not account for the real curvature (the wall of the outer cage is straight rather than curved).                
           Fpatch *= 360./5.;               
           
           if (Pstream::master()) 
           {
           
             // col 0: Time
               outS_() << U().mesh().time().value() << tab; 
                             
              // col 1: Cd             
               outS_() << Fpatch.x() << endl; 
           
           }
        
//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 


// ************************************************************************* //
