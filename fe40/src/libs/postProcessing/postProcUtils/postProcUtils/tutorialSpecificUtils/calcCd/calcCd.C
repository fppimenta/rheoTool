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

#include "calcCd.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ppUtils
{
    defineTypeNameAndDebug(calcCd, 0);
    addToRunTimeSelectionTable(ppUtil, calcCd, dictFS);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ppUtils::calcCd::calcCd
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U),
cyl_(U.mesh().boundaryMesh().findPatchID("cylinder"))  
{
  //- Create output streams 
    createFile();
   
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
void Foam::ppUtils::calcCd::createFile()
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
          outS_.reset(new OFstream(ppDir_/"Cd.txt"));    
        }      
     }
}

void Foam::ppUtils::calcCd::update()
{

if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//


           const volSymmTensorField& tau = U().mesh().lookupObject<volSymmTensorField>("tau");
           const volScalarField& p = U().mesh().lookupObject<volScalarField>("p");
           const dictionary& constDict = U().mesh().lookupObject<IOdictionary>("constitutiveProperties");
           dimensionedScalar rho_(constDict.subDict("parameters").lookup("rho"));
           dimensionedScalar etaS_(constDict.subDict("parameters").lookup("etaS"));
           dimensionedScalar etaP_(constDict.subDict("parameters").lookup("etaP"));
  
          // Compute cd
 
           volTensorField L = fvc::grad(U());

	   volSymmTensorField F = tau + symm( L + L.T() ) * etaS_ - p * symmTensor::I * rho_;

           vector Fpatch = gSum( ( -U().mesh().boundaryMesh()[cyl_].faceAreas() ) & F.boundaryField()[cyl_] )/(etaS_ + etaP_).value();
  
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
