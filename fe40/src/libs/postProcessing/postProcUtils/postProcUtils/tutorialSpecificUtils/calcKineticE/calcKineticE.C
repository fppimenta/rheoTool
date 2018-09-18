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
 
#include "calcKineticE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace ppUtils
{
    defineTypeNameAndDebug(calcKineticE, 0);
    addToRunTimeSelectionTable(ppUtil, calcKineticE, dictFS);
}
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ppUtils::calcKineticE::calcKineticE
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U) 
{
   //- Create output streams 
    createFile();   
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
  
void Foam::ppUtils::calcKineticE::createFile()
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
        }      
     }
}

void Foam::ppUtils::calcKineticE::update()
{

if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//

           const volSymmTensorField& tau = U().mesh().lookupObject<volSymmTensorField>("tau");
           const dictionary& cttP = U().mesh().lookupObject<IOdictionary>("constitutiveProperties");
           dimensionedScalar lambda_(cttP.subDict("parameters").lookup("lambda"));
           dimensionedScalar etaP_(cttP.subDict("parameters").lookup("etaP"));
 
          // Compute kinetic energy

           int nCells = U().mesh().nCells(); 
           reduce(nCells, sumOp<int>());
           
           scalar ek = ( 0.5/nCells ) * gSum( mag( U().internalField() ) * mag( U().internalField() ) );
           scalar eElast = ( 0.5/nCells ) * (lambda_.value()/etaP_.value()) * gSum( tr( tau.internalField() ) ) ;          

           if (Pstream::master()) 
           {
             
             // Time | Average kinE | Average elasticE
             
               outS_() << U().mesh().time().value() << tab 
                       << ek << tab
                       << eElast << endl;            
           }
           
//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 
 

// ************************************************************************* //
