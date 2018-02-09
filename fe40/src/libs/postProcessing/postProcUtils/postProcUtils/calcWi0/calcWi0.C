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
 
#include "calcWi0.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calcWi0, 0);
    addToRunTimeSelectionTable(ppUtil, calcWi0, dictFS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcWi0::calcWi0
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U),
cellC_(U.mesh().findCell(vector(0,0,0.5)))  
{
   //- Create output streams 
    createFile(); 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::calcWi0::createFile()
{
  // Note: ppDir_ is the directory for the output of postPorcessing data.
  // It is defined in the base class and mkdir is called there also.
    
  //- Create file for fluxes   
    if (outS_.empty())
     {     
        if (Pstream::parRun()) 
        {
           for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
              if (procI == Pstream::myProcNo())
                {
                  if (cellC_!=-1)
                   {
                    // Open the file
                    outS_.reset(new OFstream(ppDir_/"Wi0.txt"));                   
                   }
                }
            }          
        }      
        else
        {
         // Open the file
          outS_.reset(new OFstream(ppDir_/"Wi0.txt"));  
        }
     }
}

void Foam::calcWi0::update()
{

 if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//

           const dictionary& constDict = U().mesh().lookupObject<IOdictionary>("constitutiveProperties");
           dimensionedScalar lambda_(constDict.subDict("parameters").lookup("lambda"));
         
           volTensorField gradU = fvc::grad(U());
          
           if (Pstream::parRun())
   	     {
   	        for (label procI = 0; procI < Pstream::nProcs(); procI++)
	         {
                   if (procI == Pstream::myProcNo())
                     {
                       if (cellC_!=-1)
                       {
                           // col 0: Time
                           outS_()<< U().mesh().time().value() << tab; 
                             
                           // col 1: Wi0             
                           outS_()<< Foam::sqrt( gradU[cellC_].xx() * gradU[cellC_].xx() + gradU[cellC_].yx() * gradU[cellC_].xy() ) * lambda_.value()<< endl; 
                       }                   
                     }
                 }               
             }
            else
             {
              // col 0: Time
               outS_()<< U().mesh().time().value() << tab; 
                             
              // col 1: Wi0             
               outS_()<< Foam::sqrt( Foam::max(1e-12, gradU[cellC_].xx() * gradU[cellC_].xx() + gradU[cellC_].yx() * gradU[cellC_].xy() ) ) * lambda_.value() << endl; 
             }

//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 


// ************************************************************************* //
