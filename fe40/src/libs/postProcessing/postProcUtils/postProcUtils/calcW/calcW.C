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
 
#include "calcW.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calcW, 0);
    addToRunTimeSelectionTable(ppUtil, calcW, dictFS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcW::calcW
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

void Foam::calcW::createFile()
{

  // Not prepared to post-process decomposed cases (TODO) 
     
  if (Pstream::master()) 
   {            
         
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "void Foam::calcW::"
            "createFile()"
        )<< "This post-processing util is not prepared to be run in parallel. "
         << "Please disable it in fvSolution dict to keep running the case in parallel." <<
          nl << abort(FatalError);
    }
            
   }
   
  // Note: ppDir_ is the directory for the output of postPorcessing data.
  // It is defined in the base class and mkdir is called there also.
  
  //- Create file for fluxes   
    if (outS_.empty())
     {     
      // Create name of dir
        if (Pstream::master()) 
        {
          // Open the file
          outS_.reset(new OFstream(ppDir_/"DropWidth.txt"));    
        }      
     }
}

void Foam::calcW::update()
{

if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//


           const volScalarField& alpha1 = U().mesh().lookupObject<volScalarField>("alpha1");
           surfaceScalarField alpha1f = fvc::snGrad(alpha1);
           
           const volVectorField& C = U().mesh().C();
                 
           scalar rowy(C[0].y());
           bool endly(false);
           scalar maxInter(0.), inter(0.);
           
           forAll(C,idx)
            {
                if (mag(C[idx].y()-rowy)<SMALL)
                {
		      if ((C[idx].x()>0.) && (!endly))
		       {
		         if ((alpha1[idx]-0.5)*(alpha1[idx-1]-0.5)<0.)
		          {
		             scalar y1=alpha1[idx-1];
		             scalar y2=alpha1[idx];
		             scalar x1=C[idx-1].x();
		             scalar x2=C[idx].x();
		             
		             inter = (0.5-y1)/( (y2-y1)/(x2-x1) ) + x1;
		             
		             endly = true;
		             if (inter>maxInter) {maxInter = inter;}
		            
		          }
		       
		       }  
		       
		  }  
		  else
		  {
		 
		    endly = false;
		    rowy = C[idx].y();
		    		    		  
		  }       
            
            }
   
         
          if (Pstream::master()) 
           {
           
             // col 0: Time
               outS_() << U().mesh().time().value() << tab; 
                             
              // col 1: The (max) width of the drop            
               outS_() << maxInter << endl; 
           
           }
         
         
//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 

// ************************************************************************* //
