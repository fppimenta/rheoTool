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
#include "IFstream.H"
#include "OFstream.H"

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
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(dict, U),
cyl_(U.mesh().boundaryMesh().findPatchID("cylinder"))  
{
   createFile(); 
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


   
void Foam::calcW::createFile()
{

  if (outS_.empty())
   {
   
      // Create name of dir
        if (Pstream::master()) 
        {
            fileName dir;
          
            if (Pstream::parRun())
            {
                dir = U().time().path()/".."/"DropWidth.txt";
            }
            else
            {
                dir = U().time().path()/"DropWidth.txt";
            }   
            
           // Open the file
            outS_.reset(new OFstream(dir));          
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
          
          counter_ = 0;
    }

 }
              
}


// ************************************************************************* //
