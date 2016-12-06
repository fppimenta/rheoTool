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

#include "calcCd.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calcCd, 0);
    addToRunTimeSelectionTable(ppUtil, calcCd, dictFS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcCd::calcCd
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


   
void Foam::calcCd::createFile()
{

  if (outS_.empty())
   {
   
      // Create name of dir
        if (Pstream::master()) 
        {
            fileName dir;
          
            if (Pstream::parRun())
            {
                dir = U().time().path()/".."/"Cd.txt";
            }
            else
            {
                dir = U().time().path()/"Cd.txt";
            }   
            
           // Open the file
            outS_.reset(new OFstream(dir));          
        }      
   }
}

void Foam::calcCd::update()
{

 if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

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
           
           counter_ = 0;
    }

 }

}


// ************************************************************************* //
