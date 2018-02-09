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
 
#include "calcVortexL.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(calcVortexL, 0);
    addToRunTimeSelectionTable(ppUtil, calcVortexL, dictFS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcVortexL::calcVortexL
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U
)
:
ppUtil(name, dict, U), 
vortTop_(U.mesh().boundaryMesh().findPatchID("wall_vorttop")), 
vortDown_(U.mesh().boundaryMesh().findPatchID("wall_vortdown")), 
lipTop_(U.mesh().boundaryMesh().findPatchID("wall_liptop")),
lipDown_(U.mesh().boundaryMesh().findPatchID("wall_lipdown")) 
{
   //- Create output streams 
    createFile();   
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


   
void Foam::calcVortexL::createFile()
{

 // Not prepared to post-process decomposed cases (TODO) 
     
  if (Pstream::master()) 
   {            
         
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "void Foam::calcVortexL::"
            "createFile()"
        )<< "This post-processing util is not prepared to be run in parallel. "
         << "Please disable it in fvSolution dict to keep running the case in parallel." <<
          nl << abort(FatalError);
    }
            
   }

  if (outS_.empty())
   {
          // Open the file
            outS_.reset(new OFstream(ppDir_/"Xr_top.txt"));                  
   }
   
   if (outS1_.empty())
   {                        
         // Open the file
            outS1_.reset(new OFstream(ppDir_/"Xr_down.txt"));              
   }
   
   if (outS2_.empty())
   {
         // Open the file
            outS2_.reset(new OFstream(ppDir_/"Xl_top.txt"));           
   }
   
   if (outS3_.empty())
   {
         // Open the file
            outS3_.reset(new OFstream(ppDir_/"Xl_down.txt"));      
   }
}

void Foam::calcVortexL::update()
{
if (enabled_)
 {
   counter_++;
   
   if (counter_ > nEval_)
    {

//****  User-defined function STARTS here *****//


           // col 0: Time
           outS_()  << U().mesh().time().value() << tab; 
           outS1_() << U().mesh().time().value() << tab; 
           outS2_() << U().mesh().time().value() << tab; 
           outS3_() << U().mesh().time().value() << tab; 

           const volVectorField& C = U().mesh().C();
           
// Xr_top
        {
           const polyPatch& cPatchvort = U().mesh().boundaryMesh()[vortTop_];
            
          // Define reference parameters 
 
           vector refPoint(0., 4., 0.5); // Reference point to zero the vortex length
           vector refDir(1., 0., 0.); // Vector aligned with the wall

          // Compute vortex length based on the point of velocity inversion

           int index(0);
           scalar uPrev=0.0; vector CPrev(0., 0., 0.);
           vector refDirU(refDir/mag(refDir));
      
           forAll(cPatchvort, facei )       
              {
                label  curCell = cPatchvort.faceCells()[facei];
                scalar uCmp = (U()[curCell] & refDirU);
               
                if (uPrev*uCmp<0.0)
                 {
                   vector r_curCell = -uCmp * ( CPrev - C[curCell] ) / (uPrev - uCmp) + C[curCell];
  
                   // col N: Time
                   outS_() << mag( ( (r_curCell - refPoint) & refDirU ) ) << tab; 
                   
                   index++;
                 } 
       
                uPrev = uCmp;
                CPrev = C[curCell];
            } 
            
           outS_() << endl;
           
         }

// Xr_down
        {
          const polyPatch& cPatchvort = U().mesh().boundaryMesh()[vortDown_];
        
          // Define reference parameters 
 
           vector refPoint(0., -4., 0.5); // Reference point to zero the vortex length
           vector refDir(1., 0., 0.); // Vector aligned with the wall

          // Compute vortex length based on the point of velocity inversion

           int index(0);
           scalar uPrev=0.0; vector CPrev(0., 0., 0.);
           vector refDirU(refDir/mag(refDir));
      
           forAll(cPatchvort, facei )       
              {
                label  curCell = cPatchvort.faceCells()[facei];
                scalar uCmp = (U()[curCell] & refDirU);
               
                if (uPrev*uCmp<0.0)
                 {
                   vector r_curCell = -uCmp * ( CPrev - C[curCell] ) / (uPrev - uCmp) + C[curCell];
  
                   // col N: Time
                   outS1_() << mag( ( (r_curCell - refPoint) & refDirU ) ) << tab; 
                   
                   index++;
                 } 
       
                uPrev = uCmp;
                CPrev = C[curCell];
            } 
            
           outS1_() << endl;
           
         }
         
// Xl_top
        {
          const polyPatch& cPatchvort = U().mesh().boundaryMesh()[lipTop_];
        
          // Define reference parameters 
 
           vector refPoint(0., 1., 0.5); // Reference point to zero the vortex length
           vector refDir(0., 1., 0.); // Vector aligned with the wall

          // Compute vortex length based on the point of velocity inversion

           int index(0);
           scalar uPrev=0.0; vector CPrev(0., 0., 0.);
           vector refDirU(refDir/mag(refDir));
      
           forAll(cPatchvort, facei )       
              {
                label  curCell = cPatchvort.faceCells()[facei];
                scalar uCmp = (U()[curCell] & refDirU);
               
                if (uPrev*uCmp<0.0)
                 {
                   vector r_curCell = -uCmp * ( CPrev - C[curCell] ) / (uPrev - uCmp) + C[curCell];
  
                   // col N: Time
                   outS2_() << mag( ( (r_curCell - refPoint) & refDirU ) ) << tab; 
                   
                   index++;
                 } 
       
                uPrev = uCmp;
                CPrev = C[curCell];
            } 
            
           outS2_() << endl;
           
         }
         
// Xl_down
        {
          const polyPatch& cPatchvort = U().mesh().boundaryMesh()[lipDown_];
        
          // Define reference parameters 
 
           vector refPoint(0., -1., 0.5); // Reference point to zero the vortex length
           vector refDir(0., 1., 0.); // Vector aligned with the wall

          // Compute vortex length based on the point of velocity inversion

           int index(0);
           scalar uPrev=0.0; vector CPrev(0., 0., 0.);
           vector refDirU(refDir/mag(refDir));
      
           forAll(cPatchvort, facei )       
              {
                label  curCell = cPatchvort.faceCells()[facei];
                scalar uCmp = (U()[curCell] & refDirU);
               
                if (uPrev*uCmp<0.0)
                 {
                   vector r_curCell = -uCmp * ( CPrev - C[curCell] ) / (uPrev - uCmp) + C[curCell];
  
                   // col N: Time
                   outS3_() << mag( ( (r_curCell - refPoint) & refDirU ) ) << tab; 
                   
                   index++;
                 } 
       
                uPrev = uCmp;
                CPrev = C[curCell];
            } 
            
           outS3_() << endl;
           
         }
         
//****  User-defined function ENDS here *****//
    
      counter_= 1;
     
    }// if counter
    
 }// if enabled 
       
} 


// ************************************************************************* //
