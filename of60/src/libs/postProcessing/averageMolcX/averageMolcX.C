/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "IFstream.H"
#include "OFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  argList::noParallel();
  argList::validArgs.append("time");
  argList::addOption
    (
        "biased",
        "bool",
        "use biased (true) vs unbiased (false-default) average"
    );
    
  argList::addOption
    (
        "startPoint",
        "vector",
        "e.g. \"(1 0 0)\". Start point of the line"
    );
  
  argList::addOption
    (
        "endPoint",
        "vector",
        "e.g. \"(5 0 0)\". End point of the line"
    );
    
  argList::addOption
    (
        "nBins",
        "int",
        "Number of bins over the line"
    );
    
  
  #include "setRootCase.H"
  #include "createTime.H"

  // Read/parse
  bool biased(false); bool inpMeth = args.optionReadIfPresent("biased", biased);
  vector p1(0,0,0); inpMeth = args.optionReadIfPresent("startPoint", p1);
  vector p2(0,0,0); inpMeth = args.optionReadIfPresent("endPoint", p2);
  int nBins(1); inpMeth = args.optionReadIfPresent("nBins", nBins);
  
  if (nBins < 1)
   {
     FatalErrorInFunction
     << "The number of bins is less than 1. Please provide at least one bin."
     << exit(FatalError);
   }
  else
   {
     if (mag(p1-p2)>1e-20)
      {
        Info << nl << "Starting the average over the line " << p1 << " -> " << p2 << " using " << nBins << " bins." << endl;
      }
     else
      {
         FatalErrorInFunction
         << "The start " << p1 << " and end point " << p2 << " of the line seem to coincide or they are too close.\n"
         << "Please check the entries -startPoint and -endPoint (use -help to check their meaning)."
         << exit(FatalError);
      }    
   }

  // Those 2 are only used for the unbiased method
  autoPtr< List < List<scalar > > > allPosptr;  
  autoPtr< List < List<scalar > > > allStretchptr;
     
  // Build container: xmid | ymid | zmid | StretchAverg | nHits
  List<List<scalar> > M(nBins, List<scalar>(5,0.));
  
  vector dr((p2-p1)/nBins);
  scalar magdr(mag(dr));
  vector ndr(dr/magdr);
  vector p0(p1+0.5*dr);
    
  for(int i=0; i<nBins; i++)
   {
     vector pi = p0 + dr*i;
     M[i][0] = pi.x();
     M[i][1] = pi.y();
     M[i][2] = pi.z(); 
   }
    
  //- Loop over all groups 
    
  fileName ppDirRoot(runTime.path()/"rheoToolPP"/argv[1]/"moleculesStats");    
  fileNameList groups(readDir(ppDirRoot, fileName::DIRECTORY));
   
  if (groups.size()<1)
   {
      FatalErrorIn("averageMolcX")
      << "\nDirectory:  \n\n"<< ppDirRoot << nl << nl 
      <<"was not found or is empty. Check if time = " << argv[1] << " is valid.\n"
      << exit(FatalError);
   }
    
  forAll(groups, gi)
  {
    fileName cwdDir(ppDirRoot/groups[gi]);
    
    // Read data
    fileName xDir(cwdDir/"X.txt");
    fileName eTeDir(cwdDir/"Stretch.txt");
    
    IFstream xFile(xDir);
    IFstream eTeFile(eTeDir);
    
    if (!xFile.opened())
     {
        FatalErrorIn("averageMolcX")
         << "\nFile "<< xDir <<" not found.\n"
         << exit(FatalError);
     } 
     
    if (!eTeFile.opened())
     {
        FatalErrorIn("averageMolcX")
         << "\nFile "<< eTeDir <<" not found.\n"
         << exit(FatalError);
     }    
    
    string leTe;
    string lx;
    eTeFile.getLine(leTe);    
    IStringStream ieTe(leTe);
    
    // Get the number of cols (molecules)
    int nCols(0);
    double tmp(0);
    while(ieTe >> tmp)
    {
        nCols++;        
    }
    
    int nMolc(nCols-1); // First col is time
    
    int ndt(1); // Number of rows (times)
    if (!biased)
     {
        while (eTeFile.getLine(leTe))
         ndt++;

        // Clean EOF flag
        eTeFile.stdStream().clear();

        allPosptr = new List<List<scalar > >(ndt, List<scalar>(3*nMolc, 0.) );    
        allStretchptr = new List<List<scalar > >(ndt, List<scalar>(nMolc, 0.) );    
     }
    
    // Back to first line
    eTeFile.rewind();
    
    if (biased)
    {
      Info << "\nBiased average for group " << groups[gi] << ", starting from t = " << argv[1] << nl  << endl;
    }
    else
    {
      Info << "\nUnbiased average for group " << groups[gi] << ", starting from t = " << argv[1] << nl  << endl;
    }
     
    // Loop over all lines
    scalar x(0.), y(0.), z(0.);
    scalar molcLen(0.);
    int ti = 0;
    
    while (eTeFile.getLine(leTe))
    {
       IStringStream ieTe(leTe);
       ieTe >> tmp; // This is the time, ignore
       
       xFile.getLine(lx);
       IStringStream ix(lx);
       ix >> tmp; // This is the time, ignore
        
       if (biased)
       {
         // Loop over all columns (molecules)    
         for (int i=0; i<nMolc; i++)
          {
            ix >> x >> y >> z;
            ieTe >> molcLen;
                      
            // Untracked molecules have molcLen = 0
            if (molcLen>0)    
            {
              vector v(vector(x,y,z) - p1);
              scalar vtmag(ndr & v);
               
              // if vtmag < 0 the projected position of the particle
              // does not lie between p1 and p2, thus ignore 
              if (vtmag > 0)
               {
                  int bin = int(vtmag/magdr);
                
                  if (bin<nBins)
                   {
                     M[bin][3] = ( M[bin][3]  * M[bin][4]  + molcLen ) / (M[bin][4]  + 1.);
                     M[bin][4] ++;    
                   }
               }
             
             }                  
         }    
       }
       else
       { 
         // Collect all data into matrices
         int jp = 0;
         int js = 0;
 
         // Loop over all columns (molecules)    
         for (int i=0; i<nMolc; i++)
         {
            ix >> allPosptr()[ti][jp] >> allPosptr()[ti][jp+1] >> allPosptr()[ti][jp+2];
            ieTe >> allStretchptr()[ti][js];   
            jp += 3;
            js++;                
         }    
       }
     
       ti++;
    }
    
    // Processing for unbiased molecule average
    if (!biased)
     {
       // Loop over molecules
       for (int mi=0; mi<nMolc; mi++)
        {
           // Average stretch in each bin | nHits
           List<List<scalar> > Mtmp(nBins, List<scalar>(2,0.)); 
           
           for (int ti=0; ti<ndt; ti++)
            { 
              // Untracked molecules have molcLen = 0
              if (allStretchptr()[ti][mi]>0)    
               {
                 vector v( vector(allPosptr()[ti][3*mi], allPosptr()[ti][3*mi+1], allPosptr()[ti][3*mi+2]) - p1 );
                 scalar vtmag(ndr & v);
              
                 // if vtmag < 0 the projected position of the particle
                 // does not lie between p1 and p2, thus ignore 
                 if (vtmag > 0)
                  {                  
                    int bin = int(vtmag/magdr);
                    
                    if (bin<nBins)
                    {
                      Mtmp[bin][0] = (Mtmp[bin][0]*Mtmp[bin][1] + allStretchptr()[ti][mi])/( Mtmp[bin][1] + 1 );
                      Mtmp[bin][1]++;    
                    }
                  }           
               }                             
            }
            
           // Transfer the molecule average contribute to the main container.
           // Increment the counter of molecules by 1 (all molecules contribute equally->unbiased) or 0.
           for (int i = 0; i<nBins; i++)
            { 
              M[i][3] += Mtmp[i][0];
              M[i][4] += Foam::min(1, Mtmp[i][1]);   
            }
                         
        } 
      
        // After collecting all molecules divide by the nb of molecules that contributed to each bin
        for (int i = 0; i<nBins; i++)
         { 
           if (M[i][4]>0)
             M[i][3] /= M[i][4];   
         }
         
        allPosptr.clear();    
        allStretchptr.clear();          
     }
     
    // Write statistics (no need to check for the path, since this was done before)
    OFstream oFile(cwdDir/"Stretch_Xaverage.txt");
    
    for(int i=0; i<nBins; i++)
     {
       oFile << M[i][0] << tab << M[i][1] << tab << M[i][2] << tab << M[i][3] << tab << M[i][4] << nl; 
       
       // Reset to 0 averageStretch M[i][3] and counter M[i][4] for new group
       M[i][3] = 0;
       M[i][4] = 0;
     }
    oFile << endl;
  }    
  
  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << tab 
       << "ClockTime = " << runTime.elapsedClockTime() << " s"
       << nl << endl;
       
  Info << "Done.\n" << endl;
   
  return 0;
}


// ************************************************************************* //
