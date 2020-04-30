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
  
  #include "setRootCase.H"
  #include "createTime.H"
    
  //- Loop over all groups 
    
  fileName ppDirRoot(runTime.path()/"rheoToolPP"/argv[1]/"moleculesStats");    
  fileNameList groups(readDir(ppDirRoot, fileType::directory));
   
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
    fileName eTeDir(cwdDir/"Stretch.txt");
    
    IFstream eTeFile(eTeDir);
    
    if (!eTeFile.opened())
     {
        FatalErrorIn("averageMolcN")
         << "\nFile "<< eTeDir <<" not found.\n"
         << exit(FatalError);
     }    
    
    string leTe;
    string lx;
    eTeFile.getLine(leTe);    
    IStringStream ieTe(leTe);
    
    // Get the number of cols (molecules)
    int nCols(0);
    double time(0);
    while(ieTe >> time)
    {
        nCols++;        
    }
    
    int nMolc(nCols-1); // First col is time
         
    // Back to first line
    eTeFile.rewind();
    
    Info << "\nAverage for group " << groups[gi] << ", starting from t = " << argv[1] << nl  << endl;
    
    // Write stream
    OFstream oFile(cwdDir/"Stretch_Naverage.txt");
    
    // Loop over all lines 
    scalar tmp(0.); 
    while (eTeFile.getLine(leTe))
    {
       IStringStream ieTe(leTe);
       ieTe >> time; // This is the time, ignore
       oFile << time << tab;
      
       // Loop over all columns (molecules) 
       int cnt(0);   
       scalar avSum(0.);
       for (int i=0; i<nMolc; i++)
       {
         ieTe >> tmp;
                      
         // Untracked molecules have molcLen = 0
         if (tmp>0)    
         {
           avSum += tmp;
           cnt++;    
         }                  
       } 
       
       if (cnt>0)
       { 
         oFile << avSum/cnt << endl; 
       }
       else
       {
         oFile << 0. << endl; 
       }        
   }
   
  }    
 
  Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << tab 
       << "ClockTime = " << runTime.elapsedClockTime() << " s"
       << nl << endl;
       
  Info << "Done.\n" << endl;
   
  return 0;
}


// ************************************************************************* //
