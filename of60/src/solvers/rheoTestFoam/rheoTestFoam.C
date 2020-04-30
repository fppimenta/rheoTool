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

Application
    rheoTestFoam

Description
    Testing application to be used with library lconstitutiveEquations. This 
    application returns the principal components of the extra-stress tensor,
    given the velocity gradient tensor, which is defined by the user. Two modes
    are available: ramp (several shear-rates evaluated up to steady-state) and
    non-ramp (single shear-rate evaluated over time). The solver manipulates
    internally the BC in order to get the given shear-rate. The mesh to be used
    with this solver is a one-cell unitary cube.
    
    This solver is part of rheoTool.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "simpleControl.H"

#include "constitutiveModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
 
    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        for (int i=0; i<nInIter; i++)
        {
           Info<< "Iteration  " << i << nl << endl; 
            
           // --- Solve only the constitutive equation           
           constEq.correct();
        }
        
        extraStress = constEq.tauTotal();

//      runTime.write(); // Uncomment if needed
      
        if (ramp_)
        {
           scalar relError = Foam::mag(extraStress[0]-oldExtraStress)/(Foam::mag(extraStress[0]) + SMALL); 
           
           if (relError < 1e-8 || cnt > 5000)
           {
             #include "reStartOrEnd.H"
           }
               
           oldExtraStress = extraStress[0];
           cnt++;
        }
        else
        {
           #include "timeWrite.H"
        }
         
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    
    res << endl;

    return 0;
}


// ************************************************************************* //
