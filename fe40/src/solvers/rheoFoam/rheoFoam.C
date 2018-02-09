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
    rheoFoam

Description
    Transient solver for incompressible, laminar flow. Any GNF or viscoelastic
    model of library lconstitutiveEquations can be selected. Pressure-velocity
    coupling is using the SIMPLEC algorithm.
    
    This solver is part of rheoTool.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "simpleControl.H"

#include "ppUtilInterface.H"
#include "constitutiveModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    simpleControl simple(mesh);
     
    #include "createFields.H" 
    #include "createPPutil.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    // --- Time loop ---

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Inner loop iterations ---

        for (int i=0; i<nInIter; i++)
	  {

            Info<< "Inner iteration:  " << i << nl << endl; 
           
            // --- Pressure-velocity SIMPLEC corrector
            {
               // ---- Solve constitutive equation ----	
               constEq.correct(); 

               // ---- Solve U and p ----	
               #include "UEqn.H"
               #include "pEqn.H"
               
            }

            // --- Passive Scalar transport
            if (sPS)
             {
               #include "CEqn.H"
             }
         }
         
        //- Post-processing               
        postProc.update();
      
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
