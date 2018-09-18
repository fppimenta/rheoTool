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
    Transient solver for incompressible, laminar electrically-driven flow. Pressure
    forcing can also co-exist. Any GNF or viscoelastic model of library
    lconstitutiveEquations can be selected. Pressure-velocity coupling is using 
    the SIMPLEC algorithm.
    
    This solver is part of rheoTool.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "dynamicFvMesh.H"
#include "CorrectPhi.H"

#include "adjustCorrPhi.H" 

#include "ppUtilInterface.H"
#include "constitutiveModel.H"
#include "EDFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
   
int main(int argc, char *argv[])
{
    #include "postProcess.H"
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMeshDict.H"  
    #include "createDynamicFvMesh.H"  
    #include "initContinuityErrs.H" 
    
    #include "createFields.H"  
    #include "createControls.H" 
    #include "createUfIfNeeded.H"
    #include "createFvOptions.H"
    #include "createPPutil.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
 
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
 
    // --- Time loop ---

    while (simple.loop())
    {
        #include "readControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        
        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.update();

        if (mesh.changing())
         {
           // Calculate absolute flux from the mapped surface velocity
           phi = mesh.Sf() & Uf();
        
           if (correctPhi)
           {
              #include "correctPhi.H"
           }

           // Make the flux relative to the mesh motion
           fvc::makeRelative(phi, U);

           if (checkMeshCourantNo)
           {
             #include "meshCourantNo.H"
           }
         }
         
        // --- Inner loop iterations ---
        for (int i=0; i<nInIter; i++)
	  {

            Info << "Inner iteration:  " << i << nl << endl; 
           
            if (solveFluid)
            {
               // --- Pressure-velocity SIMPLEC corrector
               {
                 // ---- Solve U and p ----	
                 #include "UEqn.H"
                 #include "pEqn.H"         
               }
            
               // ---- Solve constitutive equation ----	
               constEq.correct();
            }
             
            // ---- Update electric terms ----
            elecM.correct();

            // --- Passive Scalar transport
            if (sPS)
             {
               #include "CEqn.H"
             }

         }

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
