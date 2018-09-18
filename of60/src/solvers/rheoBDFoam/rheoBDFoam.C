/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    rheoBDSFoam

Description
    Transient Brownian dynamics solver. Continuum forcing can be analytical
    or numerical. If numerical, it can be pressure-driven flow, electrically-driven
    flow and/or electrophoresis.
    
    Solver structure derived from rheoEFoam. 
    
    This solver is part of rheoTool.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "simpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "subCycle.H"
#include "adjustCorrPhi.H"
#include "ppUtilInterface.H"
#include "constitutiveModel.H"
#include "EDFModel.H"

#include "sPCloudInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    // Prevent parallel runs
    argList::noParallel();
    
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createLagrangianFields.H"
    #include "createControls.H" 
    #include "createPPutil.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
     
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        #include "readDyMControls.H"
        if (solveFluid)
         {
           #include "CourantNo.H"        
           #include "setDeltaT.H"
         }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Inner loop iterations ---
        for (int i=0; i<nInIter; i++)
        {
            if (solveFluid)
              Info << "Inner iteration:  " << i << nl << endl; 
            
            if (i==0 || moveMeshOuterCorrectors)
            {
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // --- Pressure-velocity SIMPLEC corrector
            if (solveFluid)
            {
              {
                 // ---- Solve U and p ----	
                 #include "UEqn.H"
                 #include "pEqn.H"         
              }
            
              // ---- Solve constitutive equation ----	
              constEq.correct();
            }
            
            // ---- Update electric terms ----
            if (solveElecM)
            {
              elecM.correct();
            }

            // --- Passive Scalar transport
            if (sPS)
             {
               #include "CEqn.H"
             }            
        }
        
        #include "moleculesEqns.H"

        postProc.update();
        runTime.write();
        
        // This only controls writing of continuous fields, not lagrangian ones
        if (writeContFields_)
         {
           runTime.write();
         }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
