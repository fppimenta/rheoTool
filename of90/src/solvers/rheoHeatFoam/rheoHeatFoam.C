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
    rheoEFoam

Description
    Transient solver for incompressible, laminar electrically-driven flow. Pressure
    forcing can also co-exist. Any GNF or viscoelastic model of library
    lconstitutiveEquations can be selected. Both coupled and
    segregated solvers can be runtime selected. In the segregated solver,
    pressure-velocity coupling is ensured through SIMPLEC algorithm.
    Non-isothermal solver.
    
    This solver is part of rheoTool.

\*---------------------------------------------------------------------------*/

#include "sparseMatrixSolvers.H" // Avoid namespace Foam

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "adjustCorrPhi.H"
#include "ppUtilInterface.H"
#include "constitutiveModel.H"
#include "EDFModel.H"
#include "fluidThermoModel.H"

#include "blockOperators.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createControls.H" 
    #include "createPPutil.H"
    #include "createUfIfPresent.H"
    
    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        #include "readDyMControls.H"

        if (solveFluid)
        {
          if (LTS)
          {
             #include "setRDeltaT.H"
          }
          else
          {
             #include "CourantNo.H"
             #include "setDeltaT.H"
          }
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Inner loop iterations ---
        for (int i=0; i<nInIter; i++)
        {
            Info << "Inner iteration:  " << i << nl << endl; 
            
            if (i==0 || moveMeshOuterCorrectors)
            {
                fvModels.preUpdateMesh();
                
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        #include "correctPhi.H"
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }
            
            fvModels.correct(); 
            
            if (solveCoupled)
            {            
              // Add/Solve electric equations
              elecM.correct();
                
              if (solveFluid)
              {          
                // When we enter here we know that there is no gravity,
                // because this has been checked earlier. phig = 0 and p = p_rgh
        
                #include "pUEqn.H" 
                
                // Add/solve constitutive equation 
                constEq.correct();      
                
                // Solve all coupled
                cps->solve(); 
                
                phi = fvc::flux(U) + pRC - fvc::snGrad(p_rgh)*fvc::interpolate(rAU)*mesh.magSf(); 
                      
                #include "continuityErrs.H"
                
                fvConstraints.constrain(U);
                
                // Correct Uf if the mesh is moving
                fvc::correctUf(Uf, U, phi);
                                
                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi, U); 
                
                // Update p
                p == p_rgh;

                if (p_rgh.needReference())
                {
                  p += dimensionedScalar
                  (
                    "p",
                    p.dimensions(),
                    pressureReference.refValue()
                  - getRefCellValue(p, pressureReference.refCell())
                  );
                  p_rgh = p;
                }            
              }             
            } 
            else
            {
              if (solveFluid)
              {
                #include "UEqn.H"
                
                for (int nci=0; nci<nCorrectors; nci++)
                {
                  #include "pEqn.H"
                }
                
                // ---- Solve constitutive equation ----	
                constEq.correct();
              }
              
              // ---- Update electric terms ----
              elecM.correct();
            }  
            
            // ---- Update thermo ---
            // Can be solved coupled/segregated independent of p-U
            thermo->correct(U,phi,constEq.tauTotal(),cpsT,fvModels,fvConstraints,simple.nCorrNonOrth());  
        
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
