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
    rheoFoam

Description
    Transient solver for incompressible, laminar flow. Any GNF or viscoelastic
    model of library lconstitutiveEquations can be selected. Both coupled and
    segregated solvers can be runtime selected. In the segregated solver,
    pressure-velocity coupling is ensured through SIMPLEC algorithm.
    
    This solver is part of rheoTool.

\*---------------------------------------------------------------------------*/

#include "sparseMatrixSolvers.H" // Avoid namespace Foam

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "simpleControl.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "adjustCorrPhi.H"
#include "ppUtilInterface.H"
#include "filmModel.H"
#include "fluidThermoModel.H"

#include "blockOperators.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
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

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Inner loop iterations ---
        for (int i=0; i<nInIter; i++)
        {
            Info << "Inner iteration:  " << i << nl << endl; 
            
            if (i==0 || moveMeshOuterCorrectors)
            {
              mesh.update();

              if (checkMeshCourantNo)
              {
               #include "meshCourantNo.H"
              }
            }     
            
            // ---- Solve Momentum equation                
            fvVectorMatrix UEqn
            (
               fvm::ddt(film.rho()*film.h(), U)
             ==   
               film.divTau(U)
             + g_*film.rho()*film.h()  
            );                      

            UEqn.relax();
            spSolverU->solve(UEqn);  
                     
            // Flux (not conservative in xy)
            phi = (fvc::interpolate(U) & mesh.Sf()); 

            // Correct Uf  
            fvc::correctUf(Uf, U, phi);
            
            // Update the free surface pts before next mesh motion 
            if (film.absoluteFluxNeeded())
            {
              //- Streamline method to update free-surface (phi is absolute)
              film.updateFreeSurface(i==nInIter-1, U, phi);   
              fvc::makeRelative(phi, U);
            }
            else
            {
              //- Peric's method to update free-surface (phi is relative)              
              fvc::makeRelative(phi, U);
              film.updateFreeSurface(i==nInIter-1, U, phi);               
            }
            
            // Check the flux on freeSurfaces (relative flux should be close to zero)
            film.checkFlux(U, phi);  

            // Force the zero-flux condition on free-surface
            film.forcedZeroFluxFreeSurface(true, phi);

            // ---- Solve Height equation  
            film.updateHeight(phi);             
                
            // ---- Solve constitutive equation  
            film.correctStresses(U); 
            
            // ---- Solve for temperature if non-isothermal
            film.correctThermo(U, phi); 
        }
         
        postProc.update();
        
        // Film casting specific post-processing
        film.postProcess(U);
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
