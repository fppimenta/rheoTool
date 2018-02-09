/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    rheoInterFoam

Description
    Transient solver for 2 incompressible isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach. The flow is
    laminar. Any GNF or viscoelastic model of library lconstitutiveEquations can be
    selected. The pressure-velocity coupling is runtime selectable: either SIMPLEC
    or PIMPLE (the momentum equation is always solved in both cases).
    
    Based on solver interFoam of foam-extend.
    
    This solver is part of rheoTool.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "immiscibleConstitutiveTwoPhaseMixture.H"
#include "ppUtilInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createPPutil.H"
    #include "correctPhi.H"

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    // Read extra-controls

    int    nInIter = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<int>("nInIter", 1);
    bool   simplec = mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<Switch>("SIMPLEC", true);
    bool   sPS = cttProperties.subDict("passiveScalarProperties").lookupOrDefault<Switch>("solvePassiveScalar", false);
    if (sPS) C.writeOpt() = IOobject::AUTO_WRITE;
    

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H" 
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
       for (int i=0; i<nInIter; i++)
        {
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correctAll();

            #include "UEqn.H"

            if (simplec)
             {	
               while (pimple.correct())
                {
                   #include "pEqn.H"
                } 
             }
            else
             {	
               while (pimple.correct())
                {
                   #include "pEqnP.H"
                }   
             }  
             
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
