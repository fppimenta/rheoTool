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

#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"

#include "constitutiveTwoPhaseMixture.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // --- Time loop ---

    while (runTime.run())
    {
       
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
  
        #include "alphaEqnSubCycle.H"            
        interface.correct();                     

        // --- Inner loop iterations ---

        for (int i=0; i<nInIter; i++)
	  {

            Info<< "Inner iteration:  " << i << nl << endl; 
           
            // --- Pressure-velocity SIMPLEC corrector
            {
               // ---- Solve the constitutive equation ----	
               twoPhaseProperties.correct();            
                            
               // ---- Solve U and p ----
               
               #include "UEqn.H"
               
               if (simplec)
                {	
                  #include "pEqn.H"
                }
               else
                {	
                  while (pimple.correct())
                   {
                      #include "pEqnP.H"
                   }   
                }               
            }

            // --- Passive Scalar transport
            if (sPS)
             {
               #include "CEqn.H"
             }

         }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
