/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "pimpleControl.H"
#include "IFstream.H"
#include "OFstream.H"

#include "ppUtilInterface.H"
#include "constitutiveTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pimpleControl pimple(mesh);

#   include "readGravitationalAcceleration.H"
#   include "initContinuityErrs.H"
#   include "createFields.H"
#   include "createPPutil.H"
#   include "createTimeControls.H"
#   include "correctPhi.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Inner loop iterations ---

        for (int i=0; i<nInIter; i++)
	  {
            twoPhaseProperties.correct();

#           include "alphaEqnSubCycle.H"

#           include "UEqn.H"

            // --- PISO loop
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

#           include "continuityErrs.H"

            p = pd + rho*gh;

            if (pd.needReference())
            {
                p += dimensionedScalar
                (
                    "p",
                    p.dimensions(),
                    pRefValue - getRefCellValue(p, pdRefCell)
                );
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
