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
    

Description
     

\*---------------------------------------------------------------------------*/

#include "sparseMatrixSolvers.H" // Avoid namespace Foam

#include "fvCFD.H"
#include "regionProperties.H"
#include "iCourantNo.H"
#include "fvOptions.H"
#include "pimpleMultiRegionControl.H"

#include "adjustCorrPhi.H"      
#include "ppUtilInterface.H"
#include "constitutiveModel.H"
#include "EDFModel.H"
#include "fluidThermoModel.H"
#include "solidThermoModel.H"

#include "coupledTFvPatchScalarField.H"

#include "blockOperators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "multiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"           
    #include "setInitialMultiRegionDeltaT.H"       

    while (pimples.run(runTime))
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"

        #include "multiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"            
        #include "setMultiRegionDeltaT.H"              

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PIMPLE loop
        while (pimples.loop())
        {
            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }
            
            #include "coupledTEqn.H"
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
