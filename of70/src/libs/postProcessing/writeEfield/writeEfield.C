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
    writeEfield

Description
    Reads the electric potential fields (all available), adds them, 
    computes the gradient, and writes the results as the electric field. 
    
    This utility is part of rheoTool.
   
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "EDFModel.H"
#include "timeSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
int main(int argc, char *argv[])
{
    Foam::timeSelector::addOptions();
    
    #include "setRootCase.H"
    #include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);
    #include "createMesh.H"
    #include "createFields.H"
       
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
    
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
	
	E *= 0;
	
	 IOobject psiHeader
    	 (
	       "psi",
	       runTime.timeName(),
	       mesh,
	       IOobject::MUST_READ
	 );
	 
	 bool psiHeaderOk(psiHeader.typeHeaderOk<volScalarField>(true));
	 
	 if (psiHeaderOk)
         {         
           Info<< "Reading field psi\n" << endl;
           volScalarField psi_(psiHeader, mesh); 
           
           E = -fvc::grad(psi_);          
         }
	 
	 IOobject phiEHeader
    	 (
	       "phiE",
	       runTime.timeName(),
	       mesh,
	       IOobject::MUST_READ
	 );
	 
	 bool phiEHeaderOk(phiEHeader.typeHeaderOk<volScalarField>(true));
    	
    	 if (phiEHeaderOk)
         {         
           Info<< "Reading field phiE\n" << endl;
           volScalarField phiE_(phiEHeader, mesh);
           
           E -= fvc::grad(phiE_);           
         }
         
         if ( !psiHeaderOk && !phiEHeaderOk )
         { 
            FatalErrorIn("writeEfield.C")
            << "No electric potential field (phiE or psi) was found."
            << exit(FatalError);
         }
	
	 Info<< "Writing E\n" << endl;
	 E.write();
    }
    
    Info<< "End\n" << endl;

    return 0;
}


