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
   rndPerturbation   

Description
   Local, random perturbation of fields cCation and cAnion. Names are 
   hard-coded. This utility is used in tutorial selecMembrane.  

   This utility is part of rheoTool.
   

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "EDFModel.H"
#include "Random.H"
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
 
    Random rndSeed(1234567); 
    
    forAll(timeDirs, timeI)
    {
     runTime.setTime(timeDirs[timeI], timeI);
        
     Info<< "Time = " << runTime.timeName() << nl << endl;
        
     IOobject catHeader
      (
        "cCation",
        runTime.timeName(),
	mesh,
	IOobject::MUST_READ
      );
      
      IOobject catLogHeader
      (
        "cCationLog",
        runTime.timeName(),
	mesh,
	IOobject::MUST_READ
      );
 
     if (catHeader.headerOk())
      {         
           
        Info<< "Reading field cCation\n" << endl;
	volScalarField cat
	(
	    IOobject
	    (
	        "cCation",
	        runTime.timeName(),
	        mesh,
	        IOobject::MUST_READ
	    ),
	    mesh
	);
	
	Info<< "Reading field cAnion\n" << endl;
	volScalarField ani
	(
	    IOobject
	    (
	        "cAnion",
	        runTime.timeName(),
	        mesh,
	        IOobject::MUST_READ
	    ),
	    mesh
	);

	forAll(ani, i)
	 {
	    cat[i] *= (1. + .01*( -1. + 2.*rndSeed.scalar01() ) ); 
	    ani[i] *= (1. + .01*( -1. + 2.*rndSeed.scalar01() ) );
	 }
	 
	cat.correctBoundaryConditions();
	ani.correctBoundaryConditions();
	
	cat.write();
	ani.write();
      }
     else if (catLogHeader.headerOk())
      {
        Info<< "Reading field cCation\n" << endl;
	volScalarField cat
	(
	    IOobject
	    (
	        "cCationLog",
	        runTime.timeName(),
	        mesh,
	        IOobject::MUST_READ
	    ),
	    mesh
	);
	
	Info<< "Reading field cAnion\n" << endl;
	volScalarField ani
	(
	    IOobject
	    (
	        "cAnionLog",
	        runTime.timeName(),
	        mesh,
	        IOobject::MUST_READ
	    ),
	    mesh
	);

	forAll(ani, i)
	 {	 
	   cat[i] += Foam::log( (1. + .01*( -1. + 2.*rndSeed.scalar01() ) ) ); 
	   ani[i] += Foam::log( (1. + .01*( -1. + 2.*rndSeed.scalar01() ) ) ); 
	 }
	 
	cat.correctBoundaryConditions();
	ani.correctBoundaryConditions();
	
	cat.write();
	ani.write();
      
      }
     else
      {
        Info << "Error" << endl;
      } 
     
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
