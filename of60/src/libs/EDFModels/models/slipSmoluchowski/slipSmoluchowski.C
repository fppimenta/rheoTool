/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "slipSmoluchowski.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace EDFEquations
{
    defineTypeNameAndDebug(slipSmoluchowski, 0);
    addToRunTimeSelectionTable(EDFEquation, slipSmoluchowski, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquations::slipSmoluchowski::slipSmoluchowski
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    EDFEquation(name, phi),
    phiE_
    (
        IOobject
        (
            "phiE" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    zeroFe_
    (
        IOobject
        (
            "zeroForce",
            phi.time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi.mesh(),
        dimensionedVector
        (
                "zero",
                dimVelocity*dimDensity/dimTime,
                pTraits<vector>::zero
        ),
        extrapolatedCalculatedFvPatchField<vector>::typeName
    ),
    phiEEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    maxIterPhiE_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<int>("maxIter", 50))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::EDFEquations::slipSmoluchowski::Fe() const
{
    return
    (
         zeroFe_ 
    );     
}

void Foam::EDFEquations::slipSmoluchowski::correct()
{

       scalar res=GREAT; 
       scalar iter=0;  
  
   //- Equation for the external potential (loop for the case
   //  of non-orthogonal grids)   
       while (res > phiEEqRes_ && iter < maxIterPhiE_)
          { 
          
		fvScalarMatrix phiEEqn
		(
		    fvm::laplacian(phiE_)
		);
		
		phiEEqn.relax();
		res=phiEEqn.solve().initialResidual();
		
		iter++;
          } 
}

// ************************************************************************* //
