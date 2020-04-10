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

#include "BarycentricWeights.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvcGrad.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace externalForcingInterps
{
    defineTypeNameAndDebug(BarycentricWeights, 0);
    addToRunTimeSelectionTable(externalForcingInterp, BarycentricWeights, dicteFI);
}
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalForcingInterps::BarycentricWeights::updateUex()
{
   const volVectorField& U = mesh_.lookupObject<const volVectorField>("U");
   
   Uex_ = U;
 
   if (epActive_)
    {    
     // Add electrophoretic velocity 
     const volScalarField& phiE = mesh_.lookupObject<const volScalarField>("phiE");
         
     dimensionedScalar epMobilityDim("epMobilityDim", dimensionSet(-1, 0, 2, 0, 0, 1, 0), epMobility_);
     
     Uex_ -= fvc::grad(phiE)*epMobilityDim;         
    } 
}

void Foam::externalForcingInterps::BarycentricWeights::update()
{

 if(!frozenFlow_)
 {  
   updateUex();
   interpolatorCPB_.updatePsip(Uex_);  
 }
  
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalForcingInterps::BarycentricWeights::BarycentricWeights
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
externalForcingInterp(mesh, dict),
epActive_(readBool(dict.subDict("electrophoresis").lookup("active"))),
epMobility_(epActive_ ? readScalar(dict.subDict("electrophoresis").lookup("mobility")) : 0.),
frozenFlow_(readBool(dict.subDict("externalFlow").lookup("frozenFlow"))),
Uex_(mesh.lookupObject<const volVectorField>("U")),
interpolatorCPB_(Uex_)
{
  // Correct Uex if EP
  updateUex();
  interpolatorCPB_.updatePsip(Uex_); 
}

// ************************************************************************* //
