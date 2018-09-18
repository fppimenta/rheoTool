/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "solidParticleCloud.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticleCloud::solidParticleCloud
(
    const fvMesh& mesh,
    const word& cloudName,
    bool readFields
)
:
    Cloud<solidParticle>(mesh, cloudName, false),
    mesh_(mesh),
    molcProperties_
    (
        IOobject
        (
            "moleculesControls",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    molcToDelete_(1,-1),
    wallRepX_(readScalar(molcProperties_.subDict("exclusionVolumeProperties").lookup("repulsiveDistance"))),
    isTethered_(readBool(molcProperties_.subDict("externalFlow").lookup("tethered"))),
    extForcInt_(externalForcingInterp::New(mesh, molcProperties_))
{
    if (readFields)
    {
        solidParticle::readFields(*this);
    }
       
    // Force repulsion to be positive (inward)
    wallRepX_ = Foam::mag(wallRepX_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::List<Foam::label>& Foam::solidParticleCloud::move(bool includeDrag)
{
 
 // Empty the list of molc to remove    
 molcToDelete_.clear();
  
 // May need interpolation of external forcing     
 if (includeDrag)
  {  
    extForcInt_->update();
  }   
 
 solidParticle::trackingData  td(*this, includeDrag);
 
 // Move the particles
 Cloud<solidParticle>::move(*this, td, mesh_.time().deltaTValue());
 
 return molcToDelete_;
 
}

// ************************************************************************* //
