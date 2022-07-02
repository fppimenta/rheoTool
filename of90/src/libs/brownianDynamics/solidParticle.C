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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidParticle, 0);
}

const Foam::scalar Foam::solidParticle::PREC_(1.-1e-7);
const int Foam::solidParticle::NMAXITER_(250);
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidParticle::move
(
    solidParticleCloud& cloud,
    trackingData& td,
    const scalar trackTime
)
{

// If molecule is tethered and if we get 
// the first bead, skip tracking
if ( cloud.isTethered() )
 {
  if ( ids_[1] == 0 )
   {
     td.switchProcessor = false;
     td.keepParticle = true;
     return td.keepParticle;
   }
 }

td.switchProcessor = false;
td.keepParticle = true;

const polyBoundaryMesh& pbMesh = mesh().boundaryMesh();

if (td.includeDrag())
 {     
    vector UBEV(U_); // Velocity from Brownian motion and Exclusion Volume 
    bool onWall(false);
    
    int cnt(0);
    while (td.keepParticle && !td.switchProcessor && stepFraction() < PREC_ && cnt < NMAXITER_)
    {
        if (debug)
        {
            Info<< "Time = " << mesh().time().timeName()
                << " trackTime = " << trackTime
                << " steptFraction() = " << stepFraction() << endl;
        }

        // External forcing  
        vector Uc = cloud.extForcInt()->interpolate(*this);

        U_ = UBEV+Uc;
        
        const scalar f = 1. - stepFraction();
        trackToAndHitFace(f*trackTime*U_, f, cloud, td);
      
        if (onBoundaryFace() && td.keepParticle)
        {  
           if (isA<wallPolyPatch>(pbMesh[patch()]))
            {        
               const tetIndices tetIs = this->currentTetIndices();
               vector nw = tetIs.faceTri(mesh()).normal();
               nw /= mag(nw);
    
               U_ = -nw*cloud.wallRepX()/((1. - stepFraction())*trackTime);
             
               onWall = true;  
            }           
        }
        
        // The particle hits the wall, moves back in the normal direction by wallRepX()
        // and exits the loop, even if some time still remains. Does not care about
        // colision with internal faces.  
        if (onWall)
        {
            const scalar f = 1. - stepFraction();
            track(f*trackTime*U_, f);
            break;
        }
        
        cnt++;       
    }
    
    // Rescue the particle if needed (keeps previous U_)
    if (cnt>=NMAXITER_)
     {
        const scalar f = 1. - stepFraction();
        track(f*trackTime*U_, f);
        WarningIn("solidParticle::move()")
        << "\nRescuing the particle at position " << position() << nl
        << endl;
     }
    
    if (!td.keepParticle)
     {
       // Add the corresponding molecule to the list to be deleted
       cloud.molcToDelete().append(molcID_);
     }   
 }
 else
 {    
    
    while (td.keepParticle && !td.switchProcessor && stepFraction() < PREC_)
    {
       if (debug)
        {
            Info<< "Time = " << mesh().time().timeName()
                << " trackTime = " << trackTime
                << " steptFraction() = " << stepFraction() << endl;
        }
      
       const scalar f = 1. - stepFraction();
       trackToAndHitFace(f*trackTime*U_, f, cloud, td);

       if (onBoundaryFace() && td.keepParticle)
        {  
           if (isA<wallPolyPatch>(pbMesh[patch()]))
            {      
               // When doing the spring-induced movement it is better to not 
               // add any repulsion since this can violates the spring contraint.
               // Avoiding repulsion also does not guarantees that the
               // spring constraint is kept. The bead remains on the wall.
               td.keepParticle = true;
               break;
            }           
        }
    }
    
    if (!td.keepParticle)
     {
       // Add the corresponding molecule to the list to be deleted
       cloud.molcToDelete().append(molcID_);
     }
  
 }
  
 return td.keepParticle; 
}


bool Foam::solidParticle::hitPatch(solidParticleCloud&, trackingData&)
{
    return false;
}


void Foam::solidParticle::hitProcessorPatch
(
    solidParticleCloud&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::solidParticle::hitWallPatch(solidParticleCloud& cloud, trackingData&)
{
    // Moved up;
}


void Foam::solidParticle::transformProperties(const transformer& transform)
{
    particle::transformProperties(transform);
    U_ = transform.transform(U_);
}

// ************************************************************************* //
