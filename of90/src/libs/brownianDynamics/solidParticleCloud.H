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

Class
    Foam::solidParticleCloud

Description
    A Cloud of solid particles
    
    Based on the default OpenFOAM v5.0 "solidParticleCloud" class (src/lagrangian/solidParticle/solidParticleCloud.H).

SourceFiles
    solidParticleCloudI.H
    solidParticleCloud.C
    solidParticleCloudIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticleCloud_H
#define solidParticleCloud_H

#include "Cloud.H"
#include "solidParticle.H"
#include "IOdictionary.H"
#include "externalForcingInterp.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class solidParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class solidParticleCloud
:
    public Cloud<solidParticle>
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;
        
        //- Dictionary with controls
        IOdictionary molcProperties_;
        
        //- List that contains the index of the molecules to deleted 
        List<label> molcToDelete_;
        
        //- Wall-bead distance repulsion. ONLY effective upon collision, ie, 
        // a bead close to the wall, but that does not hit the wall will 
        // not be forced to respect this distance. 
        scalar wallRepX_;
       
        //- Are the molecules tethered
        bool isTethered_;

        //- Member handling external forcing 
        autoPtr<externalForcingInterp>  extForcInt_;

public:

    //- Runtime type information
    TypeName("solidParticleCloud");

    // Constructors

        //- Construct given mesh
        solidParticleCloud
        (
            const fvMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );
        
        //- Disallow default bitwise copy construction
        solidParticleCloud(const solidParticleCloud&) = delete;


    // Member Functions

        // Access

            //- Return reference to mesh 
            inline const fvMesh& mesh() const;
            
            //- Return list of molecules to delete
            inline List<label>& molcToDelete();
            
            //- Return repulsion distance from the wall
            inline const scalar& wallRepX() const;
            
            //- Return true if tethered
            inline const bool& isTethered() const;
            
            //- Return extForcInt_
            inline const autoPtr<externalForcingInterp>& extForcInt() const;

        // Edit

            //- Move the particles. Argument controls whether to add the 
            // the drag force contribution or not
            List<label>& move(bool);
            
       // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidParticleCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
