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
    Foam::solidParticle

Description
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.
    
    Based on the default OpenFOAM "solidParticle" class (src/lagrangian/solidParticle).

SourceFiles
    solidParticleI.H
    solidParticle.C
    solidParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidParticle_H
#define solidParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidParticleCloud;
class sPCloudInterface;

// Forward declaration of friend functions and operators

class solidParticle;

Ostream& operator<<(Ostream&, const solidParticle&);

/*---------------------------------------------------------------------------*\
                           Class solidParticle Declaration
\*---------------------------------------------------------------------------*/

class solidParticle
:
    public particle
{
    // Private member data
 
        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;
        
        //- Velocity of parcel
        vector U_;
        
        //- Indexes
        Field<label> ids_;
        
        //- Molecule ID
        label  molcID_;
        
        //- Tracking precision (1 - full precision; 0 - null precision)
        //  Tracking in OpenFOAM applications uses 1, but sometimes fractions
        //  as low as 1e-15 remain to be completed and there is no need to 
        //  do an extra cycle.         
        static const scalar PREC_;
        
        //- Sometimes particles can enter a ping-pong game between the cells
        //  sharing a face, leading to an infinite while. This imposes a 
        //  maximum number of face hits (loops). It should be a high number
        //  in order to allow the cross of multiple internal faces (high dt),
        //  for the situations where it is really needed. It was shown that
        //  rescuing is mostly needed with Gradient interpolation.
        static const int NMAXITER_;


public:

    friend class Cloud<solidParticle>;
    friend class sPCloudInterface;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
       public particle::trackingData
    {
        
            const bool includeDrag_;
            
    public:

        // Constructors

            inline trackingData
            (
                const solidParticleCloud& spc,
                const bool& includeDrag_
            );

        // Member functions
            
            inline const bool& includeDrag() const;          
    };


    // Constructors

        //- Construct from components
        inline solidParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const scalar d,
            const vector& U
        );

        //- Construct from Istream
        solidParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new solidParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<solidParticle> operator()(Istream& is) const
            {
                return autoPtr<solidParticle>
                (
                    new solidParticle(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return bead velocity
            inline const vector& U() const;
            
            //- Return indices (global ID| localID | group ID)
            inline const Field<label>& ids() const;
            
            //- Return molecule ID
            inline const label& molcID() const;


        // Tracking

            //- Move
            bool move(solidParticleCloud&, trackingData&, const scalar);


    // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            bool hitPatch(solidParticleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            void hitProcessorPatch(solidParticleCloud& cloud, trackingData& td);

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(solidParticleCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            virtual void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            virtual void transformProperties(const vector& separation);


    // I-O

        static void readFields(Cloud<solidParticle>& c);

        static void writeFields(const Cloud<solidParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const solidParticle&);
};


template<>
inline bool contiguous<solidParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
