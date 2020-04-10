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

Class
    BarycentricWeights
    
Description
    The external forcing is interpolated by decomposing the cells in tets and
    using the barycentric weights on the vertices (default method of OF for
    lagrangian particles).

SourceFiles
    BarycentricWeights.C

\*---------------------------------------------------------------------------*/

#ifndef BarycentricWeights_H
#define BarycentricWeights_H

#include "externalForcingInterp.H"
#include "interpolationCellPointBary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace externalForcingInterps
{

/*---------------------------------------------------------------------------*\
                           Class BarycentricWeights Declaration
\*---------------------------------------------------------------------------*/

class BarycentricWeights
:
    public externalForcingInterp
{
    // Private data
        
        //- Is electrophoresis active
        bool epActive_;
        
        //- Electrophoretic mobility
        scalar epMobility_;

        //- Is the external numerical forcing steady in time (analytical forcing is always steady)
        bool frozenFlow_;
        
        //- Continuum field
        volVectorField Uex_;
        
        //- Interpolator
        interpolationCellPointBary<vector> interpolatorCPB_; 
     

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BarycentricWeights(const BarycentricWeights&);

        //- Disallow default bitwise assignment
        void operator=(const BarycentricWeights&);
        
        //- Update Uex 
        void updateUex();


public:

    //- Runtime type information
    TypeName("BarycentricWeights");


    // Constructors

        //- Construct from components
        BarycentricWeights
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        virtual ~BarycentricWeights()
        {}


    // Member Functions

        //- Update
        virtual void update();
        
        //- Interpolate
        virtual inline const vector interpolate(const solidParticle&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace externalForcingInterps
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BarycentricWeightsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
