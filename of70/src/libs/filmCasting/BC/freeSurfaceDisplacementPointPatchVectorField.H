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

Class
    Foam::freeSurfaceDisplacementPointPatchVectorField

Description
    Foam::freeSurfaceDisplacementPointPatchVectorField

SourceFiles
    freeSurfaceDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef freeSurfaceDisplacementPointPatchVectorField_H
#define freeSurfaceDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "SubField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class freeSurfaceDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class freeSurfaceDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data
    
       // Anchor points at inlet (they are fixed)
       List<label> anchorPoints_;
    
       // Points shared with outlet patch
       List<label> outletPoints_;
       
       // Stencil for each point (previous and next points)
       List<List<label>> stencil_;
       
       // Which method to use: Peric or streamline
       word method_;
       
       // Limiter of the advective term in the streamline equation
       word limiterFunc_;
       
       // Under-relaxation factor of the free-surface position
       scalar URF_;
       
       // Replace the normal component of velocity by the flux
       bool useFlux_;
        
       // Coord of the freeSurface points in constant/polymesh, i.e. always for t = 0
       // Used as reference by the displacament solver
       vectorField pts0_;
       
       // Point displacement: dX_ = pts_now - pts0_
       vectorField dX_;
       
       // y-coord of the freeSurface points at the previous time-step (not previous inner iter)
       scalarField y0_;
       
       // Ghost auxilliary points in Peric's method (one point per face)
       vectorField auxF_;
       
       // Motion direction in Peric's method
       vector motionDir_;
  
       
    // Private functions
       
       // Make the stencil for the model mesh provided (it is not a general algo), either sequential or parallel, but if parallel
       // decomposePar should conserve the freeSurface patch in a single processor 
       void buildStencil();
       
       // Interpolate velocity from faces to points (combining phi and U)
       void interpolateUonPoints(vectorField&, const scalarField&);
    
       // Integrate the streamline equation
       void integrate(scalarField&, const vectorField&, const vectorField&, const scalar&, scalar&, scalar&);
    
public:

    //- Runtime type information
    TypeName("freeSurfaceDisplacement");


    // Constructors

        //- Construct from patch and internal field
        freeSurfaceDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeSurfaceDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        freeSurfaceDisplacementPointPatchVectorField
        (
            const freeSurfaceDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector>> clone() const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new freeSurfaceDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        freeSurfaceDisplacementPointPatchVectorField
        (
            const freeSurfaceDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new freeSurfaceDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
            
            //- Update the freeSurface points (must be called before updateCoeffs())
            void updatePoints(const bool&,const scalarField&);
            
            //- Update the freeSurface points (must be called before updateCoeffs())
            void updatePoints(const scalarField&);
            
            //- Retrieve method
            const word& method
            (
            ) const
            {
              return method_;
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
