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
    Foam::regionCoupledAMIFvPatch

Description
    The fvPatch from regionCoupledAMIPolyPatch.
    
    NOTE: deltas and weights are computed as if there was no communication between
    regions, ie, in the same way as for a wall in a single-region mesh. This is so
    such that any regular (non-coupled) field in a multiregion mesh behaves as if 
    it was in a single region mesh. Thus, any region-coupled field needs to implement
    and use the (AMI) coupled deltas and weights.
    
    Important: non-coupled fvPatch (coupled() = false).
    
    This file is part of rheoTool.

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledAMIFvPatch_H
#define regionCoupledAMIFvPatch_H

#include "fvPatch.H"
#include "fvMesh.H"
#include "Time.H"
#include "regionCoupledAMIPolyPatch.H"
#include "regionCoupledBaseFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class regionCoupledAMIFvPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledAMIFvPatch
:
    public lduInterface,
    public fvPatch,
    public regionCoupledBaseFvPatch
{
    // Private data

        // The polyPatch
        const regionCoupledAMIPolyPatch& regionCoupledAMIPolyPatch_;

    // Private members

        
        
protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

public:

    //- Runtime type information
    TypeName(regionCoupledAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        regionCoupledAMIFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            regionCoupledBaseFvPatch
            (
                patch,
                *this
            ),
            regionCoupledAMIPolyPatch_
            (
                refCast<const regionCoupledAMIPolyPatch>(patch)
            )
        {}


    //- Destructor
    ~regionCoupledAMIFvPatch()
    {}


    // Member Functions


        // Access
            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;
            
            //- Return faceCell addressing
            virtual const labelUList& faceCells() const
            {
                return fvPatch::faceCells();
            }

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return regionCoupledAMIPolyPatch_.coupled();
            }
            
            const regionCoupledAMIPolyPatch& polyPatch() const
            {
                return regionCoupledAMIPolyPatch_;
            }
            
            //- Return regionCoupledAMIFvPatch nbr
            const regionCoupledAMIFvPatch& neighbFvPatch() const
            {
               return refCast<const regionCoupledAMIFvPatch>
               (
                 nbrFvMesh().boundary()
                 [
                    neighbPatchID()
                 ]
               );
            }


        // Interface transfer functions

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const labelUList& internalData
            ) const;

            //- Inherit initInternalFieldTransfer from lduInterface
            using lduInterface::initInternalFieldTransfer;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                labelUList& iF
            ) const
            {}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const labelUList& iF
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
