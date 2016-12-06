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
    extSTFvPatchSymmTensorField

Description

    This BC linearly extrapolates each individual tensor component of a
    symmTensorField from interior cells to patch faces. The gradient of
    each tensor component at interior cells is required for this and the
    scheme used for its computation is defined in fvSchemes through
    keyword 'extSTGrad'. The value from the previous iteration/time-step is
    used at the boundary face to compute the gradient, so that some 
    explicitness is introduced, which can be reduced by internally looping
    the BC.
    
    Usage:
    
        type            extST;
        value           uniform (0 0 0 0 0 0); //place-holder
        
    This BC is part of rheoTool.

SourceFiles
    extSTFvPatchSymmTensorField.C

\*---------------------------------------------------------------------------*/

#ifndef extSTFvPatchSymmTensorField_H
#define extSTFvPatchSymmTensorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class extTetaFvPatch Declaration
\*---------------------------------------------------------------------------*/

class extSTFvPatchSymmTensorField
:
    public fixedValueFvPatchSymmTensorField
{

public:

    //- Runtime type information
    TypeName("extST");


    // Constructors

        //- Construct from patch and internal field
        extSTFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        extSTFvPatchSymmTensorField
        (
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given extSTFvPatchSymmTensorField
        //  onto a new patch
        extSTFvPatchSymmTensorField
        (
            const extSTFvPatchSymmTensorField&,
            const fvPatch&,
            const DimensionedField<symmTensor, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchSymmTensorField> clone() const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new extSTFvPatchSymmTensorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        extSTFvPatchSymmTensorField
        (
            const extSTFvPatchSymmTensorField&,
            const DimensionedField<symmTensor, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchSymmTensorField> clone
        (
            const DimensionedField<symmTensor, volMesh>& iF
        ) const
        {
            return tmp<fvPatchSymmTensorField>
            (
                new extSTFvPatchSymmTensorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
