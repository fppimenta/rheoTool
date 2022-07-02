/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::linearExtrapolationFvPatchField

Description

    This BC linearly extrapolates each individual component of the field<Type>
    from interior cells to patch faces. The gradient of each component
    at interior cells is required for this and the scheme used for its
    computation is defined in fvSchemes through keyword 'linExtrapGrad'.
    The value from the previous iteration/time-step is used at the boundary
    face to compute the gradient, so that some explicitness is introduced
    (kind of under-relaxation), which can be reduced by looping the BC. A
    linear regression method is also available, whereby only internal cell
    values are used.
    
    Example:
      
      {
        type            linearExtrapolation;
        useRegression   false; // optional, default is false
        value           uniform (0 0 0 0 0 0);     //place-holder
      }
       
    This BC is part of rheoTool.
    
SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    linearExtrapolationFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef linearExtrapolationFvPatchField_H
#define linearExtrapolationFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class linearExtrapolationFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearExtrapolationFvPatchField
:
    public fixedValueFvPatchField<Type>
{

protected:

    // Protected data
    
    //- Name of the local mass flux
        bool useReg_;

public:

    //- Runtime type information
    TypeName("linearExtrapolation");


    // Constructors

        //- Construct from patch and internal field
        linearExtrapolationFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        linearExtrapolationFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given linearExtrapolationFvPatchField
        //  onto a new patch
        linearExtrapolationFvPatchField
        (
            const linearExtrapolationFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        linearExtrapolationFvPatchField
        (
            const linearExtrapolationFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        linearExtrapolationFvPatchField
        (
            const linearExtrapolationFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new linearExtrapolationFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "linearExtrapolationFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
