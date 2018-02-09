/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    HBprofileFvPatchVectorField

Description
    Boundary condition which imposes a fully-developped velocity profile
    for an HerschelBulkley fluid, flowing in a pipe. Allows to specify the 
    radius of the pipe (R) and the average velocity (uIn).
    The BC is used in the Aneurysm tutorial of rheoTool.
    
    Usage example:
    
       type           HBprofile;
       R              1.;
       uIn            (1. 0 0);
       value          uniform ( 0 0 0 );
    
    This BC is part of rheoTool.

SourceFiles
    HBprofileFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef HBprofileFvPatchVectorField_H
#define HBprofileFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class HBprofileFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class HBprofileFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
         
        scalar R_;    
       
        vector uIn_;
     
public:

    //- Runtime type information
    TypeName("HBprofile");


    // Constructors

        //- Construct from patch and internal field
        HBprofileFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        HBprofileFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given HBprofileFvPatchVectorField
        //  onto a new patch
        HBprofileFvPatchVectorField
        (
            const HBprofileFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new HBprofileFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        HBprofileFvPatchVectorField
        (
            const HBprofileFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new HBprofileFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

         
        scalar& R()
        {
            return R_;
        }

        vector& uIn()
        {
            return uIn_;
        }

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
