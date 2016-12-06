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
    uCosFvPatchVectorField

Description
    Boundary condition which speficies a smooth startup for a vector field.
    The BC is used in the 4:1 contraction tutorial of rheoTool.
    
    Usage example:
    
       type           uCos;
       tlim           1.;
       fac            8.;
       uav            (0.25 0 0);
       dirN           (1 0 0);
       value          uniform ( 0 0 0 );
    
    This BC is part of rheoTool.

SourceFiles
    uCosFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef uCosFvPatchVectorField_H
#define uCosFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class uCosFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class uCosFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
         
        scalar tlim_;
        
        scalar fac_;
       
        vector uav_;
        
        vector dirN_;

public:

    //- Runtime type information
    TypeName("uCos");


    // Constructors

        //- Construct from patch and internal field
        uCosFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uCosFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uCosFvPatchVectorField
        //  onto a new patch
        uCosFvPatchVectorField
        (
            const uCosFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new uCosFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uCosFvPatchVectorField
        (
            const uCosFvPatchVectorField&,
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
                new uCosFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

         
        scalar& tlim()
        {
            return tlim_;
        }

        
        scalar& fac()
        {
            return fac_;
        }

         
        vector& uav()
        {
            return uav_;
        }
        
        vector& dirN()
        {
            return dirN_;
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
