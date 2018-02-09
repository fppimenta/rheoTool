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
    boltzmannEquilibriumLogFvPatchScalarField

Description

SourceFiles
    boltzmannEquilibriumLogFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef boltzmannEquilibriumLogFvPatchScalarField_H
#define boltzmannEquilibriumLogFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boltzmannEquilibriumLogFvPatch Declaration
\*---------------------------------------------------------------------------*/

class boltzmannEquilibriumLogFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

      // Private data

        //- Ionic valence
        scalar zib_;
        
        //- Reference concentration
        scalar  c0_;
        
        //- Reference potential
        scalar psi0_;

public:

    //- Runtime type information
    TypeName("boltzmannEquilibriumLog");


    // Constructors

        //- Construct from patch and internal field
        boltzmannEquilibriumLogFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        boltzmannEquilibriumLogFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given boltzmannEquilibriumLogFvPatchScalarField
        //  onto a new patch
        boltzmannEquilibriumLogFvPatchScalarField
        (
            const boltzmannEquilibriumLogFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new boltzmannEquilibriumLogFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        boltzmannEquilibriumLogFvPatchScalarField
        (
            const boltzmannEquilibriumLogFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new boltzmannEquilibriumLogFvPatchScalarField(*this, iF)
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

#endif

// ************************************************************************* //
