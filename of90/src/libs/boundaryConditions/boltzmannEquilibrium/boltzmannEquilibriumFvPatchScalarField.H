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
    boltzmannEquilibriumFvPatchScalarField

Description
    Ionic concentration is computed assuming a Boltzmann distribution.
    The concentration is computed using the intrinsic potential read
    from the same patch. To be used with the PNP model.
 
    Example:
    
      {
         type            boltzmannEquilibrium;
	 c0              9.43798e-05;           // Reference concentration
	 psi0            0;                     // Reference intrinsic potential 
	 value           uniform 9.43798e-05;   // Place-holder
      }
      
    This BC is part of rheoTool.

SourceFiles
    boltzmannEquilibriumFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef boltzmannEquilibriumFvPatchScalarField_H
#define boltzmannEquilibriumFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boltzmannEquilibriumFvPatch Declaration
\*---------------------------------------------------------------------------*/

class boltzmannEquilibriumFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

      // Private data

        //- Ionic valence
        scalar zib_;
        
        //- Reference concentration
        autoPtr<Function1<scalar>>  c0_;
        
        //- Reference potential
         autoPtr<Function1<scalar>>  psi0_;

public:

    //- Runtime type information
    TypeName("boltzmannEquilibrium");


    // Constructors

        //- Construct from patch and internal field
        boltzmannEquilibriumFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        boltzmannEquilibriumFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given boltzmannEquilibriumFvPatchScalarField
        //  onto a new patch
        boltzmannEquilibriumFvPatchScalarField
        (
            const boltzmannEquilibriumFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        

        //- Construct as copy setting internal field reference
        boltzmannEquilibriumFvPatchScalarField
        (
            const boltzmannEquilibriumFvPatchScalarField&,
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
                new boltzmannEquilibriumFvPatchScalarField(*this, iF)
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
