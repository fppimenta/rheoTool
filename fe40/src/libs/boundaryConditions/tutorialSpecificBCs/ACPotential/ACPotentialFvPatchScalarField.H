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
    ACPotentialFvPatchScalarField

Description
    Intrinsic potential (psi) distribution AC in a conducting surface
    of arbitrary shape. To be used exclusively with 2 potentials: phiE
    and psi, where phiE induces psi's distribution.
 
    Example:
    
      {
         type            ACPotential;
	 psiF            0;                 // Bias voltage on the surface
	 value           uniform 0;         // Place-holder
      }
      
    This BC is part of rheoTool.

SourceFiles
    ACPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ACPotentialFvPatchScalarField_H
#define ACPotentialFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ACPotentialFvPatch Declaration
\*---------------------------------------------------------------------------*/

class ACPotentialFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

      // Private data

        //- Offset
        scalar offset_;
        
        //- Amplitude
        scalar amplitude_;
        
        //- Frequency
        scalar frequency_;
        
        //- Phase delay
        scalar phaseDelay_;
        
public:

    //- Runtime type information
    TypeName("ACPotential");


    // Constructors

        //- Construct from patch and internal field
        ACPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ACPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given ACPotentialFvPatchScalarField
        //  onto a new patch
        ACPotentialFvPatchScalarField
        (
            const ACPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ACPotentialFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        ACPotentialFvPatchScalarField
        (
            const ACPotentialFvPatchScalarField&,
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
                new ACPotentialFvPatchScalarField(*this, iF)
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
