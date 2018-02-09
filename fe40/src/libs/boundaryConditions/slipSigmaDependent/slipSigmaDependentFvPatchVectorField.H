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
    slipSigmaDependentFvPatchVectorField

Description
    Slip velocity following the Helmoltz-Smoluchowski theory, for
    conductivity-dependent zeta-potentials.
    It considers: E=-grad(phiE). 
    
    Example:
    
      {
        type            slipSigmaDependent;
	elecMobility0   3e-08;                // reference electroosmotic mobility
	sigma0          4;                    // reference conductivity
	m               -0.3;                 // power-law exponent
	value           uniform (0 0 0);      // place-holder 
      }
      
    This BC is part of rheoTool.
    
SourceFiles
    slipSigmaDependentFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef slipSigmaDependentFvPatchVectorField_H
#define slipSigmaDependentFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class slipSigmaDependentFvPatch Declaration
\*---------------------------------------------------------------------------*/

class slipSigmaDependentFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private data
      
       //- Reference electric mobility (=U/E)      
        scalar elecM_;
        
       //- Reference conductivity      
        scalar sigma0_;
        
       //- Coefficient for the power scaling of zeta      
        scalar m_; 

public:

    //- Runtime type information
    TypeName("slipSigmaDependent");


    // Constructors

        //- Construct from patch and internal field
        slipSigmaDependentFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        slipSigmaDependentFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given slipSigmaDependentFvPatchVectorField
        //  onto a new patch
        slipSigmaDependentFvPatchVectorField
        (
            const slipSigmaDependentFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new slipSigmaDependentFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        slipSigmaDependentFvPatchVectorField
        (
            const slipSigmaDependentFvPatchVectorField&,
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
                new slipSigmaDependentFvPatchVectorField(*this, iF)
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
