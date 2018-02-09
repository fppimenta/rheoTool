/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
 
\*---------------------------------------------------------------------------*/

#ifndef zeroIonicFluxFvPatchScalarFields_H
#define zeroIonicFluxFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class zeroIonicFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class zeroIonicFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data
        
      //- Ionic valence
        scalar zib_;
       
public:

    //- Runtime type information
    TypeName("zeroIonicFlux");


    // Constructors

        //- Construct from patch and internal field
        zeroIonicFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroIonicFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given zeroIonicFluxFvPatchScalarField onto
        //  a new patch
        zeroIonicFluxFvPatchScalarField
        (
            const zeroIonicFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        zeroIonicFluxFvPatchScalarField
        (
            const zeroIonicFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new zeroIonicFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        zeroIonicFluxFvPatchScalarField
        (
            const zeroIonicFluxFvPatchScalarField&,
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
                new zeroIonicFluxFvPatchScalarField(*this, iF)
            );
        }


     
        // Member functions
    
        //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
         
        //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );
              
            virtual tmp<Field<Foam::scalar> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;   
  
     
            virtual tmp<Field<Foam::scalar> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;
            
                
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
