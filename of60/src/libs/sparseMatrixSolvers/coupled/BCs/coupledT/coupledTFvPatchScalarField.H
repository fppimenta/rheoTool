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
    coupledTFvPatchScalarField

Description
    Boundary condition for temperature to ensure implicit/explicit continuity
    of DIFFUSIVE fluxes between different mesh regions. To be used in multi-region
    solvers. It implements and uses its own deltas and weights, which are computed
    in the same way as for a cyclic AMI fvPatch. Besides this, both gradientInternalCoeffs() 
    and gradientBoundaryCoeffs() are re-implemented to behave as a cylicAMI for diffusion. 
    ValueBoundaryCoeffs() and valueInternalCoeffs() are re-implemented also but are 
    not used in practice, since u.n = 0 is assumed on the patch, ie, no convective
    fluxes of energy across the patch.
    
    Note 1: when the BC is explicit, the temperature on the patch is computed based on the formulas
    for diffusive fluxes continuity and the gradientCoeffs() are the same as for any
    fixedValueFvPatchScalarField.
    
    Note 2: this field BC evaluates coupled() as false.
    
    Important limitation: u.n = 0 is assumed on the patch (no fluid crosses the interface).
    
    This file is part of rheoTool.

\*---------------------------------------------------------------------------*/

#ifndef coupledTFvPatchScalarField_H
#define coupledTFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "regionCoupledBaseFvPatch.H"
#include "regionCoupledAMIFvPatch.H"
#include "regionCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coupledTFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class coupledTFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

private:

    // Private data
    
    //- Name of the temperature field in the neighbouring region
    word nbrFieldName_;
    
    //- Name of the thermal conductivity in this mesh region
    word kappaName_;
    
    //- Is there contact resistance
    bool isContactResistance_;
    
    //- Value of the contact resistance coefficient
    scalar hres_;
    
    //- Is the BC being used within a coupled solver (in which case the 
    // BC is implicit)?
    bool isCoupledSystem;
      
    // Private functions

    //- Return the patch-normal component of the deltas (<=> cyclicAMI)
    tmp<scalarField> deltan() const;

    //- Return the neighbour patch-normal component of the deltas (<=> cyclicAMI)
    tmp<scalarField> nbrDeltan() const;
    
    //- Return weights (w = Nf/(Nf+Pf)) (<=> cyclicAMI)
    tmp<scalarField> weights() const;
    
    //- Return delta (P to N) vectors across coupled patch (<=> cyclicAMI)
    tmp<vectorField> delta() const;
    
    //- Return deltaCoeffs() (<=> cyclicAMI)
    tmp<scalarField> deltaCoeffs() const;
        
public:

    //- Runtime type information
    TypeName("coupledT");


    // Constructors

        //- Construct from patch and internal field
        coupledTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given coupledTFvPatchScalarField
        //  onto a new patch
        coupledTFvPatchScalarField
        (
            const coupledTFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coupledTFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        coupledTFvPatchScalarField
        (
            const coupledTFvPatchScalarField&,
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
                new coupledTFvPatchScalarField(*this, iF)
            );
        }


    // Member functions
    
        // Access functions
            
            //- Name of coupled field 
            word nbrFieldName() const;
            
            //- Name of coupled region
            word nbrMeshName() const;

        // Evaluation functions
        
            //- The implicit behaviour of the BC must be set by the
            // solver using it (explicit is default behaviour)
            void setFlagCoupled(bool flag)
            {
              isCoupledSystem = flag;
            }
        
            //- Return patch-normal gradient
            virtual tmp<Field<scalar>> snGrad() const;

            //- Return patch-normal gradient for coupled-patches
            //  using the deltaCoeffs provided
            virtual tmp<Field<scalar>> snGrad
            (
                const scalarField& deltaCoeffs
            ) const;
             
            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
            
            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<scalar>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
