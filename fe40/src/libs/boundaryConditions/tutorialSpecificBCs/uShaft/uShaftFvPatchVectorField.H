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
    uShaftFvPatchVectorField

Description

    Derived from navierSlip BC.
      
    This BC is part of rheoTool.

SourceFiles
    uShaftFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef uShaftFvPatchVectorField_H
#define uShaftFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uShaftFvPatch Declaration
\*---------------------------------------------------------------------------*/

class uShaftFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private data
    
       //- Model name
       word model_;
      
       //- Exponent of tau (non-linar Navier slip)         
       scalar m_;
         
       //- Slip factor (non-linar Navier slip) 
       scalar knl_;
       
       //- Exponent of tau (slip TT model)         
       scalar alpha_;
         
       //- Slip factor (slip TT model) 
       scalar beta_;
       
       //- Solvent viscosity for VE models
       scalar etaS_;
       
       //- Under-relaxation factor for uws (only static meshes)
       scalar URF_;
       
       //- Is there two phases (rheoInterFoam)
       bool isTwoPhaseFlow_;
       
       //- Is the wall moving (only for dynamic meshes)
       bool isMovingWall_;
       
       //- Sinusoidal motion parameters
       vector offset_;
       vector amplitude_;
       scalar frequency_;
       scalar delay_;
       
       //- This a pointer to the total extra-stress tensor (solvent + polymer).
       //  Shared by all objects, in order to avoid N calls to tauTotal(),
       //  retrieving the same values. 1st CPU, 2nd RAM. Not a very elegant
       //  solution; may need to be revisited.
       static autoPtr<volSymmTensorField> tauTotalPtr;
       
       //- Updates the total extra-stress tensor when needed 
       static void updateTauTotalPtr(const volSymmTensorField&);
       
       //- Caller to static updateTauTotalPtr() 
       void updateTau();
       
       //- Computes uws according to the selected model
       void calcUws(vectorField&, vectorField&);
    

public:

    //- Runtime type information
    TypeName("uShaft");


    // Constructors

        //- Construct from patch and internal field
        uShaftFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uShaftFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given uShaftFvPatchVectorField
        //  onto a new patch
        uShaftFvPatchVectorField
        (
            const uShaftFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new uShaftFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uShaftFvPatchVectorField
        (
            const uShaftFvPatchVectorField&,
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
                new uShaftFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
            
            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
