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
    navierSlipFvPatchVectorField

Description

    Implements fully-explicit slip boundary conditions for U. The fluid
    rheology can be generic and the mesh can be either static or moving.
    There is also support for two-phase flows (untested).
    
    Example:
    
      {
        type            navierSlip;
        value           uniform (0 0 0); 
        
        model           nonLinearNavierSlip; // model name
        m               1;                   // parameter of selected model 
        knl             5e-7;                // parameter of selected model 
        URF             1.;                  // under-relaxation factor
        isTwoPhaseFlow  false;               // is it a two phase flow 
        isMovingWall    true;                // is that patch/wall moving (dynamic mesh) 
      }
      
    This BC is part of rheoTool.

SourceFiles
    navierSlipFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef navierSlipFvPatchVectorField_H
#define navierSlipFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class navierSlipFvPatch Declaration
\*---------------------------------------------------------------------------*/

class navierSlipFvPatchVectorField
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
    TypeName("navierSlip");


    // Constructors

        //- Construct from patch and internal field
        navierSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        navierSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given navierSlipFvPatchVectorField
        //  onto a new patch
        navierSlipFvPatchVectorField
        (
            const navierSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new navierSlipFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        navierSlipFvPatchVectorField
        (
            const navierSlipFvPatchVectorField&,
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
                new navierSlipFvPatchVectorField(*this, iF)
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
