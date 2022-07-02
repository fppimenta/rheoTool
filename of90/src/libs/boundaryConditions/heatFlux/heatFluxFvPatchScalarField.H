/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::heatFluxFvPatchScalarField

Description
    
    This file is part of rheoTool.

\*---------------------------------------------------------------------------*/

#ifndef heatFluxFvPatchScalarField_H
#define heatFluxFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class heatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class heatFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Current time index (used for updating)
        label curTimeIndex_;
        
        //- Heat power [W]
        scalar Q_;

        //- Heat flux [W/m2]
        scalarField q_;

        //- Heat transfer coefficient [W/m2K]
        scalarField h_;

        //- Ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;
        
        //- Surface temperature for radiation [K]
        autoPtr<Function1<scalar>> Ts_;

        //- Relaxation for the wall temperature (thermal inertia)
        scalar relaxation_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivity_;
        
        //- Thermal conductivity field name
        word kappaName_;
        
        //- Aux coeffs
        scalarField alpha;
        
        //- Aux coeffs
        scalarField beta;
 

public:

    //- Runtime type information
    TypeName("heatFlux");


    // Constructors

        //- Construct from patch and internal field
        heatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        heatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given heatFluxFvPatchScalarField onto
        //  a new patch
        heatFluxFvPatchScalarField
        (
            const heatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        heatFluxFvPatchScalarField
        (
            const heatFluxFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        heatFluxFvPatchScalarField
        (
            const heatFluxFvPatchScalarField&,
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
                new heatFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions
 
        //- Update the patch pressure gradient field
        virtual void updateCoeffs();
        
        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> valueInternalCoeffs
        (
          const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual tmp<scalarField> valueBoundaryCoeffs
        (
          const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<scalarField> gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the
        //  evaluation of the gradient of this patchField
        virtual tmp<scalarField> gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
 
#endif

// ************************************************************************* //
