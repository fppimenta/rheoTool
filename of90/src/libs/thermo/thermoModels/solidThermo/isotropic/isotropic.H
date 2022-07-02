/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    isotropic

Description
    Thermo model considering an isotropic, temperature-dependent thermal conductivity.
    Radiation is optional. 
    This file is part of rheoTool.
    
\*---------------------------------------------------------------------------*/

#ifndef isotropic_H
#define isotropic_H

#include "radiationModel.H"
#include "noRadiation.H"
#include "solidThermoModel.H"
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidThermoModels
{ 

/*---------------------------------------------------------------------------*\
                        Class isotropic Declaration
\*---------------------------------------------------------------------------*/

class isotropic
:
    public solidThermoModel
{
private:

    // Private data
 
        //- Thermal conductivity
        volScalarField kappa_;
        
        //- Density x specific heat
        dimensionedScalar rhoCp_;
        
        //- Coefficients for the thermal conductivity
        dimensionedScalar k0_;
        dimensionedScalar k1_;
        dimensionedScalar k2_;
        dimensionedScalar k3_;
        
        //- Radiation model
        autoPtr<radiationModel> radiation_;
         
    // Private Member Functions

        //- Disallow default bitwise copy construct
        isotropic(const isotropic&);

        //- Disallow default bitwise assignment
        void operator=(const isotropic&);
        
protected:

   // Protected Members
   
public:

    //- Runtime type information
    TypeName("isotropic");

    // Constructors

        //- Construct from components
        isotropic
        (
            const word& name,
            const fvMesh& mesh
        );

    // Destructor
        virtual ~isotropic()
        {}

    // Member Functions
                 
        //- Return diffusion number (solids only)
        virtual void DiffNo
        (
           const fvMesh& mesh,
           const Time& runTime,
           scalar& maxDiffNo,
           scalar& meanDiffNo
        ) const;
         
        
        //- Update thermo fields
        virtual void correct
        ( 
          autoPtr<coupledSolver>& cpsT,          
          const fvModels&,
          fvConstraints&,
          int nNonOrtoC
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidThermoModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
