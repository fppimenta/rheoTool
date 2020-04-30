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
    boussinesq

Description
    Thermo model enabling the boussinesq approximation for density, which remains
    optional. Thermal conductivity is a temperature-dependent coefficient, but
    isotropic. Radiation is optional.
    This file is part of rheoTool.
    
\*---------------------------------------------------------------------------*/

#ifndef boussinesq_H
#define boussinesq_H

#include "radiationModel.H"
#include "noRadiation.H"
#include "UniformDimensionedField.H"
#include "uniformDimensionedFields.H"
#include "fluidThermoModel.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermoModels
{ 

/*---------------------------------------------------------------------------*\
                        Class boussinesq Declaration
\*---------------------------------------------------------------------------*/

class boussinesq
:
    public fluidThermoModel
{
private:

    // Private data
 
        //- Thermal conductivity
        volScalarField kappa_;
        
        //- Thermal expansion coeff
        dimensionedScalar beta_;
        
        //- Reference temperature for thermal expansion
        dimensionedScalar TRef_;
        
        //- Enable natural convection?
        bool isNatConvection_;
        
        //- Enable viscous dissipation?
        bool isViscDissipation_;
        
        //- Radiation model
        autoPtr<radiationModel> radiation_;
        
        //- Density x specific heat
        dimensionedScalar rhoCp_;
        
        //- Gravity field
        uniformDimensionedVectorField g_;
        
        //- Reference heigth
        uniformDimensionedScalarField hRef_;
        
        //- Coefficients for the thermal conductivity
        dimensionedScalar k0_;
        dimensionedScalar k1_;
        dimensionedScalar k2_;
        dimensionedScalar k3_;
         
    // Private Member Functions

        //- Disallow default bitwise copy construct
        boussinesq(const boussinesq&);

        //- Disallow default bitwise assignment
        void operator=(const boussinesq&);
    
        
protected:

   // Protected Members
   
public:

    //- Runtime type information
    TypeName("boussinesq");

    // Constructors

        //- Construct from components
        boussinesq
        (
            const word& name,
            const fvMesh& mesh
        );

    // Destructor
        virtual ~boussinesq()
        {}

    // Member Functions
         
        //- Return boussinesq density
        virtual tmp<volScalarField> rhok() const
        {
          if (isNatConvection_) 
            return 1. - beta_*(T() - TRef_);
          else
            return beta_*T()*0.;
        }
        
        //- Is gravity used in the model? Evaluation is based upon the gravity field: if
        // g == (0,0,0), then there is no gravity and isNatConvection_ = false. 
        virtual bool hasGravity() const
        {
          if (isNatConvection_) 
            return true;
          else
            return false;
        }    
        
        //- Return g.h for models that have hasGravity() = true
        virtual tmp<volScalarField> gh() const
        {
           dimensionedScalar ghRef(- mag(g_)*hRef_);
           return (g_ & mesh().C()) - ghRef;
        }
        
        //- Return g.hf for models that have hasGravity() = true
        virtual tmp<surfaceScalarField> ghf() const
        {
           dimensionedScalar ghRef(- mag(g_)*hRef_);
           return (g_ & mesh().Cf()) - ghRef;
        }       
                     
        //- Update thermo fields
        virtual void correct
        (          
          const volVectorField& U,
          const surfaceScalarField& phi,
          const volSymmTensorField& tau,
          autoPtr<coupledSolver>& cpsT,
          fv::options&,
          int nNonOrtoC
        );
        
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidThermoModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
