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
    frozenTemperature

Description
    Thermo model with frozen temperature field, ie, T field must be provided but is
    not evolved in time/space. 
    This file is part of rheoTool.
    
\*---------------------------------------------------------------------------*/

#ifndef frozenTemperature_H
#define frozenTemperature_H

#include "radiationModel.H"
#include "noRadiation.H"
#include "fluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermoModels
{ 

/*---------------------------------------------------------------------------*\
                        Class frozenTemperature Declaration
\*---------------------------------------------------------------------------*/

class frozenTemperature
:
    public fluidThermoModel
{
private:

    // Private data
          
    // Private Member Functions

        //- Disallow default bitwise copy construct
        frozenTemperature(const frozenTemperature&);

        //- Disallow default bitwise assignment
        void operator=(const frozenTemperature&);
        
        
protected:

   // Protected Members
   
public:

    //- Runtime type information
    TypeName("frozenTemperature");

    // Constructors

        //- Construct from components
        frozenTemperature
        (
            const word& name,
            const fvMesh& mesh
        );

    // Destructor
        virtual ~frozenTemperature()
        {}

    // Member Functions
            
        //- Is gravity used in the model?        
        virtual bool hasGravity() const
        {
          return false;
        } 
        
        //- Update thermo fields
        virtual void correct
        (          
          const volVectorField& U,
          const surfaceScalarField& phi,
          const volSymmTensorField& tau,
          autoPtr<coupledSolver>& cpsT,
          fvModels&,
          fvConstraints&,
          int nNonOrtoC
        );
      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidThermoModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
