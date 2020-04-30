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
    frozenTemperatureMixture

Description
    Thermo model with frozen temperature field, ie, T field must be provided but is
    not evolved in time/space. 
    This file is part of rheoTool.
    
\*---------------------------------------------------------------------------*/

#ifndef frozenTemperatureMixture_H
#define frozenTemperatureMixture_H

#include "radiationModel.H"
#include "noRadiation.H"
#include "fluidThermoMixtureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fluidThermoMixtureModels
{ 

/*---------------------------------------------------------------------------*\
                        Class frozenTemperatureMixture Declaration
\*---------------------------------------------------------------------------*/

class frozenTemperatureMixture
:
    public fluidThermoMixtureModel
{
private:

    // Private data
          
    // Private Member Functions

        //- Disallow default bitwise copy construct
        frozenTemperatureMixture(const frozenTemperatureMixture&);

        //- Disallow default bitwise assignment
        void operator=(const frozenTemperatureMixture&);
        
        
protected:

   // Protected Members
   
public:

    //- Runtime type information
    TypeName("frozenTemperature");

    // Constructors

        //- Construct from components
        frozenTemperatureMixture
        (
            const word& name,
            const fvMesh& mesh,
            const word& phase1,
            const word& phase2
        );

    // Destructor
        virtual ~frozenTemperatureMixture()
        {}

    // Member Functions
            
        //- Update thermo        
        virtual void correct
        (          
          const volVectorField& U,
          const surfaceScalarField& rhoPhi,
          const volSymmTensorField& tau,
          const volScalarField& alpha1,
          const volScalarField& rho,
          fv::options&,
          int nNonOrtoC
        );
      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidThermoMixtureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
