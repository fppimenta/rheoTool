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
    PowerLaw

SourceFiles
    PowerLaw.C

\*---------------------------------------------------------------------------*/

#ifndef PowerLaw_H
#define PowerLaw_H

#include "constitutiveEq.H"
#include "thermoFunction.H" 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class PowerLaw Declaration
\*---------------------------------------------------------------------------*/

class PowerLaw
:
    public constitutiveEq
{
    // Private data

         // Model constants

            //- Density
            dimensionedScalar rho_;
  
            //- Lower bound for viscosity 
            dimensionedScalar etaMin_;

            //- Upper bound for viscosity 
            dimensionedScalar etaMax_;

            //- Consistency index
            dimensionedScalar k_;

            //- Behavior index
            dimensionedScalar n_;

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        //- Viscosity field
        volScalarField eta_;

        //- Thermofunction for temperature dependence of viscosity
        autoPtr<thermoFunction> thermoEtaPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        PowerLaw(const PowerLaw&);

        //- Disallow default bitwise assignment
        void operator=(const PowerLaw&);

protected:
   
       //- Return the viscosity field
       virtual const volScalarField& eta() const
       {
          return eta_;
       }
       
       // Is the model prepared to work in non-isothermal conditions
       virtual bool hasThermo() const
       {
         return true;
       }
      
public:

    //- Runtime type information
    TypeName("PowerLaw");

    // Constructors

        //- Construct from components
        PowerLaw
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~PowerLaw()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const
        {
            return tau_;
        }

        //- Return the density
        virtual const dimensionedScalar rho() const
        {
            return rho_;
        }

        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const
        {
          return true;
        };

        //- Correct the viscoelastic stress: alpha is the color function in two phase-flows
        virtual void correct
        (
          const volScalarField* alpha = nullptr,
          const volTensorField* gradU = nullptr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
