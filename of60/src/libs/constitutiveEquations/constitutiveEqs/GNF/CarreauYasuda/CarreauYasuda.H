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
    CarreauYasuda

SourceFiles
    CarreauYasuda.C

\*---------------------------------------------------------------------------*/

#ifndef CarreauYasuda_H
#define CarreauYasuda_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class CarreauYasuda Declaration
\*---------------------------------------------------------------------------*/

class CarreauYasuda
:
    public constitutiveEq
{
    // Private data

        // Model constants

            //- Density
            dimensionedScalar rho_;
  
            //- Zero strain rate viscosity 
            dimensionedScalar eta0_;

            //- Infinite strain rate viscosity 
            dimensionedScalar etaInf_;

            //- Characteristic time-scale
            dimensionedScalar k_;

            //- Dimensionaless parameters
            dimensionedScalar a_;
            dimensionedScalar n_;

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        //- Viscosity field
        volScalarField eta_;

        


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CarreauYasuda(const CarreauYasuda&);

        //- Disallow default bitwise assignment
        void operator=(const CarreauYasuda&);
        
protected:
   
     //- Return the viscosity field
     virtual const volScalarField& eta() const
       {
          return eta_;
       }

public:

    //- Runtime type information
    TypeName("CarreauYasuda");

    // Constructors

        //- Construct from components
        CarreauYasuda
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~CarreauYasuda()
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

        //- Correct the viscoelastic stress
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
