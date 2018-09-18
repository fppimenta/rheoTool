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
    RRM
    TODO: this model should not be inside folder EVP
SourceFiles
    RRM.C

\*---------------------------------------------------------------------------*/

#ifndef RRM_H
#define RRM_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class RRM Declaration
\*---------------------------------------------------------------------------*/

class RRM
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;
        
        //- Orientation tensor
        volSymmTensorField S_;
        
        //- Rod length normalized with the initial length
        volScalarField Lstar_;

        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;
            
            //- Poymer viscosity (ONLY used for stabilization!)
            dimensionedScalar etaP_;
 
            //- Rs constant of proportionality
            dimensionedScalar lambdaS_;
            
            //- Product n0*Kb*t (elastic modulus)
            dimensionedScalar G0_;
            
            //- Rotational diffusion coefficient
            dimensionedScalar Dr0_;
            
            //- Ln initial aspect ratio
            dimensionedScalar m_;
            
            //- Lmax parameter
            dimensionedScalar alpha_;
            
            //- Lmax parameter
            dimensionedScalar beta_;
            
            //- Constant of alignment-induced growth
            dimensionedScalar k_;
                       
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        RRM(const RRM&);

        //- Disallow default bitwise assignment
        void operator=(const RRM&);
        
        //- Approximation of K:<uuuu>
        inline volSymmTensorField KDDFourthMoment
        ( 
          const volSymmTensorField& S,
          const volSymmTensorField& D
        )
        {
          return (1./5.) * symm( (S & D) + (D & S) - (S & S & D) - (D & S & S) + 2.*(S & D & S) + 3.*S*(S && D)  );
        }
        
protected:

       //- Return the solvent viscosity
       virtual const dimensionedScalar etaS() const
       {
          return etaS_;
       }
      
       //- Return the polymeric viscosity
       virtual const dimensionedScalar etaP() const
       {
          return etaP_;
       }
        
 
public:

    //- Runtime type information
    TypeName("RRM");

    // Constructors

        //- Construct from components
        RRM
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~RRM()
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
          return false;
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
