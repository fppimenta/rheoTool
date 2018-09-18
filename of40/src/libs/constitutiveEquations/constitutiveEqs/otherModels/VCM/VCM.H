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
    VCM

SourceFiles
    VCM.C

\*---------------------------------------------------------------------------*/

#ifndef VCM_H
#define VCM_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class VCM Declaration
\*---------------------------------------------------------------------------*/

class VCM
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;
        
        //- Dimensionless number density of A (long chains)
        volScalarField nA_;
        
        //- Dimensionless number density of B (short chains)
        volScalarField nB_;
        
        //- Configuration tensor A
        volSymmTensorField A_;
        
        //- Configuration tensor B
        volSymmTensorField B_;

        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity (only for stabilization purposes)
            dimensionedScalar etaP_;

            //- Relaxation time of A
            dimensionedScalar lambdaA_;
            
            //- Diffusion coefficient of A
            dimensionedScalar DA_;
            
            //- Diffusion coefficient of B
            dimensionedScalar DB_;
            
            //- Parameter of the breakage event equation (\xi*\mu)
            dimensionedScalar chi_;
            
            //- Dimensionless equilibrium breakage rate of A
            dimensionedScalar cAEq_;
            
            //- Dimensionless equilibrium breakage rate of B
            dimensionedScalar cBEq_;
            
            //- Plateau elastic modulus 
            dimensionedScalar G0_;
            
            //- Ratio of relaxation times (lambdaB/lambdaA) 
            dimensionedScalar eps_;
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        VCM(const VCM&);

        //- Disallow default bitwise assignment
        void operator=(const VCM&);
        
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
    TypeName("VCM");

    // Constructors

        //- Construct from components
        VCM
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~VCM()
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
