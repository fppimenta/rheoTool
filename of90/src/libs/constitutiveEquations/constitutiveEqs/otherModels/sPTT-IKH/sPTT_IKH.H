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
    sPTT_IKH

SourceFiles
    sPTT_IKH.C

\*---------------------------------------------------------------------------*/

#ifndef sPTT_IKH_H
#define sPTT_IKH_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class sPTT_IKH Declaration
\*---------------------------------------------------------------------------*/

class sPTT_IKH
:
    public constitutiveEq 
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;
        
        //- Back strain
        volSymmTensorField A_;
        
        //- Structure parameter     
        volScalarField lambda_;
         
        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity (absent in original model)
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity (when lambda = 1)
            dimensionedScalar eta0_;
            
            //- Viscosity only used for stabilization purposes
            dimensionedScalar etaStab_;
            
            //- Elastic modulus
            dimensionedScalar G_;
            
            //- Yield-stress (when lambda = 1)
            dimensionedScalar k0_;
            
            //- PTT parameter
            dimensionedScalar eps_;
            
            //- Back-stress modulus (when lambda = 1)
            dimensionedScalar C0_;
            
            //- Back stress constitutive parameter
            dimensionedScalar q_;
            
            //- Lambda equation parameters
            dimensionedScalar k1_;          
            dimensionedScalar k2_;
            dimensionedScalar k3_;
            dimensionedScalar n1_;           
            dimensionedScalar n2_;            
            dimensionedScalar n3_;
            
            //- Scaling indexes 
            dimensionedScalar m1_;           
            dimensionedScalar m2_;
            dimensionedScalar m3_;
           
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        sPTT_IKH(const sPTT_IKH&);

        //- Disallow default bitwise assignment
        void operator=(const sPTT_IKH&);
        
protected:

       //- Return the solvent viscosity
       virtual const dimensionedScalar etaS() const
       {
          return etaS_;
       }
      
       //- Return the polymeric viscosity
       virtual const dimensionedScalar etaP() const
       {
          return etaStab_;
       }
 
public:

    //- Runtime type information
    TypeName("sPTT-IKH");

    // Constructors

        //- Construct from components
        sPTT_IKH
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~sPTT_IKH()
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
