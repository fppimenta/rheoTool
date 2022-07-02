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
    ML_IKH

SourceFiles
    ML_IKH.C

\*---------------------------------------------------------------------------*/

#ifndef ML_IKH_H
#define ML_IKH_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class ML_IKH Declaration
\*---------------------------------------------------------------------------*/

class ML_IKH
:
    public constitutiveEq 
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;
        
        //- Back strain
        volSymmTensorField A_;
         
        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;

            //- List of weights for of ML model 
            List<scalar> C_;
            
            //- List of Di parameters for each lambda
            List<scalar> D_; 
            
            //- Exponent of ML model       
            dimensionedScalar a_;
            
            //- Exponent of ML model  
            dimensionedScalar n_;
            
            //- Exponent of ML model 
            dimensionedScalar b_;
            
            //- Yield stress
            dimensionedScalar ky_;
            
            //- Hardening coefficient
            dimensionedScalar kh_;
            
            //- Parameter of ML model
            dimensionedScalar k1_;
            
            //- Parameter of ML model
            dimensionedScalar k2_;
            
            //- Parameter of ML model
            dimensionedScalar k3_;
            
            //- Material constant
            dimensionedScalar q_;
            
            //- Relaxation time in structured state
            dimensionedScalar lambdaE_;
            
            //- Is \phi stress-controlled (rate-controlled otherwise)
            bool isSC_;
            
            //- Vector of solved-for directions in mesh
            vector  dims_;
            
            //- Number of solved-for directions in mesh 
            scalar  nDims;
            
            //- Identity tensor consistent with nDims
            dimensionedSymmTensor Itensor;
            
        // Each lambda of the ML model    
        class lambda
          {   
              
            //- Model parameters
            
              // lambdaI field     
                volScalarField lambdaI_;
              
            // Private Member Functions

              //- Disallow default bitwise copy construct
                lambda(const lambda&);

              //- Disallow default bitwise assignment
                void operator=(const lambda&);
                
              public:
              
               // Constructor               
                lambda
                (
                  const word& name,
                  const volVectorField& U,
                  const surfaceScalarField& phi,
                  const dictionary& dict
                );
              
               // Destructor
                ~lambda() {}
             
             // Public Member Functions 
             
             // Access to lambdaI   
              volScalarField& lambdaI()
               {
                   return lambdaI_;
               }        
          };
          
        //- List of species
        PtrList<lambda> lambdas_;   
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        ML_IKH(const ML_IKH&);

        //- Disallow default bitwise assignment
        void operator=(const ML_IKH&);
        
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
    TypeName("ML-IKH");

    // Constructors

        //- Construct from components
        ML_IKH
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~ML_IKH()
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
