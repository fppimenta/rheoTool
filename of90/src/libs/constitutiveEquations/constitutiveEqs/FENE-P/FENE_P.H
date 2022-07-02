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
    FENE_P

SourceFiles
    FENE_P.C

\*---------------------------------------------------------------------------*/

#ifndef FENE_P_H
#define FENE_P_H

#include "constitutiveEq.H"
#include "thermoFunction.H" 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class FENE_P Declaration
\*---------------------------------------------------------------------------*/

class FENE_P
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;
        
        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;

            //- Molecular adimensional extensibility
            dimensionedScalar L2_;

            //- Relaxation time
            dimensionedScalar lambda_;
    
    //- Should the equation be solved in tau 
        bool solveInTau_;
        
     //- Should the modified (simplified) formulation be used
     //  (this flag is only valid if solveInTau_ = true)
        bool modifiedForm_;
        
     //- Conformation tensor
        volSymmTensorField A_;
        
     //- Convected parameter needing register for models in tau
     //  varf_ = 1/f for tau models or varf_ = f otherwise  
        volScalarField varf_;
    
     //- Thermofunction for temperature dependence of lambda
        autoPtr<thermoFunction> thermoLambdaPtr_;
        
     //- Thermofunction for temperature dependence of viscosity
        autoPtr<thermoFunction> thermoEtaPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        FENE_P(const FENE_P&);

        //- Disallow default bitwise assignment
        void operator=(const FENE_P&);
        
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
       
       // Return etaS corrected for temperature
       virtual tmp<volScalarField> etaSThermo() const
       {
          return thermoEtaPtr_->createField(etaS_);
       }
       
       // Return etaP corrected for temperature
       virtual tmp<volScalarField> etaPThermo() const
       {
          return thermoEtaPtr_->createField(etaP_);
       }
       
       // Is the model prepared to work in non-isothermal conditions
       virtual bool hasThermo() const
       {
         return true;
       }
 
public:

    //- Runtime type information
    TypeName("FENE-P");


    // Constructors

        //- Construct from components
        FENE_P
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~FENE_P()
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
