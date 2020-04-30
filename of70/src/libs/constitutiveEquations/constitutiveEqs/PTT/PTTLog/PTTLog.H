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
    PTTLog

SourceFiles
    PTTLog.C

\*---------------------------------------------------------------------------*/

#ifndef PTTLog_H
#define PTTLog_H

#include "constitutiveEq.H"
#include "thermoFunction.H" 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class PTTLog Declaration
\*---------------------------------------------------------------------------*/

class PTTLog
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        volSymmTensorField tau_;

        //- Logarithm of the conformation tensor
        volSymmTensorField theta_;
 
        //- Eigenvalues of theta diagonalization (values on the tensor diagonal)
        volTensorField eigVals_;

        //- Eigenvectors of theta diagonalization
        volTensorField eigVecs_;
 
        // Model constants

            //- Density
            dimensionedScalar rho_;

            //- Solvent viscosity
            dimensionedScalar etaS_;

            //- Zero shear rate polymer viscosity
            dimensionedScalar etaP_;
 
            //- Extensibility parameter
            dimensionedScalar epsilon_;
            
            //- Slip parameter
            dimensionedScalar zeta_;

            //- Relaxation time
            dimensionedScalar lambda_;
            
            //- Parameter of the general PTT
            dimensionedScalar alpha_;
            
            //- Parameter of the general PTT
            dimensionedScalar beta_;
            
            //- Relative error while approaching the Mittag-Leffler function
            scalar MLrtol_;
            
            //- Max number of iterations while approaching the Mittag-Leffler function
            int MLmaxIter_;
            
            //- List of values of gamma function to save time
            List<scalar> gammaFunValues_;
            
            //- Thermofunction for temperature dependence of lambda
            autoPtr<thermoFunction> thermoLambdaPtr_;
        
            //- Thermofunction for temperature dependence of viscosity
            autoPtr<thermoFunction> thermoEtaPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        PTTLog(const PTTLog&);

        //- Disallow default bitwise assignment
        void operator=(const PTTLog&);
        
public:
        
    //- Function types
    enum PTTFunctions
    {
      pfLinear,
      pfExpt,
      pfGen
    };       
        
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
       
       // Enum names for PTT function type
       static const NamedEnum<PTTFunctions, 3> PTTFunctionNames_;
       
       //- Type of destruction function
       PTTFunctions PTTFunction_;
 
public:

    //- Runtime type information
    TypeName("PTTLog");

    // Constructors

        //- Construct from components
        PTTLog
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~PTTLog()
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
