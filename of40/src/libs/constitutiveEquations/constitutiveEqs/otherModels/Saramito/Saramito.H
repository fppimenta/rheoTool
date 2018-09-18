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
    Saramito
    
SourceFiles
    Saramito.C

\*---------------------------------------------------------------------------*/

#ifndef Saramito_H
#define Saramito_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class Saramito Declaration
\*---------------------------------------------------------------------------*/

class Saramito
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

            //- Extensibility parameter
            dimensionedScalar epsilon_;
            
            //- Slip parameter
            dimensionedScalar zeta_;

            //- Relaxation time
            dimensionedScalar lambda_;
            
            //- Yield stress
            dimensionedScalar tau0_;
            
            //- Power-law exponent
            dimensionedScalar n_;
            
            //- Consistency index (forced to be = etaP for n = 1)
            dimensionedScalar k_;
            
            //- Vector of solved-for directions in mesh
            vector  dims_;
            
            //- Number of solved-for directions in mesh 
            scalar  nDims;
            
            //- Identity tensor consistent with nDims
            dimensionedSymmTensor Itensor;
            
            //- Write second invariant of deviatoric stresses ?
            Switch writeII_;
            
            //- PTT function to be used
            label funcPTT;
 

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Saramito(const Saramito&);

        //- Disallow default bitwise assignment
        void operator=(const Saramito&);
        
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
    TypeName("Saramito");

    // Constructors

        //- Construct from components
        Saramito
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~Saramito()
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
