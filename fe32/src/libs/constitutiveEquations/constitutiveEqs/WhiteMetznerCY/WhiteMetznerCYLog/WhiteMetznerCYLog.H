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
    WhiteMetznerCYLog

SourceFiles
    WhiteMetznerCYLog.C

\*---------------------------------------------------------------------------*/

#ifndef WhiteMetznerCYLog_H
#define WhiteMetznerCYLog_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WhiteMetznerCYLog Declaration
\*---------------------------------------------------------------------------*/

class WhiteMetznerCYLog
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

            //- Relaxation time
            dimensionedScalar lambda_;
            
            //- Fitted parameter for lambda
            dimensionedScalar m_;
        
             //- Fitted parameter for etaP
            dimensionedScalar n_;

            //- Fitted parameter for etaP
            dimensionedScalar K_;

            //- Fitted parameter for lambda
            dimensionedScalar L_;

            //- Fitted parameter for etaP
            dimensionedScalar a_;

            //- Fitted parameter for lambda
            dimensionedScalar b_;

        //- Include or not the U-tau coupling
        bool uTauCoupling_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WhiteMetznerCYLog(const WhiteMetznerCYLog&);

        //- Disallow default bitwise assignment
        void operator=(const WhiteMetznerCYLog&);

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
       
       //- Return true if there is coupling Utau
       virtual bool coupledUtau() const
       {
          return uTauCoupling_;
       }
       
public:

    //- Runtime type information
    TypeName("WhiteMetznerCYLog");

    // Constructors

        //- Construct from components
        WhiteMetznerCYLog
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~WhiteMetznerCYLog()
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
