/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::constitutiveTwoPhaseMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    constitutiveTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/


#ifndef constitutiveTwoPhaseMixture_H
#define constitutiveTwoPhaseMixture_H

#include "constitutiveEq.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class constitutiveTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class constitutiveTwoPhaseMixture
:
    public IOdictionary
{
protected:

    // Protected data
    
        word phase1Name_;
        word phase2Name_;
        
        autoPtr<constitutiveEq> phase1_;
        autoPtr<constitutiveEq> phase2_;
        
        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;
        
        const volScalarField& alpha1_;
        
   // Polymeric extra-stress
        volSymmTensorField tauMF_;
   
   //- Get phase name (backward compatibility)
        word getPhaseName(const word& key) const;
        
public:

    // Constructors

        //- Construct from components
        constitutiveTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word& alpha1Name = "alpha1"
        );


    //- Destructor
    virtual ~constitutiveTwoPhaseMixture()
    {}

    // Member Functions
    
        const word& phase1Name() const
        {
            return phase1Name_;
        }

        const word& phase2Name() const
        {
            return phase2Name_;
        }

        //- Return const-access to phase1 constitutiveModel
        const constitutiveEq& phase1() const
        {
            return phase1_();
        }

       //- Return const-access to phase2 constitutiveModel
        const constitutiveEq& phase2() const
        {
            return phase2_();
        }
        
        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Update the constitutiveEqs
        void correct()
        {
            phase1_->correct();
            phase2_->correct();
        }
        
        //- Polymeric extra-stress tensor
        const tmp<volSymmTensorField> tauMF() const
        {
            return tauMF_;
        }
      
        //- Return the multiphase coupling term for the momentum equation
        tmp<fvVectorMatrix> divTauMF(volVectorField& U);
        
        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
