/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::constitutiveTwoPhaseMixture

Description
    A two-phase GNF or VE transportModel

SourceFiles
    constitutiveTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef constitutiveTwoPhaseMixture_H
#define constitutiveTwoPhaseMixture_H

#include "constitutiveEq.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class constitutiveTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class constitutiveTwoPhaseMixture
:
    public IOdictionary,
    public twoPhaseMixture
{
protected:

    // Protected data

         // Protected data
        
        autoPtr<constitutiveEq> phase1_;
        autoPtr<constitutiveEq> phase2_;
        
        dimensionedScalar rho1_;
        dimensionedScalar rho2_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;
        
   // Polymeric extra-stress
        volSymmTensorField tauMF_;


public:

    TypeName("constitutiveTwoPhaseMixture");


    // Constructors

        //- Construct from components
        constitutiveTwoPhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~constitutiveTwoPhaseMixture()
    {}


    // Member Functions

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
        
        const volVectorField& U() const
        {
          return U_;
        }

        //- Update the constitutiveEqs
        void correct()
        {
            phase1_->correct(&alpha1_);
            phase2_->correct(&alpha2_);
        }
        
        //- Polymeric extra-stress tensor
        const tmp<volSymmTensorField> tauMF() const
        {
            return tauMF_;
        }
      
        //- Return the multiphase coupling term for the momentum equation
        tmp<fvVectorMatrix> divTauMF(volVectorField& U);
        
        //- Return the multiphase coupling term for the momentum equation, accounting for thermal effects
        tmp<fvVectorMatrix> divTauMFThermo(volVectorField& U);
        
        //- Total extra-stress tensor (polymeric + solvent)
        tmp<volSymmTensorField> tauTotalMF() const;
        
        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
