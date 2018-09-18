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
    constitutiveModel

Description
    Base class which can be directly called by applications to use any of 
    the models derived from class constitutivEq.

    Based on the libviscoelasticTransportModels library written by Jovani
    L. Favero. 

\*---------------------------------------------------------------------------*/

#ifndef constitutiveModel_H
#define constitutiveModel_H

#include "constitutiveEq.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class constitutiveModel Declaration
\*---------------------------------------------------------------------------*/

class constitutiveModel
:
    public IOdictionary
{
private:

    // Private data

        //- Constitutive Eq
        autoPtr<constitutiveEq> eqPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constitutiveModel(const constitutiveModel&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveModel&);


public:

    //- Runtime type information
    TypeName("constitutiveModel");


    // Constructors

        
        //- Construct from components
        constitutiveModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    // Destructor

        virtual ~constitutiveModel()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const;
        
        //- Return the total stress tensor (solvent + polymeric)
        virtual tmp<volSymmTensorField> tauTotal() const;
        
        //- Return the density
        virtual const dimensionedScalar rho() const;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;
        
        //- Is GNF or VE
        virtual bool isGNF() const;
       
        //- Correct the viscoelastic stress
        virtual void correct();

        //- Read constitutiveProperties dictionary
        virtual bool read();

         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
