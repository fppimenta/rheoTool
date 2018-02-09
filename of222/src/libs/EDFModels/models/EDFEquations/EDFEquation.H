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

 

\*---------------------------------------------------------------------------*/

#ifndef EDFEquation_H
#define EDFEquation_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class  Declaration
\*---------------------------------------------------------------------------*/

class EDFEquation
{
private:

    // Private data

        //- Name
        const word name_;
 
        //- Reference to face flux field
        const surfaceScalarField& phi_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        EDFEquation(const EDFEquation&);

        //- Disallow default bitwise assignment
        void operator=(const EDFEquation&);

public:

    //- Static members
    
      //- Vacuum permittivity
        static const dimensionedScalar epsilonK_;
        
      //- Avogrado's number
        static const dimensionedScalar AK_;
        
      //- Elementary charge
        static const dimensionedScalar eK_;
        
      //- Boltzmann constant
        static const dimensionedScalar kbK_;
        
      //- Faraday constant
        static const dimensionedScalar FK_;

protected:

     //- Will check if phiE_ exists in the time Dir: true if yes, false otherwise
        bool checkForPhiE(const word& name, const surfaceScalarField& phi);
        
public:

    //- Runtime type information
    TypeName("EDFEquation");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            EDFEquation,
            dictionary,
            (
                const word& name,
                const surfaceScalarField& phi,
                const dictionary& dict
            ),
            (name, phi, dict)
        );


    // Constructors

        //- Construct from components
        EDFEquation
        (
            const word& name,
            const surfaceScalarField& phi
        );


    // Selectors

        //- Return a reference to the selected viscoelastic law
        static autoPtr<EDFEquation> New
        (
            const word& name,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~EDFEquation()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the face flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return the electric force for the momentum equation
        virtual tmp<volVectorField> Fe() const = 0;

        //- Update the electric force
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
