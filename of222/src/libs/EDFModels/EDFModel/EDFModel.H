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

#ifndef EDFModel_H
#define EDFModel_H

#include "EDFEquation.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class EDFModel Declaration
\*---------------------------------------------------------------------------*/

class EDFModel
:
    public IOdictionary
{
private:

    // Private data

        //- Viscoelastic law
        autoPtr<EDFEquation> modelPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        EDFModel(const EDFModel&);

        //- Disallow default bitwise assignment
        void operator=(const EDFModel&);


public:

    //- Runtime type information
    TypeName("EDFModel");


    // Constructors

        //- Construct from components
        EDFModel
        (            
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~EDFModel()
        {}


    // Member Functions

        //- Return the electric force for the momentum equation
        virtual tmp<volVectorField> Fe() const;

        //- Update the electric force
        virtual void correct();

        //- Read the dictionary
        virtual bool read();
         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
