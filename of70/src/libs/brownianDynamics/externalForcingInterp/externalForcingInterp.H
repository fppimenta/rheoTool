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
    externalForcingInterp
    
Description
    Base class for external forcing interpolation.
   
SourceFiles
    externalForcingInterp.C

\*---------------------------------------------------------------------------*/

#ifndef externalForcingInterp_H
#define externalForcingInterp_H

#include "volFields.H"
#include "solidParticle.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class externalForcingInterp Declaration
\*---------------------------------------------------------------------------*/

class externalForcingInterp
{
private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        externalForcingInterp(const externalForcingInterp&);

        //- Disallow default bitwise assignment
        void operator=(const externalForcingInterp&);
        
protected:

    // Protected Member Functions
       //- Reference to fvMesh
       const fvMesh& mesh_;  
        
public:

    //- Runtime type information
    TypeName("externalForcingInterp");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            externalForcingInterp,
            dicteFI,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        externalForcingInterp
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference  
        static autoPtr<externalForcingInterp> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~externalForcingInterp() 
        {}


    // Member Functions

        //- Update
        virtual void update() = 0;
        
        //- Interpolate
        virtual const vector interpolate(const solidParticle&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
