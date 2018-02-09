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
    calcKineticE
    
Description
    This util computes and writes both the volume-averaged kinetic and
    elastic energy in a domain composed of identical cells. It is used
    in the cavity tutorial of rheoTool.  
    
    This class is part of rheoTool.

SourceFiles
    calcKineticE.C

\*---------------------------------------------------------------------------*/

#ifndef calcKineticE_H
#define calcKineticE_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcKineticE Declaration
\*---------------------------------------------------------------------------*/

class calcKineticE
:
    public ppUtil
{
    // Private data

    

    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcKineticE(const calcKineticE&);

        //- Disallow default bitwise assignment
        void operator=(const calcKineticE&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcKineticE");


    // Constructors

        //- Construct from components
        calcKineticE
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcKineticE()
        {}


    // Member Functions

        //- Correct the viscoelastic stress
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
