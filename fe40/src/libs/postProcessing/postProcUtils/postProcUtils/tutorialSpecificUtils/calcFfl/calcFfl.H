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
    calcFfl
    
Description
    This util computes and writes the reaction force exerted by the fluid 
    on the piston and shaft. This util is used in the fluidDamper tutorial
    provided in rheoTool.  
    
    This class is part of rheoTool.

SourceFiles
    calcFfl.C

\*---------------------------------------------------------------------------*/


#ifndef calcFfl_H
#define calcFfl_H

#include "ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ppUtils
{

/*---------------------------------------------------------------------------*\
                           Class calcFfl Declaration
\*---------------------------------------------------------------------------*/

class calcFfl
:
    public ppUtil
{
    // Private data

        label shaft_;
        label piston_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcFfl(const calcFfl&);

        //- Disallow default bitwise assignment
        void operator=(const calcFfl&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcFfl");


    // Constructors

        //- Construct from components
        calcFfl
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcFfl()
        {}


    // Member Functions

        //- Correct the viscoelastic stress
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ppUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
