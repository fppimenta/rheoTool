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
    calcWi0
    
Description
    This util computes and writes the local Weissenberg number at point 
    (0,0,0.5) of the mesh. This corresponds to the geometrical center of
    the cross-slot used in the respective tutorial. The point coordinates
    can be changed in the constructor.   
    
    This class is part of rheoTool.

SourceFiles
    calcWi0.C

\*---------------------------------------------------------------------------*/

#ifndef calcWi0_H
#define calcWi0_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcWi0 Declaration
\*---------------------------------------------------------------------------*/

class calcWi0
:
    public ppUtil
{
    // Private data

        label cellC_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcWi0(const calcWi0&);

        //- Disallow default bitwise assignment
        void operator=(const calcWi0&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcWi0");


    // Constructors

        //- Construct from components
        calcWi0
        (
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcWi0()
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
