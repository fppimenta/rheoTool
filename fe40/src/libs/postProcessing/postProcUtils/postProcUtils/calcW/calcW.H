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
    calcW
    
Description
    This util computes and writes the width of a splashing drop. The idea is
    to loop over all the rows of the mesh, find the region were alpha1 jump the 
    0.5 barrier and extrapolate to obtain that exact position of the interface
    The code is not generic, depends on the mesh topo (also assuming symmetry)
    and needs to be improved to give smoother profiles. 
  
    This class is part of rheoTool.

SourceFiles
    calcW.C

\*---------------------------------------------------------------------------*/


#ifndef calcW_H
#define calcW_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcW Declaration
\*---------------------------------------------------------------------------*/

class calcW
:
    public ppUtil
{
    // Private data
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcW(const calcW&);

        //- Disallow default bitwise assignment
        void operator=(const calcW&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcW");


    // Constructors

        //- Construct from components
        calcW
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcW()
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
