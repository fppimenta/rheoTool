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
    calcCd
    
Description
    This util computes and writes the drag coefficient of patch named 'cylinder',
    whose axis is z aligned and for a flow in the x direction. This util is
    used in the cylinder tutorial provided in rheoTool.  
    
    This class is part of rheoTool.

SourceFiles
    calcCd.C

\*---------------------------------------------------------------------------*/


#ifndef calcCd_H
#define calcCd_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcCd Declaration
\*---------------------------------------------------------------------------*/

class calcCd
:
    public ppUtil
{
    // Private data

        label cyl_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcCd(const calcCd&);

        //- Disallow default bitwise assignment
        void operator=(const calcCd&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcCd");


    // Constructors

        //- Construct from components
        calcCd
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcCd()
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
