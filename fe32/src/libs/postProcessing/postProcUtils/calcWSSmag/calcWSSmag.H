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
    calcWSSmag
    
Description
    This util computes and writes the wall shear stress magnitude for any GNF
    model. The interior field is kept as zero. The util is used in the aneurysm
    tutorial of rheoTool.   
    
    This class is part of rheoTool.

SourceFiles
    calcWSSmag.C

\*---------------------------------------------------------------------------*/

#ifndef calcWSSmag_H
#define calcWSSmag_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcWSSmag Declaration
\*---------------------------------------------------------------------------*/

class calcWSSmag
:
    public ppUtil
{
    // Private data

        volScalarField WSSmag_;
    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcWSSmag(const calcWSSmag&);

        //- Disallow default bitwise assignment
        void operator=(const calcWSSmag&);
        

public:

    //- Runtime type information
    TypeName("calcWSSmag");


    // Constructors

        //- Construct from components
        calcWSSmag
        (
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcWSSmag()
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
