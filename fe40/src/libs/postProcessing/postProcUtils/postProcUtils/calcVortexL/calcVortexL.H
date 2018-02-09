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
    calcVortexL
    
Description
    This util computes and writes the coordinates of all the points potentially
    delimiting a vortex. This is done for the two corner vortices and the two
    lip vortices potentially present in the contraction tutorial of rheoTool.
    
    Warning: not prepared to be used in parallel runs!   
    
    This class is part of rheoTool.

SourceFiles
    calcVortexL.C

\*---------------------------------------------------------------------------*/

#ifndef calcVortexL_H
#define calcVortexL_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcVortexL Declaration
\*---------------------------------------------------------------------------*/

class calcVortexL
:
    public ppUtil
{
    // Private data

        label vortTop_;
        label vortDown_;
        label lipTop_;
        label lipDown_;
        
      // OStream
        autoPtr<OFstream> outS1_;
        autoPtr<OFstream> outS2_;
        autoPtr<OFstream> outS3_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcVortexL(const calcVortexL&);

        //- Disallow default bitwise assignment
        void operator=(const calcVortexL&);
        
        //- Create file to write
        void createFile();


public:

    //- Runtime type information
    TypeName("calcVortexL");


    // Constructors

        //- Construct from components
        calcVortexL
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcVortexL()
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
