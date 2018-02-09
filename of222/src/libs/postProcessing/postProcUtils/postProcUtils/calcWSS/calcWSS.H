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
    calcWSS
    
Description
    This util computes and writes the wall shear stress magnitude for any GNF
    or viscoelastic model. The interior field is kept as zero. The user must
    provide the following entries (I guess names are suggestive):
    
     - isViscoelastic: true/false
     - includeSolventStresses: true/false
     - includePolymericStresses: true/false
     
    Note that the last two entries have no meaning if isViscoelastic = false 
    (for GNF there are no polymeric stresses). 
    
    This class is part of rheoTool.

SourceFiles
    calcWSS.C

\*---------------------------------------------------------------------------*/

#ifndef calcWSS_H
#define calcWSS_H

#include "../ppUtil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class calcWSS Declaration
\*---------------------------------------------------------------------------*/

class calcWSS
:
    public ppUtil
{
    // Private data
    
        //- Is it a viscoelastic fluid
        Switch  isVE_;
        
        //- Account for Newtonian stresses
        Switch incSolv_;
        
        //- Account for polymeric stresses
        Switch  incPoly_;
        
        //- Solvent viscosity for the case of VE
        dimensionedScalar etaSWW_;
        
        //- Wall shear stress magnitude
        volScalarField WSSmag_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        calcWSS(const calcWSS&);

        //- Disallow default bitwise assignment
        void operator=(const calcWSS&);
        

public:

    //- Runtime type information
    TypeName("calcWSS");


    // Constructors

        //- Construct from components
        calcWSS
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U
        );

    // Destructor

        virtual ~calcWSS()
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
