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
    Gradient
    
Description
    The external forcing is interpolated cellwise using its cell-centered
    value and gradient. The gradient discretization scheme can be defined
    by the user.

SourceFiles
    Gradient.C

\*---------------------------------------------------------------------------*/

#ifndef Gradient_H
#define Gradient_H

#include "externalForcingInterp.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace externalForcingInterps
{

/*---------------------------------------------------------------------------*\
                           Class Gradient Declaration
\*---------------------------------------------------------------------------*/

class Gradient
:
    public externalForcingInterp
{
    // Private data
    
        //- Is electrophoresis active
        bool epActive_;
        
        //- Electrophoretic mobility
        scalar epMobility_;

        //- Is the external numerical forcing steady in time (analytical forcing is always steady)
        bool frozenFlow_;
        
        //- Continuum field
        volVectorField Uex_;
        
        //- Gradient of continuum field
        volTensorField gradUex_;
        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        Gradient(const Gradient&);

        //- Disallow default bitwise assignment
        void operator=(const Gradient&);
        
        //- Update Uex 
        void updateUex();

public:

    //- Runtime type information
    TypeName("Gradient");


    // Constructors

        //- Construct from components
        Gradient
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Destructor

        virtual ~Gradient()
        {}


    // Member Functions

        //- Update
        virtual void update();
        
        //- Interpolate
        virtual inline const vector interpolate(const solidParticle&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace externalForcingInterps
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GradientI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
