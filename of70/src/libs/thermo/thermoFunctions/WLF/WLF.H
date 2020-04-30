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
    WLF

Description
    Implements Williams–Landel–Ferry (WLF) function.
    This file is part of rheoTool.

\*---------------------------------------------------------------------------*/

#ifndef WLF_H
#define WLF_H

#include "thermoFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermoFunctions
{ 

/*---------------------------------------------------------------------------*\
                           Class WLF Declaration
\*---------------------------------------------------------------------------*/

class WLF
:
    public thermoFunction
{
    // Private data

         // Model WLFs

           //- Numerator
           dimensionedScalar c1_;
           
           //- Denominator
           dimensionedScalar c2_;

           //- Reference temperature
           dimensionedScalar T0_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        WLF(const WLF&);

        //- Disallow default bitwise assignment
        void operator=(const WLF&);

protected:
   
    
      
public:

    //- Runtime type information
    TypeName("WLF");

    // Constructors

        //- Construct from components
        WLF
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~WLF()
        {}

    // Member Functions

        //- Given the field for a physical parameter, make it equal to the 
        // dimensionedScalar multiplied by the scaling function
        virtual void multiply
        (
          volScalarField&,
          const dimensionedScalar&
        ) const;
        
        //- Multiply the given field by the scaling function
        virtual void multiply
        (
          volScalarField&
        ) const;
        
        //- Create a field with its internal and boundary fields set equal to the 
        // dimensionedScalar passed and multiply it by the scaling function
        virtual tmp<volScalarField> createField
        (
          const dimensionedScalar&
        ) const;
        
        //- Compare two thermoFunctions for equality 
        virtual bool isItEqualTo
        (
          const thermoFunction& tf
        ) const
        {
           bool fl(false);
           
           if ( isType<WLF>(tf) )
           {
             const WLF& ar = refCast<const WLF>(tf);
             fl = (this->c1_.value() == ar.c1_.value() && this->c2_.value() == ar.c2_.value() && this->T0_.value() == ar.T0_.value());     
           } 
                       
           return fl;     
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermoFunctions 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
