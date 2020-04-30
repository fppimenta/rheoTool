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
    Constant

Description
    Implements an isothermal function. Any parameter given to this function
    remains unmodified.
    This file is part of rheoTool.

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "thermoFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermoFunctions
{ 

/*---------------------------------------------------------------------------*\
                           Class Constant Declaration
\*---------------------------------------------------------------------------*/

class Constant
:
    public thermoFunction
{
    // Private data

         // Model Constant


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Constant(const Constant&);

        //- Disallow default bitwise assignment
        void operator=(const Constant&);

protected:
   
    
      
public:

    //- Runtime type information
    TypeName("Constant");

    // Constructors

        //- Construct from components
        Constant
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~Constant()
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
           
           if ( isType<Constant>(tf) )
           {
             fl = true;     
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
