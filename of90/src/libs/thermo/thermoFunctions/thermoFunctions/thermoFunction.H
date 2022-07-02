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
    thermoFunction

Description
    Base class for a set of functions to scale physical parameters (viscosity
    and relaxation time) using the temperature field. Temperature is not owned
    by the class, but initialized based on a lookup to db. Must also work with 
    isothermal solvers.
    If no type is specified by the user for the thermo function, the isothermal
    one is automatically selected (this is needed to ensure compatibility with
    isothermal solvers). 
    This file is part of rheoTool.
    
\*---------------------------------------------------------------------------*/

#ifndef thermoFunction_H
#define thermoFunction_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "NamedEnum.H"
#include "runTimeSelectionTables.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class thermoFunction Declaration
\*---------------------------------------------------------------------------*/

class thermoFunction
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermoFunction(const thermoFunction&);

        //- Disallow default bitwise assignment
        void operator=(const thermoFunction&);
        
        
protected:

   // Protected Members
   
   //- Pointer to temperature field
        const volScalarField* T_;
        
   
public:

    //- Runtime type information
    TypeName("thermoFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermoFunction,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Constructors

        //- Construct from components
        thermoFunction
        (
            const word& name,
            const fvMesh& mesh
        );


    // Selectors

        //- Return a reference  
        static autoPtr<thermoFunction> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~thermoFunction()
        {}
  
    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the temperature
        const volScalarField& T() const
        {
            return *T_;
        }
        
        //- Compare two thermoFunctions for equality 
        virtual bool isItEqualTo
        (
          const thermoFunction&
        ) const = 0;
        
        //- Given the field for a physical parameter, make it equal to the 
        // dimensionedScalar multiplied by the scaling function
        virtual void multiply
        (
          volScalarField&,
          const dimensionedScalar&
        ) const = 0;
        
        //- Multiply the given field by the scaling function
        virtual void multiply
        (
          volScalarField&
        ) const = 0;
        
        //- Create a field with its internal and boundary fields set equal to the 
        // dimensionedScalar passed and multiply it by the scaling function
        virtual tmp<volScalarField> createField
        (
          const dimensionedScalar&
        ) const = 0;
         
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
