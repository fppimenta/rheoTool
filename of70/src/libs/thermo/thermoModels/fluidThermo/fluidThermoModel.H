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
    fluidThermoModel

Description
    Base class for a set of thermodynamic models for fluid regions. All the models
    are temperature-based.
    This file is part of rheoTool.
\*---------------------------------------------------------------------------*/

#ifndef fluidThermoModel_H
#define fluidThermoModel_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "fvOptions.H"
#include "NamedEnum.H"
#include "runTimeSelectionTables.H"
#include "extrapolatedCalculatedFvPatchField.H"

#include "coupledSolver.H" 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidThermoModel Declaration
\*---------------------------------------------------------------------------*/

class fluidThermoModel
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        fluidThermoModel(const fluidThermoModel&);

        //- Disallow default bitwise assignment
        void operator=(const fluidThermoModel&);
        
        
protected:

   // Protected Members 
   
    //- Temperature field
    volScalarField T_;
        
    //- Dictionary
    IOdictionary dict_;       
   
public:

    //- Runtime type information
    TypeName("fluidThermoModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidThermoModel,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh 
            ),
            (name, mesh)
        );


    // Constructors

        //- Construct from components
        fluidThermoModel
        (
            const word& name,
            const fvMesh& mesh
        );


    // Selectors

        //- Return a reference  
        static autoPtr<fluidThermoModel> New
        (
            const word& name,
            const fvMesh& mesh
        );


    // Destructor

        virtual ~fluidThermoModel()
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
            return T_;
        }
        
        //- Return dictionary
        const dictionary& dict() const
        {
            return dict_;
        }
        
        //- FIXME
        virtual tmp<volScalarField> rhok() const
        {
          NotImplemented;
          return tmp<volScalarField>();
        }
        
        //- Is gravity used in the model?
        virtual bool hasGravity() const = 0;
        
        //- Return g.h for models that have hasGravity() = true
        virtual tmp<volScalarField> gh() const
        {
          NotImplemented;
          return tmp<volScalarField>();
        }
        
        //- Return g.hf for models that have hasGravity() = true
        virtual tmp<surfaceScalarField> ghf() const
        {
          NotImplemented;
          return tmp<surfaceScalarField>();
        }   
        
        //- Update temperature  
        virtual void correct
        (
          const volVectorField& U,
          const surfaceScalarField& phi,
          const volSymmTensorField& tau,
          autoPtr<coupledSolver>& cpsT,
          fv::options&,
          int nNonOrtoC
        ) = 0;
       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
