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
    fluidThermoMixtureModel

Description
    Base class for a set of thermodynamic models for two-phase incompressible
    flows. All the models are temperature-based. The same temperature field is
    shared by the two phases. 
    This file is part of rheoTool.
    
\*---------------------------------------------------------------------------*/

#ifndef fluidThermoMixtureModel_H
#define fluidThermoMixtureModel_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "NamedEnum.H"
#include "runTimeSelectionTables.H"
#include "extrapolatedCalculatedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidThermoMixtureModel Declaration
\*---------------------------------------------------------------------------*/

class fluidThermoMixtureModel
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        fluidThermoMixtureModel(const fluidThermoMixtureModel&);

        //- Disallow default bitwise assignment
        void operator=(const fluidThermoMixtureModel&);
        
        
protected:

   // Protected Members 
   
    //- Temperature field
    volScalarField T_;
        
    //- Dictionary
    IOdictionary dict_;       
   
public:

    //- Runtime type information
    TypeName("fluidThermoMixtureModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidThermoMixtureModel,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh, 
                const word& phase1,
                const word& phase2
            ),
            (name, mesh, phase1, phase2)
        );


    // Constructors

        //- Construct from components
        fluidThermoMixtureModel
        (
            const word& name,
            const fvMesh& mesh,
            const word& phase1,
            const word& phase2
        );


    // Selectors

        //- Return a reference  
        static autoPtr<fluidThermoMixtureModel> New
        (
            const word& name,
            const fvMesh& mesh,
            const word& phase1,
            const word& phase2
        );


    // Destructor

        virtual ~fluidThermoMixtureModel()
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
        
        //- Return specific heat
        virtual tmp<volScalarField> cp
        ( 
          const volScalarField& alpha1 
        ) const
        {
          NotImplemented;
          return tmp<volScalarField>();
        };        
        
        //- Return dictionary
        const dictionary& dict() const
        {
            return dict_;
        }
               
        //- Update thermo fields 
        virtual void correct
        (
          const volVectorField& U,
          const surfaceScalarField& rhoPhi,
          const volSymmTensorField& tau,
          const volScalarField& alpha1,
          const volScalarField& rho,
          fvModels&,
          fvConstraints&,
          int nNonOrtoC
        ) = 0;
       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
