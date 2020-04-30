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
    solidThermoModel

Description
    Base class for a set of thermodynamic models for solid regions. All the models
    are temperature-based.
    This file is part of rheoTool.
    
\*---------------------------------------------------------------------------*/

#ifndef solidThermoModel_H
#define solidThermoModel_H

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
                        Class solidThermoModel Declaration
\*---------------------------------------------------------------------------*/

class solidThermoModel
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to mesh
        const fvMesh& mesh_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidThermoModel(const solidThermoModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidThermoModel&);
        
        
protected:

   // Protected Members 
   
    //- Temperature field
    volScalarField T_;
        
    //- Dictionary
    IOdictionary dict_;       
   
public:

    //- Runtime type information
    TypeName("solidThermoModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidThermoModel,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh 
            ),
            (name, mesh)
        );


    // Constructors

        //- Construct from components
        solidThermoModel
        (
            const word& name,
            const fvMesh& mesh
        );


    // Selectors

        //- Return a reference  
        static autoPtr<solidThermoModel> New
        (
            const word& name,
            const fvMesh& mesh
        );


    // Destructor

        virtual ~solidThermoModel()
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
        
        //- Return diffusion number 
        virtual void DiffNo
        (
           const fvMesh& mesh,
           const Time& runTime,
           scalar& maxDiffNo,
           scalar& meanDiffNo
        ) const
        {
          NotImplemented;
        }      
        
        //- Update thermo fields  
        virtual void correct
        ( 
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
