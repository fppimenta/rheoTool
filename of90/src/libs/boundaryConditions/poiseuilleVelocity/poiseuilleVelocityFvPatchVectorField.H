/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    poiseuilleVelocityFvPatchVectorField

Description
    Analytical velocity profile for laminar pressure-driven flow
    through rectangular (ducts) and circular cross-sections (pipes).
    Valid for Newtonian fluids, although it also fits other types of fluids.
    
    Reference: Frank White (1991) Viscous Fluid Flow, 2nd Ed., McGraw Hill.
    
    Examples:
    
     -> for a 3D duct 
      {
	type            poiseuilleVelocity;
	geometryType    Duct3D;
	uMean           1.; 
	dirA            (0 1 0);
	dirB            (0 0 1);
	wA              1.5;
	wB              1.;
	center          (0 -0.5 0);
	
	value           uniform (0 0 0); // Place-holder for paraview
      }
      
     -> for a 2D duct 
      {
	type            poiseuilleVelocity;
	geometryType    Duct2D;
	uMean           1.; 
	dirA            (0 1 0);
	wA              1.5;
	center          (0 -0.5 0);
	
	value           uniform (0 0 0); // Place-holder for paraview
      }
      
     -> for a pipe 
      {
	type            poiseuilleVelocity;
	geometryType    Circular;
	uMean           0.4;
	R               1.5e-3;
	center          (0 0 0);
	
	value           uniform (0 0 0); // Place-holder for paraview
      }
      
    This BC is part of rheoTool.

SourceFiles
    poiseuilleVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef poiseuilleVelocityFvPatchVectorField_H
#define poiseuilleVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class  Declaration
\*---------------------------------------------------------------------------*/

class poiseuilleVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Options for the patch geommetry
    enum geomOpt
    {
      goDuct3D,
      goDuct2D,
      goCircular
    };

private:

    // Private data
    
    // Center of the velocity profile (not necessarily equal to the geometric center of the patch)
    vector center_;
    
    // First tangential direction for duct3D and the single tangential direction or duct2D
    vector dirA_;
    
    // Second tangential direction for duct3D (normal to dirA_)
    vector dirB_;
    
    // Duct half-width in dirA_
    scalar wA_;
    
    // Duct half-width in dirB_
    scalar wB_;
    
    // (Optional) Max iters for the series used in 3D duct
    label imax_;
    
    // Mean velocity through the patch
    scalar uMean_;
    
    // Radius of circular inlet
    scalar R_;
    
    // Velocity profile (optional - could be directly assigned to operator== upon initEvaluate)
    vectorField U_;
    
    // Names of enum
    static const NamedEnum<geomOpt, 3> geomOptionNames_;
    
    // The geometric option
    geomOpt geomOpts_;
    
    // Functions
    
    // Make the velocity profile for rectangular patch (flow through a 3D duct)
    void computeUDuct3D();
    
    // Make the velocity profile for planar patch (flow between parallel plates)
    void computeUDuct2D();
    
    // Make the velocity profile for circular patch (flow inside pipe)
    void computeUCircular();              

public:

    //- Runtime type information
    TypeName("poiseuilleVelocity");


    // Constructors

        //- Construct from patch and internal field
        poiseuilleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        poiseuilleVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given poiseuilleVelocityFvPatchVectorField
        //  onto a new patch
        poiseuilleVelocityFvPatchVectorField
        (
            const poiseuilleVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct as copy setting internal field reference
        poiseuilleVelocityFvPatchVectorField
        (
            const poiseuilleVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new poiseuilleVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
