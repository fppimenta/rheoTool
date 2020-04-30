/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  7.0                                   |
|   \\  /    A nd           | Website:  https://openfoam.org                  |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    elecNorth
    {
	type            zeroGradient;      
    }
   
    elecSouth
    {
	type            zeroGradient;        
    }
    
    cylinder
    {
	type            zeroGradient; 
    }

    "wall.*"
    {
	type            fixedValue;
	value           uniform  0;    
    }
    
    frontAndBack
    {
	type            empty;     
    }
}

// ************************************************************************* //
