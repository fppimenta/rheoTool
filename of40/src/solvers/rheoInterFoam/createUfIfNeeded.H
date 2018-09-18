/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
Global
    createUf
Description
    Creates and initialises the velocity field Uf if required.
\*---------------------------------------------------------------------------*/

// This is an hack to check if mesh is changing. At this point, we can't
// check it with mesh.changing() because update() was still no called.
// Thus read the dict and consider that the mesh changes if its type is not
// staticFvMesh. 

IOdictionary dyMDict
(
 IOobject
 (
   "dynamicMeshDict",
   runTime.constant(),
   mesh,
   IOobject::MUST_READ,
   IOobject::NO_WRITE,
   false
 )
);

word  dyMType = dyMDict.lookup("dynamicFvMesh");

autoPtr<surfaceVectorField>  Uf
(
    dyMType == "staticFvMesh"
    ? 
    NULL
    :
    new surfaceVectorField
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U)
    )
);

// ************************************************************************* //
