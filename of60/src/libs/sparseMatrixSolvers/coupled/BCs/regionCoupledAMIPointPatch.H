/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::regionCoupledAMIPointPatch

Description
    Place holder for now.
    
    This file is part of rheoTool.
 
\*---------------------------------------------------------------------------*/

#ifndef regionCoupledAMIPointPatch_H
#define regionCoupledAMIPointPatch_H

#include "facePointPatch.H"
#include "regionCoupledAMIPolyPatch.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionCoupledAMIPointPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledAMIPointPatch
:
    public facePointPatch
{
public:

    //- Runtime type information
    TypeName(regionCoupledAMIPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        regionCoupledAMIPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        ):
             facePointPatch(patch, bm)
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
