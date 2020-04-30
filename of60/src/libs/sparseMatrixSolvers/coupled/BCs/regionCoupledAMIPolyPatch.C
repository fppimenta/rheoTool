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

\*---------------------------------------------------------------------------*/

#include "regionCoupledAMIPolyPatch.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, regionCoupledAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, regionCoupledAMIPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::regionCoupledAMIPolyPatch::regionCoupledAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    regionCoupledBase(static_cast<const polyPatch&>(*this))
{}


Foam::regionCoupledAMIPolyPatch::regionCoupledAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    regionCoupledBase(*this, dict)
{}


Foam::regionCoupledAMIPolyPatch::regionCoupledAMIPolyPatch
(
    const regionCoupledAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    regionCoupledBase(*this, pp)
{}


Foam::regionCoupledAMIPolyPatch::regionCoupledAMIPolyPatch
(
    const regionCoupledAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    regionCoupledBase(*this, pp)
{}


Foam::regionCoupledAMIPolyPatch::regionCoupledAMIPolyPatch
(
    const regionCoupledAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    regionCoupledBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupledAMIPolyPatch::~regionCoupledAMIPolyPatch()
{
    regionCoupledBase::clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupledAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initGeometry(pBufs);
}


void Foam::regionCoupledAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::regionCoupledAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    regionCoupledBase::clearGeom();
}


void Foam::regionCoupledAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::regionCoupledAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    regionCoupledBase::clearGeom();
}


void Foam::regionCoupledAMIPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    regionCoupledBase::write(os);
}


// ************************************************************************* //
