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

#include "regionCoupledAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledAMIFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, regionCoupledAMIFvPatch, polyPatch);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::regionCoupledAMIFvPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::vectorField> Foam::regionCoupledAMIFvPatch::delta() const
{
    // Use patch-normal delta for all non-coupled BCs
    const vectorField nHat(nf());
    return nHat*(nHat & (Cf() - Cn()));
}



Foam::tmp<Foam::labelField> Foam::regionCoupledAMIFvPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::regionCoupledAMIFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    if (neighbFvPatch().sameRegion())
    {
        return neighbFvPatch().patchInternalField(iF);
    }
    else
    {
        return tmp<labelField>(new labelField(iF.size(), 0));

    }

    return tmp<labelField>(nullptr);
}


// ************************************************************************* //
