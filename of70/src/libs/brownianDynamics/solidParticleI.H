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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::solidParticle::trackingData::trackingData
(
    const solidParticleCloud& spc,
    const bool& includeDrag
)
:
    particle::trackingData(spc),
    includeDrag_(includeDrag)
{}


inline Foam::solidParticle::solidParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const scalar d,
    const vector& U
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    U_(U)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
inline const bool&
Foam::solidParticle::trackingData::includeDrag() const
{
    return includeDrag_;
}


inline const Foam::vector& Foam::solidParticle::U() const
{
    return U_;
}

inline const Foam::Field<Foam::label>& Foam::solidParticle::ids() const
{
    return ids_;
}

inline const Foam::label& Foam::solidParticle::molcID() const
{
    return molcID_;
}

// ************************************************************************* //
