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
    
\*---------------------------------------------------------------------------*/

#include "externalForcingInterp.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(externalForcingInterp, 0);
    defineRunTimeSelectionTable(externalForcingInterp, dicteFI);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

externalForcingInterp::externalForcingInterp
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
mesh_(mesh)
{}

autoPtr<externalForcingInterp> externalForcingInterp::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word typeName = dict.subDict("externalFlow").lookup("interpolation");

    Info<< "Selecting interpolation type: " << typeName << endl;

    dicteFIConstructorTable::iterator cstrIter =
        dicteFIConstructorTablePtr_->find(typeName);

    if (cstrIter == dicteFIConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "externalForcingInterp::New(const fvMesh& mesh, const dictionary& dict)"
        )   << "Unknown interpolation type " << typeName
            << endl << endl
            << "Valid interpolation types are :" << endl
            << dicteFIConstructorTablePtr_->toc()
            << exit(FatalError);
    }
 
    return autoPtr<externalForcingInterp>(cstrIter()(mesh, dict));
}


} //End namespace

// ************************************************************************* //
