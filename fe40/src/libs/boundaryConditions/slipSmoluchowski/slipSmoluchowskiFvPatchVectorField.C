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

#include "slipSmoluchowskiFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"
#include "EDFEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slipSmoluchowskiFvPatchVectorField::
slipSmoluchowskiFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    elecM_()
{}

Foam::slipSmoluchowskiFvPatchVectorField::
slipSmoluchowskiFvPatchVectorField
(
    const slipSmoluchowskiFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    elecM_(ptf.elecM_)   
{}

Foam::slipSmoluchowskiFvPatchVectorField::
slipSmoluchowskiFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    elecM_(readScalar(dict.lookup("elecMobility")))
{
    fvPatchField<vector>::operator=
    (
        vectorField("value", dict, p.size())
    );
    
}
    
Foam::slipSmoluchowskiFvPatchVectorField::
slipSmoluchowskiFvPatchVectorField
(
    const slipSmoluchowskiFvPatchVectorField& tppsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tppsf, iF),
    elecM_(tppsf.elecM_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slipSmoluchowskiFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const volScalarField& phiE_ =
        db().lookupObject<volScalarField>("phiE");
        
    volVectorField Ef(-fvc::grad(phiE_));
     
    vectorField::operator=( elecM_ * Ef.boundaryField()[patch().index()] );
         
    fixedValueFvPatchVectorField::updateCoeffs();
}

 
void Foam::slipSmoluchowskiFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("elecMobility")
        << elecM_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        slipSmoluchowskiFvPatchVectorField
    );
}

// ************************************************************************* //
