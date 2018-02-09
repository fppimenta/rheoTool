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

#include "inducedPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "EDFEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inducedPotentialFvPatchScalarField::
inducedPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    psiF_()
{}

Foam::inducedPotentialFvPatchScalarField::
inducedPotentialFvPatchScalarField
(
    const inducedPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    psiF_(ptf.psiF_().clone().ptr())
{}

Foam::inducedPotentialFvPatchScalarField::
inducedPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    psiF_(DataEntry<scalar>::New("psiF", dict))
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}
    
Foam::inducedPotentialFvPatchScalarField::
inducedPotentialFvPatchScalarField
(
    const inducedPotentialFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    psiF_(tppsf.psiF_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inducedPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const scalar t = this->db().time().timeOutputValue();
    const scalar psiF = psiF_->value(t);

    const fvPatchField<scalar>& phiE =
        patch().lookupPatchField<volScalarField, scalar>("phiE");

    // Face-averaged external potential on the patch
    scalar avPsi( gSum( patch().magSf() * phiE) / gSum(patch().magSf()) );
 
    operator == (-phiE + psiF + avPsi) ;

    fixedValueFvPatchScalarField::updateCoeffs();
}

 
void Foam::inducedPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    psiF_->writeData(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inducedPotentialFvPatchScalarField
    );
}

// ************************************************************************* //
