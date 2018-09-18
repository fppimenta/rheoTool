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

#include "ACPotentialFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ACPotentialFvPatchScalarField::
ACPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    offset_(),
    amplitude_(),
    frequency_(),
    phaseDelay_()
{}

Foam::ACPotentialFvPatchScalarField::
ACPotentialFvPatchScalarField
(
    const ACPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    offset_(ptf.offset_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    phaseDelay_(ptf.phaseDelay_)
{}

Foam::ACPotentialFvPatchScalarField::
ACPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    offset_(readScalar(dict.lookup("offset"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    phaseDelay_(readScalar(dict.lookup("phaseDelay")))
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
    
    updateCoeffs();
}
    
Foam::ACPotentialFvPatchScalarField::
ACPotentialFvPatchScalarField
(
    const ACPotentialFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    offset_(tppsf.offset_),
    amplitude_(tppsf.amplitude_),
    frequency_(tppsf.frequency_),
    phaseDelay_(tppsf.phaseDelay_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ACPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
     
    const scalar& t = this->db().time().timeOutputValue(); 
      
     
    operator == (offset_ + amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*t + phaseDelay_)) ;

    fixedValueFvPatchScalarField::updateCoeffs();
}

 
void Foam::ACPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("offset")
        << offset_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("phaseDelay")
        << phaseDelay_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        ACPotentialFvPatchScalarField
    );
}

// ************************************************************************* //
