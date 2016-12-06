/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "HBprofileFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

HBprofileFvPatchVectorField::HBprofileFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    R_(1),
    uIn_(0, 0, 0)
{}


HBprofileFvPatchVectorField::HBprofileFvPatchVectorField
(
    const HBprofileFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    R_(ptf.R_),
    uIn_(ptf.uIn_)
{}


HBprofileFvPatchVectorField::HBprofileFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    R_(readScalar(dict.lookup("R"))),
    uIn_(dict.lookup("uIn"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


HBprofileFvPatchVectorField::HBprofileFvPatchVectorField
(
    const HBprofileFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    R_(fcvpvf.R_),
    uIn_(fcvpvf.uIn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void HBprofileFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const dictionary& constDict = db().lookupObject<IOdictionary>("constitutiveProperties");
    dimensionedScalar n_(constDict.subDict("parameters").lookup("n"));
    
    scalar n(n_.value());

    const vectorField& x = patch().Cf();

    scalarField r = mag(x);
      
    vectorField::operator=( ((3.*n+1.)/(n+1.))*( 1. - pow(r/R_, (n+1.)/n)  ) * uIn_ );
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void HBprofileFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("R")
        << R_ << token::END_STATEMENT << nl;
    os.writeKeyword("uIn")
        << uIn_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, HBprofileFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
