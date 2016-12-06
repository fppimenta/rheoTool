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

#include "uCosFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

uCosFvPatchVectorField::uCosFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    tlim_(1),
    fac_(1),
    uav_(0, 0, 0),
    dirN_(1, 0, 0)
{}


uCosFvPatchVectorField::uCosFvPatchVectorField
(
    const uCosFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tlim_(ptf.tlim_),
    fac_(ptf.fac_),
    uav_(ptf.uav_),
    dirN_(ptf.dirN_)
{}


uCosFvPatchVectorField::uCosFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    tlim_(readScalar(dict.lookup("tlim"))),
    fac_(readScalar(dict.lookup("fac"))),
    uav_(dict.lookup("uav")),
    dirN_(dict.lookup("dirN"))
{
  updateCoeffs();
}


uCosFvPatchVectorField::uCosFvPatchVectorField
(
    const uCosFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    tlim_(fcvpvf.tlim_),
    fac_(fcvpvf.fac_),
    uav_(fcvpvf.uav_),
    dirN_(fcvpvf.dirN_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void uCosFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
     
    vector Ut = uav_;
 
    if (t<=tlim_)
    {
        Ut = ( ( (1 - Foam::cos( 3.1415926535897932 * t) ) / fac_) * dirN_);
    }
      
    vectorField::operator=(Ut);
    
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void uCosFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("tlim")
        << tlim_ << token::END_STATEMENT << nl;
    os.writeKeyword("fac")
        << fac_ << token::END_STATEMENT << nl;
    os.writeKeyword("uav")
        << uav_ << token::END_STATEMENT << nl;
    os.writeKeyword("dirN")
        << dirN_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, uCosFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
