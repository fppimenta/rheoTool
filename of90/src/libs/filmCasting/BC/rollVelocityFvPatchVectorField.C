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

#include "rollVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::rollVelocityFvPatchVectorField::
rollVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}

Foam::rollVelocityFvPatchVectorField::
rollVelocityFvPatchVectorField
(
    const rollVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    hT_(ptf.hT_),
    uStart_(ptf.uStart_),
    uEnd_(ptf.uEnd_),
    isZTBC_(ptf.isZTBC_)  
{}

Foam::rollVelocityFvPatchVectorField::
rollVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    hT_(readScalar(dict.lookup("hT"))),
    uStart_(readScalar(dict.lookup("uStart"))),
    uEnd_(readScalar(dict.lookup("uEnd"))),
    isZTBC_(false)
{
 if (dict.found("value"))
 {
   fvPatchField<vector>::operator=
   (
     vectorField("value", dict, p.size())
   );
 }
 else
 {
   fvPatchField<vector>::operator=
   (
     vectorField(p.size(), vector(uStart_,0.,0.))
   );
 }    
}
    
Foam::rollVelocityFvPatchVectorField::
rollVelocityFvPatchVectorField
(
    const rollVelocityFvPatchVectorField& tppsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tppsf, iF),
    hT_(tppsf.hT_),
    uStart_(tppsf.uStart_),
    uEnd_(tppsf.uEnd_),
    isZTBC_(tppsf.isZTBC_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rollVelocityFvPatchVectorField::updateCoeffs()
{
   if (updated())
   {
     return;
   }  

   const scalar& t = this->db().time().timeOutputValue();
     
   const dictionary& cttDict = db().lookupObject<IOdictionary>("constitutiveProperties");
   isZTBC_ = readBool(cttDict.subDict("filmProperties").lookup("isTanZeroTraction"));
   
   scalar dU = uEnd_ - uStart_;
   scalar uNow;
   
   if (t<hT_)
   {
     uNow = uStart_ + dU*(1. - Foam::cos(M_PI*t/hT_))/2.;
   }
   else
   {
     uNow = uEnd_;
   }

   if (!isZTBC_)
   {
     vectorField::operator= (vector(uNow, 0, 0));
   }
   else
   {
     // This is rigorous for inelastic fluids, but an approximation for viscoelastic fluids
     vectorField uf(this->patchInternalField());
     forAll(uf, i)
       uf[i].x() = uNow;

     vectorField::operator=(uf);
   }     
   
   fixedValueFvPatchVectorField::updateCoeffs();
}

 
void Foam::rollVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "hT", hT_);
    writeEntry(os, "uStart", uStart_);
    writeEntry(os, "uEnd", uEnd_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        rollVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
