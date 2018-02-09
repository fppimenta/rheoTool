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

#include "slipSigmaDependentFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"
#include "EDFEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slipSigmaDependentFvPatchVectorField::
slipSigmaDependentFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    elecM_(),
    sigma0_(),
    m_(0.)
{}

Foam::slipSigmaDependentFvPatchVectorField::
slipSigmaDependentFvPatchVectorField
(
    const slipSigmaDependentFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    elecM_(ptf.elecM_,false),
    sigma0_(ptf.sigma0_,false),   
    m_(ptf.m_)   
{}

Foam::slipSigmaDependentFvPatchVectorField::
slipSigmaDependentFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    elecM_(Function1<scalar>::New("elecMobility0", dict)),
    sigma0_(Function1<scalar>::New("sigma0", dict)),
    m_(dict.lookupOrDefault("m",0.))
{
    fvPatchField<vector>::operator=
    (
        vectorField("value", dict, p.size())
    );
}
    
Foam::slipSigmaDependentFvPatchVectorField::
slipSigmaDependentFvPatchVectorField
(
    const slipSigmaDependentFvPatchVectorField& tppsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tppsf, iF),
    elecM_(tppsf.elecM_,false),
    sigma0_(tppsf.sigma0_,false),
    m_(tppsf.m_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slipSigmaDependentFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const scalar t = this->db().time().timeOutputValue();
    const scalar elecM = elecM_->value(t);
    const scalar sigma0 = sigma0_->value(t);
    
    const volScalarField& phiE_ =
        db().lookupObject<volScalarField>("phiE");
        
     const volScalarField& sigma_ =
        db().lookupObject<volScalarField>("sigma");
        
    volVectorField Ef(-fvc::grad(phiE_));
     
    vectorField::operator=( elecM * Ef.boundaryField()[patch().index()] * Foam::pow(sigma_.boundaryField()[patch().index()]/sigma0, m_));
         
    fixedValueFvPatchVectorField::updateCoeffs();
}

 
void Foam::slipSigmaDependentFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("m") << m_ << token::END_STATEMENT << nl;
    elecM_->writeData(os);
    sigma0_->writeData(os);
    writeEntry("value", os);
    
}    


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        slipSigmaDependentFvPatchVectorField
    );
}

// ************************************************************************* //
