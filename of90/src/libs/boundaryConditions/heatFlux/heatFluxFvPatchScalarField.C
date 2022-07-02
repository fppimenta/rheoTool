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

#include "heatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatFluxFvPatchScalarField::heatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1),
    Q_(0),
    Ta_(),
    Ts_(),
    emissivity_(0),
    kappaName_(" ")
{}


Foam::heatFluxFvPatchScalarField::heatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1),
    Q_(0),
    Ta_(),
    Ts_(),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0)),
    kappaName_(dict.lookupOrDefault<word>("kappaName","kappa"))
{
    dict.lookup("Q") >> Q_;
    q_ = scalarField("q", dict, p.size());
    h_ = scalarField("h", dict, p.size());
    Ta_ = Function1<scalar>::New("Ta", dict);
    Ts_ = Function1<scalar>::New("Ts", dict);

    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::heatFluxFvPatchScalarField::heatFluxFvPatchScalarField
(
    const heatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(p, iF),
    curTimeIndex_(-1),
    Q_(ptf.Q_),
    Ta_(ptf.Ta_, false),
    Ts_(ptf.Ts_, false),
    emissivity_(ptf.emissivity_),
    kappaName_(ptf.kappaName_)
{
    // Map gradient. Set unmapped values and overwrite with mapped ptf
    gradient() = 0.0;
    mapper(gradient(), ptf.gradient());
    
   // q_.setSize(mapper.size());
   // q_.map(ptf.q_, mapper);
    mapper(q_, ptf.q_);
    
   // h_.setSize(mapper.size());
   // h_.map(ptf.h_, mapper);
    mapper(h_, ptf.h_);
    
    // Evaluate the value field from the gradient if the internal field is valid
    if (notNull(iF) && iF.size())
    {
        scalarField::operator=
        (
            // patchInternalField() + gradient()/patch().deltaCoeffs()
            // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
            // which fails for AMI patches for some mapping operations
            patchInternalField() + gradient()*(patch().nf() & patch().delta())
        );
    }
    else
    {
        // Enforce mapping of values so we have a valid starting value. This
        // constructor is used when reconstructing fields
        mapper(*this, ptf);
    }
}


Foam::heatFluxFvPatchScalarField::heatFluxFvPatchScalarField
(
    const heatFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    curTimeIndex_(-1),
    Q_(wbppsf.Q_),
    q_(wbppsf.q_),
    h_(wbppsf.h_),
    Ta_(wbppsf.Ta_, false),
    Ts_(wbppsf.Ts_, false),
    emissivity_(wbppsf.emissivity_),
    kappaName_(wbppsf.kappaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
void Foam::heatFluxFvPatchScalarField::updateCoeffs()
{
  if (updated())
  {
    return;
  }

  scalarField Tc(this->patchInternalField());
  scalarField delta(1./patch().deltaCoeffs());
  const scalar Ta = Ta_->value(this->db().time().timeOutputValue());
  const scalar Ts = Ts_->value(this->db().time().timeOutputValue());
  
  const fvPatchScalarField& k = this->patch().lookupPatchField<volScalarField, scalar>(kappaName_);
  
  alpha = -h_ - emissivity_*sigma.value()*pow3(*this);
  beta = h_*Ta + emissivity_*sigma.value()*pow4(Ts) + q_ + Q_/gSum(patch().magSf());
  
  curTimeIndex_ = this->db().time().timeIndex();

  gradient() = (alpha*Tc+beta)/(k-delta*alpha);
  
  fixedGradientFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::heatFluxFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
   const fvPatchScalarField& k = this->patch().lookupPatchField<volScalarField, scalar>(kappaName_);
   
   return ( 1./(1.-alpha/(k*patch().deltaCoeffs())) );
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::heatFluxFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
   const fvPatchScalarField& k = this->patch().lookupPatchField<volScalarField, scalar>(kappaName_);
   
   return (beta/(k*patch().deltaCoeffs()-alpha));
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::heatFluxFvPatchScalarField::gradientInternalCoeffs() const
{
   const fvPatchScalarField& k = this->patch().lookupPatchField<volScalarField, scalar>(kappaName_);
   
   return (alpha/(k-alpha/patch().deltaCoeffs()));
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::heatFluxFvPatchScalarField::gradientBoundaryCoeffs() const
{
   const fvPatchScalarField& k = this->patch().lookupPatchField<volScalarField, scalar>(kappaName_);
   
   return (beta/(k-alpha/patch().deltaCoeffs()));
}

void Foam::heatFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    
    os.writeKeyword("Q") << Q_ << token::END_STATEMENT << nl;
    writeEntry(os, "q", q_);
    writeEntry(os, "h", h_);
    writeEntry(os, Ta_());
    writeEntry(os, Ts_());
    os.writeKeyword("emissivity") << emissivity_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappaName") << kappaName_ << token::END_STATEMENT << nl;
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        heatFluxFvPatchScalarField
    );
}


// ************************************************************************* //
