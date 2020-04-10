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

#include "boltzmannEquilibriumFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "EDFEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boltzmannEquilibriumFvPatchScalarField::
boltzmannEquilibriumFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    zib_(1.),
    c0_(),
    psi0_()
{}

Foam::boltzmannEquilibriumFvPatchScalarField::
boltzmannEquilibriumFvPatchScalarField
(
    const boltzmannEquilibriumFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    zib_(ptf.zib_),
    c0_(ptf.c0_, false),
    psi0_(ptf.psi0_, false)
{}

Foam::boltzmannEquilibriumFvPatchScalarField::
boltzmannEquilibriumFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    zib_(1),
    c0_(Function1<scalar>::New("c0", dict)),
    psi0_(Function1<scalar>::New("psi0", dict))
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
    
    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    PtrList<entry> specEntries_(elecDict.subDict("parameters").lookup("species"));
    
    forAll (specEntries_, specI)
    {    
       if ( specEntries_[specI].keyword() == this->internalField().name() )
        { 
          dimensionedScalar zid_(specEntries_[specI].dict().lookup("z"));
          zib_ = zid_.value();
          break;
        }
    } 
}
    
Foam::boltzmannEquilibriumFvPatchScalarField::
boltzmannEquilibriumFvPatchScalarField
(
    const boltzmannEquilibriumFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    zib_(tppsf.zib_),
    c0_(tppsf.c0_, false),
    psi0_(tppsf.psi0_,false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boltzmannEquilibriumFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const scalar t = this->db().time().timeOutputValue();
    const scalar c0 = c0_->value(t);
    const scalar psi0 = psi0_->value(t);

    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    dimensionedScalar T_(elecDict.subDict("parameters").lookup("T"));
    
    const fvPatchField<scalar>& psib_ = patch().lookupPatchField<volScalarField, scalar>("psi");
    
    scalar eK_(Foam::EDFEquation::eK_.value());
    scalar kbK_(Foam::EDFEquation::kbK_.value());
    
    scalar alpha = eK_ * zib_ / ( kbK_ * T_.value() );
     
    scalarField::operator=( Foam::max(0., c0 * Foam::exp(-alpha*(psib_ - psi0)) ) );
     
    fixedValueFvPatchScalarField::updateCoeffs();
}

 
void Foam::boltzmannEquilibriumFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    c0_->writeData(os);
    psi0_->writeData(os);
    //writeEntry("value", os);
    writeEntry(os,"value");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        boltzmannEquilibriumFvPatchScalarField
    );
}

// ************************************************************************* //
