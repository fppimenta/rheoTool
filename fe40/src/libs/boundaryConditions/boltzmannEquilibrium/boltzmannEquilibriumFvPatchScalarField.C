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
    zib_(),
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
    c0_(ptf.c0_),
    psi0_(ptf.psi0_)
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
    zib_(0),
    c0_(readScalar(dict.lookup("c0"))),
    psi0_(readScalar(dict.lookup("psi0")))
{
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
    
    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    PtrList<entry> specEntries_(elecDict.subDict("parameters").lookup("species"));
    
    forAll (specEntries_, specI)
    {    
       if ( specEntries_[specI].keyword() == this->dimensionedInternalField().name() )
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
    c0_(tppsf.c0_),
    psi0_(tppsf.psi0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boltzmannEquilibriumFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
     

    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    dimensionedScalar T_(elecDict.subDict("parameters").lookup("T"));
    
    const fvPatchField<scalar>& psib_ = patch().lookupPatchField<volScalarField, scalar>("psi");
    
    scalar eK_(Foam::EDFEquation::eK_.value());
    scalar kbK_(Foam::EDFEquation::kbK_.value());
    
    scalar alpha = eK_ * zib_ / ( kbK_ * T_.value() );
     
    scalarField::operator=( Foam::max(0., c0_ * Foam::exp(-alpha*(psib_ - psi0_)) ) );
     
    fixedValueFvPatchScalarField::updateCoeffs();
}

 
void Foam::boltzmannEquilibriumFvPatchScalarField::write(Ostream& os) const
{

    fvPatchScalarField::write(os);
    os.writeKeyword("c0")
        << c0_ << token::END_STATEMENT << nl;
    os.writeKeyword("psi0")
        << psi0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
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
