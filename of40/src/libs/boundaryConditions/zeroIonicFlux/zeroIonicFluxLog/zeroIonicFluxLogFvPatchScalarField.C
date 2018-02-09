/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "zeroIonicFluxLogFvPatchScalarField.H"
#include "EDFEquation.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroIonicFluxLogFvPatchScalarField::zeroIonicFluxLogFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    zib_()
{}


Foam::zeroIonicFluxLogFvPatchScalarField::zeroIonicFluxLogFvPatchScalarField
(
    const zeroIonicFluxLogFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    zib_(ptf.zib_)
{}


Foam::zeroIonicFluxLogFvPatchScalarField::zeroIonicFluxLogFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    zib_(0)
{
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
    
    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    PtrList<entry> specEntries_(elecDict.subDict("parameters").lookup("species"));
    
    string fieldN(this->internalField().name());
    
    fieldN.erase(fieldN.end()-3, fieldN.end()); // Remove "Log" suffix
    
    forAll (specEntries_, specI)
    {    
       if ( specEntries_[specI].keyword() == fieldN )
        { 
          dimensionedScalar zid_(specEntries_[specI].dict().lookup("z"));
          zib_ = zid_.value();
          break;
        }
    } 
}


Foam::zeroIonicFluxLogFvPatchScalarField::zeroIonicFluxLogFvPatchScalarField
(
    const zeroIonicFluxLogFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    zib_(wbppsf.zib_)
{}


Foam::zeroIonicFluxLogFvPatchScalarField::zeroIonicFluxLogFvPatchScalarField
(
    const zeroIonicFluxLogFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    zib_(wbppsf.zib_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroIonicFluxLogFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
        
    }
    
    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    dimensionedScalar T_(elecDict.subDict("parameters").lookup("T"));
    
    const volScalarField& psi_ =
        db().lookupObject<volScalarField>("psi");

    const volScalarField& phiE_ =
        db().lookupObject<volScalarField>("phiE");
        
    scalarField Epatch( (psi_+phiE_)().boundaryField()[patch().index()].snGrad() );
   
    scalar eK_(Foam::EDFEquation::eK_.value());
    scalar kbK_(Foam::EDFEquation::kbK_.value());
          
    gradient() = - eK_ * zib_ * Epatch / (kbK_ * T_.value());      

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::zeroIonicFluxLogFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        zeroIonicFluxLogFvPatchScalarField
    );
}

// ************************************************************************* //
