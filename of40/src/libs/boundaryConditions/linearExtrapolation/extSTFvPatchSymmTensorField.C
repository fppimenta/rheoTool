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

#include "extSTFvPatchSymmTensorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extSTFvPatchSymmTensorField::
extSTFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(p, iF)
{}


Foam::extSTFvPatchSymmTensorField::
extSTFvPatchSymmTensorField
(
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchSymmTensorField(p, iF)
{
    fvPatchSymmTensorField::operator=(symmTensorField("value", dict, p.size()));
}


Foam::extSTFvPatchSymmTensorField::
extSTFvPatchSymmTensorField
(
    const extSTFvPatchSymmTensorField& ptf,
    const fvPatch& p,
    const DimensionedField<symmTensor, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchSymmTensorField(ptf, p, iF, mapper)
{}


Foam::extSTFvPatchSymmTensorField::
extSTFvPatchSymmTensorField
(
    const extSTFvPatchSymmTensorField& pivpvf,
    const DimensionedField<symmTensor, volMesh>& iF
)
:
    fixedValueFvPatchSymmTensorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void Foam::extSTFvPatchSymmTensorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
 
     word sTensorFName( this->internalField().name() );

     const fvMesh& mesh = patch().boundaryMesh().mesh();  
 
     const volSymmTensorField& taufield = db().objectRegistry::lookupObject<volSymmTensorField>(sTensorFName);  
  
     symmTensorField tau_av = symmTensorField( patch().size(), pTraits<symmTensor>::zero ) ;  

     for (int i=0; i<6 ;i++) 
     {

        volVectorField gradT = fvc::grad(taufield.component(i), "extSTGrad");

        forAll(patch(), facei )        
          {
 
             vector r_face = patch().Cf()[facei];  

             label cell1 = patch().faceCells()[facei];

             vector r_cell1 = mesh.C()[ cell1 ];  

             vector fc1 = (r_face - r_cell1);

             tau_av[facei].component(i) = taufield[ cell1 ].component(i) + (gradT[ cell1 ] & fc1 );
           } 
 
     }
     
     symmTensorField::operator=(tau_av);
 
     fixedValueFvPatchSymmTensorField::updateCoeffs();
    

}
void Foam::extSTFvPatchSymmTensorField::write(Ostream& os) const
{
    fvPatchSymmTensorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchSymmTensorField,
        extSTFvPatchSymmTensorField
    );
}


// ************************************************************************* //
