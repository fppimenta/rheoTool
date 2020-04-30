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

#include "coupledTFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledTFvPatchScalarField::
coupledTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    nbrFieldName_(" "),
    kappaName_(" "),
    isContactResistance_(false),
    hres_(0),
    isCoupledSystem(false)
{}

Foam::coupledTFvPatchScalarField::
coupledTFvPatchScalarField
(
    const coupledTFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    nbrFieldName_(ptf.nbrFieldName_),
    kappaName_(ptf.kappaName_),
    isContactResistance_(ptf.isContactResistance_),
    hres_(ptf.hres_),
    isCoupledSystem(ptf.isCoupledSystem)
{}

Foam::coupledTFvPatchScalarField::
coupledTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    nbrFieldName_(dict.lookupOrDefault<word>("nbrFieldName","T")),
    kappaName_(dict.lookupOrDefault<word>("kappaName","kappa")),
    isContactResistance_(readBool(dict.lookup("isContactResistance"))),
    isCoupledSystem(false)
{
    if (!isA<regionCoupledAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "' not type '" << regionCoupledAMIFvPatch::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
  
    if (isContactResistance_)        
      hres_ = readScalar(dict.lookup("hres"));
    
    fvPatchField<scalar>::operator=
    (
        scalarField("value", dict, p.size())
    );
}
    
Foam::coupledTFvPatchScalarField::
coupledTFvPatchScalarField
(
    const coupledTFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    nbrFieldName_(tppsf.nbrFieldName_),
    kappaName_(tppsf.kappaName_),
    isContactResistance_(tppsf.isContactResistance_),
    hres_(tppsf.hres_),
    isCoupledSystem(tppsf.isCoupledSystem)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::coupledTFvPatchScalarField::deltan() const
{
  const regionCoupledAMIFvPatch& fvp = refCast<const regionCoupledAMIFvPatch>(this->patch());
  
  return (fvp.nf() & (fvp.Cf()-fvp.Cn()));
}



Foam::tmp<Foam::scalarField> Foam::coupledTFvPatchScalarField::nbrDeltan() const
{  
  const regionCoupledAMIFvPatch& fvp = refCast<const regionCoupledAMIFvPatch>(this->patch());
  
  const regionCoupledAMIFvPatch& nbrPatch = fvp.neighbFvPatch();

  tmp<scalarField> tnbrDeltan;
  
  tnbrDeltan = fvp.polyPatch().interpolate(nbrPatch.nf() & (nbrPatch.Cf()-nbrPatch.Cn()) );
                
  return tnbrDeltan;  
}


Foam::tmp<Foam::scalarField> Foam::coupledTFvPatchScalarField::weights() const
{
  const scalarField deltan(this->deltan());
  const scalarField nbrDeltan(this->nbrDeltan());
  
  tmp<scalarField> tw( new scalarField(deltan.size(), 0.));
  scalarField& w = tw.ref();

  forAll(deltan, facei)
  {
    scalar di = deltan[facei];
    scalar dni = nbrDeltan[facei];

    w[facei] = dni/(di + dni);
  }
  
  return tw;
} 



Foam::tmp<Foam::vectorField> Foam::coupledTFvPatchScalarField::delta() const
{
  const regionCoupledAMIFvPatch& fvp = refCast<const regionCoupledAMIFvPatch>(this->patch());
  
  const regionCoupledAMIFvPatch& nbrPatch = fvp.neighbFvPatch();
 
  const vectorField patchD(fvp.Cf()-fvp.Cn());

  tmp<vectorField> tnbrPatchD;
  tnbrPatchD = fvp.polyPatch().interpolate(nbrPatch.Cf()-nbrPatch.Cn());
  
  const vectorField& nbrPatchD = tnbrPatchD();

  tmp<vectorField> tpdv(new vectorField(patchD.size()));
  vectorField& pdv = tpdv.ref();
 
  forAll(patchD, facei)
  {
    const vector& ddi = patchD[facei];
    const vector& dni = nbrPatchD[facei];
    
    // Translational transformation would be here, but the AMI in regionCoupledBase
    // does not allow such transformations. Would need re-writting class regionCoupledBase.
    
    // pdv[facei] = ddi - transform(forwardT()[0], dni);
    pdv[facei] = ddi - dni;
  }
           
  return tpdv;
}



Foam::tmp<Foam::scalarField> Foam::coupledTFvPatchScalarField::deltaCoeffs() const
{
 // Ortho 
 //return  1.0/mag(this->delta());
  
 // Non-ortho 
 tmp<vectorField> deltas(this->delta());  
 return 1.0/max(this->patch().nf() & deltas(), 0.05*mag(deltas()));
}


Foam::word Foam::coupledTFvPatchScalarField::nbrMeshName() const
{
   const regionCoupledBaseFvPatch& rcb = refCast<const regionCoupledBaseFvPatch>(this->patch());    
   const regionCoupledBase& rcbpp = rcb.regionCoupledPatch();
   
   return rcbpp.nbrRegionName();
}

Foam::word Foam::coupledTFvPatchScalarField::nbrFieldName() const
{
   return nbrFieldName_;
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::coupledTFvPatchScalarField::snGrad() const
{   
   const regionCoupledBase& rcbpp = refCast<const regionCoupledBaseFvPatch>(this->patch()).regionCoupledPatch();
   
   word nbrMeshName = rcbpp.nbrRegionName();   
   const fvMesh& thisMesh = this->patch().boundaryMesh().mesh();
   const polyMesh& nbrMesh = thisMesh.time().lookupObject<polyMesh>(nbrMeshName);
   const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[rcbpp.neighbPatchID()];
   scalarField Tc = rcbpp.interpolate(nbrPatch.lookupPatchField<volScalarField, scalar>(nbrFieldName_).patchInternalField());
   
   return this->deltaCoeffs() * (Tc - this->patchInternalField());
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::coupledTFvPatchScalarField::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    // Won't be called because fvpf is non-coupled
    return snGrad();
}


void Foam::coupledTFvPatchScalarField::updateCoeffs()
{
   if (this->updated())
   {
       return;
   }
   
   int oldTag = UPstream::msgType();
   UPstream::msgType() = oldTag+1;
      
   const regionCoupledBaseFvPatch& rcb = refCast<const regionCoupledBaseFvPatch>(this->patch());    
   const regionCoupledBase& rcbpp = rcb.regionCoupledPatch();
   
   //- Lookup conductivity in current mesh
  
   const fvMesh& thisMesh = this->patch().boundaryMesh().mesh();
   
   word thisMeshName = thisMesh.name();
  
   volScalarField& thisKappa = thisMesh.lookupObjectRef<volScalarField>(kappaName_);
   
   fvPatchField<scalar>& kpf = thisKappa.boundaryFieldRef()[this->patch().index()]; 
   
   scalarField kp = kpf.patchInternalField(); 
   
   //- Lookup conductivity in neighb mesh
    
   const polyMesh& nbrMesh = thisMesh.time().lookupObject<polyMesh>(rcbpp.nbrRegionName());
   
   const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[rcbpp.neighbPatchID()];
   
   scalarField kcf = nbrPatch.lookupPatchField<volScalarField, scalar>(kappaName_);
   
   scalarField kc = nbrPatch.lookupPatchField<volScalarField, scalar>(kappaName_).patchInternalField();
   
   kcf = rcbpp.interpolate(kcf);
   kc = rcbpp.interpolate(kc);
   
   //- Lookup temperature and dx in current mesh
   scalarField Tp = this->patchInternalField(); 
   const scalarField dxp(this->deltan());
   
   //- Lookup temperature and dx in neighb mesh
   scalarField Tc = rcbpp.interpolate(nbrPatch.lookupPatchField<volScalarField, scalar>(nbrFieldName_).patchInternalField()); 
   const scalarField dxc(this->nbrDeltan());
   
   scalarField w(this->weights());
   
   // Note: when writting the 2 equations for the flux and temperature at the interface,
   // we can either assume that the face k is equal to the cell k, or that the face k is
   // computed locally according to the interface temperature. The first approach is adopted
   // for implicit coupling, whereas the second one for explicit coupling. Whether using one
   // or the other, we must be consistent in the values of k/kf used to compute Ti, Tci and Tpi.  
   
   if (!isContactResistance_)
   {  
     scalarField Ti(Tp.size(), 0.);
     
     if (isCoupledSystem)
     {
       // Set interface conductivity
       kpf = 1./( w/kc + (1.0 - w)/kp );
       
       // Set interface temperature 
       Ti = ( kc*Tc*dxp + kp*Tp*dxc ) / ( kc*dxp + kp*dxc );
     }
     else
     {
       // The conductivity is calculated in the thermo-model as a function of the temperature.
       // If thermal conductivity is independent of temperature, then this is equivalent 
       // to applying a zero-gradient condition to conductivity.
       
       // Set interface temperature 
       Ti = ( kcf*Tc*dxp + kpf*Tp*dxc ) / ( kcf*dxp + kpf*dxc );
     }
    
     this->operator==(Ti);  
   }
   else 
   {
     scalar h(hres_);
 
     scalarField Tci(Tp.size(), 0.);
     scalarField Tpi(Tp.size(), 0.);
     
     if (isCoupledSystem)
     { 
       // Set interface conductivity
       kpf = 1./( w/kc + (1.0 - w)/kp + this->deltaCoeffs()/h);
       
       // Set interface temperature    
       Tci = ( h*(kc*Tc*dxp+kp*Tp*dxc) + kc*Tc*kp ) / ( h*(kc*dxp+kp*dxc) + kc*kp );     
       Tpi = ( h*dxp*Tci + kp*Tp )/( h*dxp + kp );  
     }
     else
     {
       // The conductivity is calculated in the thermo-model as a function of the temperature.
       // If thermal conductivity is independent of temperature, then this is equivalent 
       // to applying a zero-gradient condition to conductivity.
       
       // Set interface temperature    
       Tci = ( h*(kcf*Tc*dxp+kpf*Tp*dxc) + kcf*Tc*kpf ) / ( h*(kcf*dxp+kpf*dxc) + kcf*kpf );     
       Tpi = ( h*dxp*Tci + kpf*Tp )/( h*dxp + kpf ); 
     }
 
     this->operator==(Tpi);         
   }
   
   fixedValueFvPatchField<scalar>::updateCoeffs();
  
   UPstream::msgType() = oldTag;
}


// Same value and gradient coeffs as a regular cylicAMIFvPatchField 

// Note: only gradientCoeffs are supposed to be used by this BC, in the 
// Laplacian operator. valueCoeffs are used in the div() operator but are
// supposed to be associated to a zero flux through the face (no contribution),
// as expected for a solid-fluid interface.
// If convection is to be used in a coupled way, then this BC must be of coupled
// type to work correctly with gaussDefCmpwConvectionScheme schemes (not sure about
// gaussConvectionScheme; see the source files).

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::coupledTFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{    
    scalarField thisW(this->weights());
    return scalar(pTraits<scalar>::one)*thisW;
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::coupledTFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    scalarField thisW(this->weights());
    return scalar(pTraits<scalar>::one)*(1.0 - thisW);
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::coupledTFvPatchScalarField::gradientInternalCoeffs() const
{
  if (isCoupledSystem)
  {
    return -pTraits<scalar>::one*this->deltaCoeffs();
  }
  else
  {
    return -pTraits<scalar>::one*this->patch().deltaCoeffs();
  }  
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::coupledTFvPatchScalarField::gradientBoundaryCoeffs() const
{
  if (isCoupledSystem)
  {
    return -this->gradientInternalCoeffs();
  }
  else
  {
    return this->patch().deltaCoeffs()*(*this);
  }
}

void Foam::coupledTFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("nbrFieldName") << nbrFieldName_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappaName") << kappaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("isContactResistance") << isContactResistance_ << token::END_STATEMENT << nl;
    os.writeKeyword("hres") << hres_ << token::END_STATEMENT << nl;
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coupledTFvPatchScalarField
    );
}

// ************************************************************************* //
