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

#include "uShaftFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "constitutiveModel.H"
#include "constitutiveTwoPhaseMixture.H"
 
// * * * * * * * * * * * * * * * * Static Data  * * * * * * * * * * * * * * //
 
Foam::autoPtr<Foam::volSymmTensorField>
Foam::uShaftFvPatchVectorField::tauTotalPtr(NULL);

void Foam::uShaftFvPatchVectorField::
updateTauTotalPtr(const Foam::volSymmTensorField& tau)
{
  tauTotalPtr.reset(new volSymmTensorField(tau));
  tauTotalPtr().rename("tauTotalPtr");
} 
         
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uShaftFvPatchVectorField::
uShaftFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    model_("nonLinearNavierSlip"),
    m_(0.),
    knl_(0.),
    alpha_(0.),
    beta_(0.),
    etaS_(0.),
    URF_(1.),
    isTwoPhaseFlow_(false),
    isMovingWall_(false),
    offset_(vector::zero),
    amplitude_(vector::zero),
    frequency_(0.),
    delay_(0)
{}

Foam::uShaftFvPatchVectorField::
uShaftFvPatchVectorField
(
    const uShaftFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    model_(ptf.model_),
    m_(ptf.m_),   
    knl_(ptf.knl_),
    alpha_(ptf.alpha_),   
    beta_(ptf.beta_),
    etaS_(ptf.etaS_),
    URF_(ptf.URF_),
    isTwoPhaseFlow_(ptf.isTwoPhaseFlow_),
    isMovingWall_(ptf.isMovingWall_),
    offset_(ptf.offset_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    delay_(ptf.delay_)   
{}

Foam::uShaftFvPatchVectorField::
uShaftFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    model_(dict.lookup("model")),
    m_(model_ == "nonLinearNavierSlip" ? readScalar(dict.lookup("m")) : 0.),
    knl_(model_ == "nonLinearNavierSlip" ? readScalar(dict.lookup("knl")) : 0.),
    alpha_(model_ == "slipTT" ? readScalar(dict.lookup("alpha")) : 0.),
    beta_(model_ == "slipTT" ? readScalar(dict.lookup("beta")) : 0.),
    etaS_(0.),
    URF_(readScalar(dict.lookup("URF"))),
    isTwoPhaseFlow_(readBool(dict.lookup("isTwoPhaseFlow"))),
    isMovingWall_(readBool(dict.lookup("isMovingWall"))),
    offset_(dict.lookup("offset")),
    amplitude_(dict.lookup("amplitude")),
    frequency_(readScalar(dict.lookup("frequency"))),
    delay_(readScalar(dict.lookup("delay")))
{
    if (!(model_ == "nonLinearNavierSlip" ||
          model_ == "slipTT")
       )
    {
        FatalErrorIn("Foam::uShaftFvPatchVectorField::uShaftFvPatchVectorField\n")
        << "\nUnknown model <" << model_ <<"> specified in the Navier slip BC for field "
        << this->dimensionedInternalField().name() << ".\n"
        << "\nAvailable models are:\n"
        << "\n. nonLinearNavierSlip" <<"\n. slipTT"  
        << abort(FatalError);
    }
    
    if (m_<0. || beta_<0. || alpha_<0. )
     {
        FatalErrorIn("Foam::uShaftFvPatchVectorField::uShaftFvPatchVectorField\n")
        << "\nError defining negative value of m and/or alpha and/or beta in the Navier slip BC for field "
        << this->dimensionedInternalField().name() << ".\n"
        << "\nPlease ensure that m>0, alpha>0 and beta>0.\n"
        << abort(FatalError);
     }
    
    fvPatchField<vector>::operator=
    (
        vectorField("value", dict, p.size())
    );
}
    
Foam::uShaftFvPatchVectorField::
uShaftFvPatchVectorField
(
    const uShaftFvPatchVectorField& tppsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tppsf, iF),
    model_(tppsf.model_),
    m_(tppsf.m_),
    knl_(tppsf.knl_),
    alpha_(tppsf.alpha_),
    beta_(tppsf.beta_),
    etaS_(tppsf.etaS_),
    URF_(tppsf.URF_),
    isTwoPhaseFlow_(tppsf.isTwoPhaseFlow_),
    isMovingWall_(tppsf.isMovingWall_),
    offset_(tppsf.offset_),
    amplitude_(tppsf.amplitude_),
    frequency_(tppsf.frequency_),
    delay_(tppsf.delay_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::uShaftFvPatchVectorField::updateTau()
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();   
  
    // Update the total stress
    if (isTwoPhaseFlow_) 
    {     
      constitutiveTwoPhaseMixture& constEq_ = const_cast<constitutiveTwoPhaseMixture&>
      (
        db().lookupObject<constitutiveTwoPhaseMixture>("constitutiveProperties")
      );
               
      if (mesh.time().timeIndex() == mesh.time().startTimeIndex()+1) 
      {
        uShaftFvPatchVectorField::updateTauTotalPtr(constEq_.tauTotalMF()());
        return;
      }
      
      // Only update if needed 
      if (!uShaftFvPatchVectorField::tauTotalPtr().upToDate("U"))
      {      
        uShaftFvPatchVectorField::updateTauTotalPtr(constEq_.tauTotalMF()());
      }   
    } 
    else
    {      
      constitutiveModel& constEq_ = const_cast<constitutiveModel&>
      (
        db().lookupObject<constitutiveModel>("constitutiveProperties")
      );
                         
      if (mesh.time().timeIndex() == mesh.time().startTimeIndex()+1) 
      {
        uShaftFvPatchVectorField::updateTauTotalPtr(constEq_.tauTotal()());
        return;
      }
      
      // Only update if needed      
      if (!uShaftFvPatchVectorField::tauTotalPtr().upToDate("U"))
      { 
        uShaftFvPatchVectorField::updateTauTotalPtr(constEq_.tauTotal()());
      }     
    }
}

void Foam::uShaftFvPatchVectorField::calcUws
(
  vectorField& uws,
  vectorField& n
)
{

  label patchi = this->patch().index();
    
  // Total extra-stress
  symmTensorField tauP = tauTotalPtr().boundaryField()[patchi]; 
    
  if (model_ == "nonLinearNavierSlip")
   {
     if (m_ == 1)
     {
        uws = -knl_ * ( (tauP&n) - n*((tauP&n)&n) );
     }
     else
     {
        vectorField tt = (tauP&n) - n*((tauP&n)&n);
        uws = -knl_ * Foam::pow(Foam::mag(tt), m_) * tt/(mag(tt) + 1e-20);
     }
   }
  else if (model_ == "slipTT")
   {
      vectorField tt = (tauP&n) - n*((tauP&n)&n);
      uws = -alpha_ * tt / Foam::sqrt(1.-mag(tt)*beta_);
   } 

}

void Foam::uShaftFvPatchVectorField::updateCoeffs()
{

    if (updated())
    {
        return;
    }
          
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();  
 
    // Since we can not initEvaluate due to possible undefinition of tauTotal(),
    // we need to force initialization here 
    if (mesh.time().timeIndex() == mesh.time().startTimeIndex()+1)  
    {      
       updateTau();
    }
 
    vectorField uws = vectorField( this->patch().size(), pTraits<vector>::zero );
    vectorField n = this->patch().nf();
    
    calcUws(uws, n);
  
    // Note: under-relaxation is only applied to static meshes
    if (mesh.moving() && isMovingWall_)
    {
        const scalar& t = this->db().time().timeOutputValue();
        vectorField::operator=(uws + offset_ + amplitude_*Foam::cos(2.*M_PI*frequency_*t + delay_));
    }
    else
    {
        if (URF_>0 && URF_<1)
         {
           vectorField::operator=( *this * (1.-URF_) + URF_*uws );
         }
        else
         {
           vectorField::operator=( uws );
         }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::uShaftFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{     
    updateTau();
    
    if (!updated())
    {
        updateCoeffs();
    }
 
    fixedValueFvPatchVectorField::evaluate();  
}

 
void Foam::uShaftFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("model") << model_ << token::END_STATEMENT << nl;
    
    if (model_ == "nonLinearNavierSlip")
     {
       os.writeKeyword("m") << m_ << token::END_STATEMENT << nl;
       os.writeKeyword("knl") << knl_ << token::END_STATEMENT << nl;
     }
    else if (model_ == "slipTT")
     {
       os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
       os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;
     }
   
    os.writeKeyword("URF") << URF_ << token::END_STATEMENT << nl;
    os.writeKeyword("isTwoPhaseFlow") << isTwoPhaseFlow_ << token::END_STATEMENT << nl;
    os.writeKeyword("isMovingWall") << isMovingWall_ << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude") <<amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("delay") << delay_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        uShaftFvPatchVectorField
    );
}

// ************************************************************************* //
