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

#include "poiseuilleVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

 
template<>
const char* Foam::NamedEnum
<
  Foam::poiseuilleVelocityFvPatchVectorField::geomOpt,
  3
>::names[] =
{
  "Duct3D",
  "Duct2D",
  "Circular"
};
  
const Foam::NamedEnum<Foam::poiseuilleVelocityFvPatchVectorField::geomOpt, 3> Foam::poiseuilleVelocityFvPatchVectorField::geomOptionNames_;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poiseuilleVelocityFvPatchVectorField::
poiseuilleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}

Foam::poiseuilleVelocityFvPatchVectorField::
poiseuilleVelocityFvPatchVectorField
(
    const poiseuilleVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    center_(ptf.center_),  
    dirA_(ptf.dirA_),  
    dirB_(ptf.dirB_), 
    wA_(ptf.wA_), 
    wB_(ptf.wB_), 
    imax_(ptf.imax_), 
    uMean_(ptf.uMean_),
    R_(ptf.R_), 
    U_(ptf.U_),     
    geomOpts_(ptf.geomOpts_)
{}
 
Foam::poiseuilleVelocityFvPatchVectorField::
poiseuilleVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    U_(p.size(), vector::zero),
    geomOpts_(geomOptionNames_.read(dict.lookup("geometryType")))  
{
    switch (geomOpts_) 
    {     
     case goDuct3D : 
      center_ = dict.lookup("center");
      dirA_ = dict.lookup("dirA"); dirA_ /= mag(dirA_);
      dirB_ = dict.lookup("dirB"); dirB_ /= mag(dirB_);
      wB_ = readScalar(dict.lookup("wB"));
      imax_ = dict.lookupOrDefault<label>("maxIter", 1000);
      wA_ = readScalar(dict.lookup("wA"));
      uMean_ = readScalar(dict.lookup("uMean"));
      computeUDuct3D();
      break;
     
     case goDuct2D :
      center_ = dict.lookup("center");
      dirA_ = dict.lookup("dirA"); dirA_ /= mag(dirA_);
      wA_ = readScalar(dict.lookup("wA"));
      uMean_ = readScalar(dict.lookup("uMean"));
      computeUDuct2D();
      break;
      
     case goCircular :
      center_ = dict.lookup("center");
      R_ = readScalar(dict.lookup("R"));
      uMean_ = readScalar(dict.lookup("uMean"));
      computeUCircular();
      break;    
    }
     
    fvPatchField<vector>::operator=
    (
        U_
    );   
}
    
Foam::poiseuilleVelocityFvPatchVectorField::
poiseuilleVelocityFvPatchVectorField
(
    const poiseuilleVelocityFvPatchVectorField& tppsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tppsf, iF),
    center_(tppsf.center_),  
    dirA_(tppsf.dirA_),  
    dirB_(tppsf.dirB_), 
    wA_(tppsf.wA_), 
    wB_(tppsf.wB_), 
    imax_(tppsf.imax_), 
    uMean_(tppsf.uMean_),
    R_(tppsf.R_), 
    U_(tppsf.U_),     
    geomOpts_(tppsf.geomOpts_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::poiseuilleVelocityFvPatchVectorField::computeUDuct3D()
{
  // Compute dp/dx based on Q
  scalar cumsum(0.);
  for (int i=1; i<imax_; i+=2)
  {
    cumsum += Foam::tanh(i*M_PI*wB_/(2.*wA_))/(i*i*i*i*i);
  }
 
  scalar Q = 4.*wA_*wB_*uMean_;
  
  scalar gradp = Q/((4./3.)*wB_*wA_*wA_*wA_*(1.-(192.*wA_*cumsum)/(wB_*M_PI*M_PI*M_PI*M_PI*M_PI)));
  
  const fvPatch& p = patch();
  const polyPatch& pp = p.patch();
  const vectorField& C(p.Cf());
  vectorField S(pp.faceAreas()); 
  S /= mag(S);
  label sz(S.size());
  reduce(sz, sumOp<label>()); 
  vector n = -gSum(S)/sz;
  n /= mag(n);
   
  forAll(C, i)
  {
    scalar y = (dirA_ & (C[i] - center_));
    scalar z = (dirB_ & (C[i] - center_)); 
    scalar cumsum(0.);
    for (int j=1; j<imax_; j+=2)
    {
      scalar rat;      
      if (0.5*j*M_PI*wB_/wA_ > 700)
      {
        rat = 0;        
      }
      else
      {
        rat = Foam::cosh(0.5*j*M_PI*z/wA_)/Foam::cosh(0.5*j*M_PI*wB_/wA_);        
      }   
      cumsum += Foam::pow(-1, (j-1)/2) * (1. - rat) * Foam::cos(0.5*j*M_PI*y/wA_)/(j*j*j);           
    }  
    
    U_[i] = n*(16.*wA_*wA_*cumsum*gradp/(M_PI*M_PI*M_PI));
  }
  Info << "Poiseuille BC at " << pp.name() << endl;
  Info << "Average U: " << gSum(U_*mag(pp.faceAreas()))/gSum(mag(pp.faceAreas())) << tab 
       << "Max |U|: "   << gMax(mag(U_)) << nl << endl; 
}

void Foam::poiseuilleVelocityFvPatchVectorField::computeUDuct2D()
{
  const fvPatch& p = patch();
  const polyPatch& pp = p.patch();
  const vectorField& C(p.Cf());
  vectorField S(pp.faceAreas()); 
  S /= mag(S);
  label sz(S.size());
  reduce(sz, sumOp<label>()); 
  vector n = -gSum(S)/sz;
  n /= mag(n);
  
  forAll(C, i)
  {
     scalar y = (dirA_ & (C[i] - center_));
     U_[i] = n * 1.5 * uMean_ * (1. - (y/wA_)*(y/wA_)); 
  } 
  Info << "Poiseuille BC at " << pp.name() << endl;
  Info << "Average U: " << gSum(U_*mag(pp.faceAreas()))/gSum(mag(pp.faceAreas())) << tab 
       << "Max |U|: "   << gMax(mag(U_)) << nl << endl; 
}

void Foam::poiseuilleVelocityFvPatchVectorField::computeUCircular()
{
  const fvPatch& p = patch();
  const polyPatch& pp = p.patch();
  const vectorField& C(p.Cf());
  vectorField S(pp.faceAreas()); 
  S /= mag(S);
  label sz(S.size());
  reduce(sz, sumOp<label>()); 
  vector n = -gSum(S)/sz;
  n /= mag(n);
  
  forAll(C, i)
  {
     scalar r(mag(C[i]-center_));
     U_[i] = n * 2.*uMean_ * (1. - (r/R_)*(r/R_)); 
  } 
  Info << "Poiseuille BC at " << pp.name() << endl;
  Info << "Average U: " << gSum(U_*mag(pp.faceAreas()))/gSum(mag(pp.faceAreas())) << tab 
       << "Max |U|: "   << gMax(mag(U_)) << nl << endl; 
}

void Foam::poiseuilleVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    } 
   
    vectorField::operator=( U_ );
         
    fixedValueFvPatchVectorField::updateCoeffs();
}

 
void Foam::poiseuilleVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("center") << center_ << token::END_STATEMENT << nl;
    os.writeKeyword("uMean") << uMean_ << token::END_STATEMENT << nl;
    os.writeKeyword("geometryType") << geomOptionNames_[geomOpts_] << token::END_STATEMENT << nl;
    
    
    if (geomOpts_ == goDuct3D)
    {
      os.writeKeyword("dirA") << dirA_ << token::END_STATEMENT << nl;
      os.writeKeyword("dirB") << dirB_ << token::END_STATEMENT << nl;
      os.writeKeyword("wA") << wA_ << token::END_STATEMENT << nl;
      os.writeKeyword("wB") << wB_ << token::END_STATEMENT << nl;
      os.writeKeyword("imax") << imax_ << token::END_STATEMENT << nl;
    }
    else if (geomOpts_ == goDuct2D)
    {
      os.writeKeyword("dirA") << dirA_ << token::END_STATEMENT << nl;
      os.writeKeyword("wA") << wA_ << token::END_STATEMENT << nl;
    }
    else if (geomOpts_ == goCircular)
    {
      os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;
    }
    
    writeEntry(os,"value",*this);
}
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        poiseuilleVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
