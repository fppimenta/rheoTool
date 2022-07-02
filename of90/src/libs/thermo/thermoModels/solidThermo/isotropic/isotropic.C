/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "isotropic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace solidThermoModels
 {
   defineTypeNameAndDebug(isotropic, 0);
   addToRunTimeSelectionTable(solidThermoModel, isotropic, dictionary);
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermoModels::isotropic::isotropic
(
  const word& name,
  const fvMesh& mesh
)
:
  solidThermoModel(name, mesh),
  kappa_
  (
    IOobject
    (
      "kappa",
      mesh.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("kt",dimensionSet(1, 1, -3, -1, 0, 0, 0),0.),
    extrapolatedCalculatedFvPatchField<scalar>::typeName
  ),
  rhoCp_(dict_.lookup("rhoCp")),
  k0_(dimensionedScalar("k0",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.)),
  k1_(dimensionedScalar("k1",dimensionSet(1, 1, -3, -2, 0, 0, 0), 0.)),
  k2_(dimensionedScalar("k2",dimensionSet(1, 1, -3, -3, 0, 0, 0), 0.)),  
  k3_(dimensionedScalar("k3",dimensionSet(1, 1, -3, -4, 0, 0, 0), 0.)),
  radiation_(radiationModel::New(T_))
{
  // Fill-in kappaCoeffs
  scalar tmpV;
  dict_.subDict("kappaCoeffs").lookup("a0") >> tmpV; k0_.value() = tmpV;  
  dict_.subDict("kappaCoeffs").lookup("a1") >> tmpV; k1_.value() = tmpV;
  dict_.subDict("kappaCoeffs").lookup("a2") >> tmpV; k2_.value() = tmpV;
  dict_.subDict("kappaCoeffs").lookup("a3") >> tmpV; k3_.value() = tmpV;
  
  // Initialize kappa from T
  kappa_ = k0_ 
         + k1_*T_
         + k2_*T_*T_
         + k3_*T_*T_*T_;
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
 
void Foam::solidThermoModels::isotropic::DiffNo
(
   const fvMesh& mesh,
   const Time& runTime,
   scalar& maxDiffNo,
   scalar& meanDiffNo
) const
{

 surfaceScalarField kapparhoCpbyDelta
 (
    sqr(mesh.surfaceInterpolation::deltaCoeffs())
    *fvc::interpolate(kappa_)/rhoCp_
 );

 maxDiffNo = max(kapparhoCpbyDelta).value()*runTime.deltaTValue();
 meanDiffNo = average(kapparhoCpbyDelta).value()*runTime.deltaTValue();
 
 Info<< "Region: " << mesh.name() << " Diffusion Number mean: " << meanDiffNo
        << " max: " << maxDiffNo << endl;
}
 

void Foam::solidThermoModels::isotropic::correct
(
  autoPtr<coupledSolver>& cpsT,
  const fvModels& fvModel,
  fvConstraints& fvConstraint,
  int nNonOrtoC
)
{

// Only performs non-ortho correctors loop if not coupled
int nNonOC(cpsT.empty()?max(1,nNonOrtoC+1):1);

for (int jj = 0; jj<nNonOC; jj++)
{
  // Update kappa. kappa(T) = k0 + k1*T + k2*T^2 + k3*T^3
  if (k1_.value() != 0 || k2_.value() != 0 || k3_.value() != 0 )
  { 
    kappa_ = k0_ 
           + k1_*T_
           + k2_*T_*T_
           + k3_*T_*T_*T_;
  }   
  
  // TEqn
  fvScalarMatrix TEqn
  (
      fvm::ddt(rhoCp_, T_)
    - fvm::laplacian(kappa_, T_)
    ==
      rhoCp_*radiation_->ST(rhoCp_, T_) 
  );
  
  if (fvModel.addsSupToField(T_.name()))
  { 
    // Only the dimensions of rhoCpField are used in fvOptions. The value
    // itself is not used and can be any.
    volScalarField rhoCpField
    (
     IOobject
     (
      "tmpOpt",
      mesh().time().constant(),
      mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
     ),
     mesh(),
     rhoCp_,
     extrapolatedCalculatedFvPatchField<scalar>::typeName
    );
     
    TEqn -= fvModel.source(rhoCpField, T_);
  }
   
  TEqn.relax();

  fvConstraint.constrain(TEqn);

  if (cpsT.empty())
  {
    TEqn.solve(); 
    fvConstraint.constrain(T_); 
  }
  else
  {
    cpsT().insertEquation
    (
      T_.name(), 
      T_.name(),
      TEqn
    ); 
  }
  
  radiation_->correct();
  
  // Info<< "Min/max T: " << min(T_).value() << ' ' << max(T_).value() << endl;
}

}
 
// ************************************************************************* //
