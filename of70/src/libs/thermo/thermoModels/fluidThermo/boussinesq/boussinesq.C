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

#include "boussinesq.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace fluidThermoModels
 {
   defineTypeNameAndDebug(boussinesq, 0);
   addToRunTimeSelectionTable(fluidThermoModel, boussinesq, dictionary);
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermoModels::boussinesq::boussinesq
(
  const word& name,
  const fvMesh& mesh
)
:
  fluidThermoModel(name, mesh),
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
    dimensionedScalar("kt", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.),
    extrapolatedCalculatedFvPatchField<scalar>::typeName
  ),
  beta_("beta", dimless/dimTemperature, .0), 
  TRef_("TRef", dimTemperature, .0), 
  isNatConvection_(true),
  isViscDissipation_(readBool(dict_.lookup("enableViscousDissipation"))),  
  radiation_(radiationModel::New(T_)),
  rhoCp_(dict_.lookup("rhoCp")),
  g_
  (
   IOobject
   (
     "g",
     mesh.time().constant(),
     mesh,
     IOobject::MUST_READ,
     IOobject::NO_WRITE
   )
  ),
  hRef_
  (
   IOobject
   (
     "hRef",
     mesh.time().constant(),
     mesh,
     IOobject::READ_IF_PRESENT,
     IOobject::NO_WRITE
   ),
   dimensionedScalar("hRef", dimLength, 0)
  ),
  k0_(dimensionedScalar("k0",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.)),
  k1_(dimensionedScalar("k1",dimensionSet(1, 1, -3, -2, 0, 0, 0), 0.)),
  k2_(dimensionedScalar("k2",dimensionSet(1, 1, -3, -3, 0, 0, 0), 0.)),  
  k3_(dimensionedScalar("k3",dimensionSet(1, 1, -3, -4, 0, 0, 0), 0.))
{
  if (mag(g_.value()) == 0)
    isNatConvection_ = false;

  // Natural convection
  if (isNatConvection_)
  {
     dict_.lookup("beta") >> beta_;
     dict_.lookup("TRef") >> TRef_;
  }
  
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

void Foam::fluidThermoModels::boussinesq::correct
(
  const volVectorField& U,
  const surfaceScalarField& phi,
  const volSymmTensorField& tau,
  autoPtr<coupledSolver>& cpsT,
  fv::options& fvOptions,
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
    + fvm::div(phi*rhoCp_, T_)
    - fvm::laplacian(kappa_, T_)
    ==
      rhoCp_*radiation_->ST(rhoCp_, T_) 
  );
  

  if (fvOptions.appliesToField(T_.name()))
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
     
    TEqn -= fvOptions(rhoCpField, T_);
  }
  
  if (isViscDissipation_)
   TEqn -= (tau && fvc::grad(U));
  
  TEqn.relax();

  fvOptions.constrain(TEqn);
 
  if (cpsT.empty())
  {
    TEqn.solve(); 
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

  fvOptions.correct(T_);
  
  //Info<< "Min/max T: " << min(T_).value() << ' ' << max(T_).value() << endl;  
}

}
 

// ************************************************************************* //
