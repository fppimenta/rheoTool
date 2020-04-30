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

#include "isotropicMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace fluidThermoMixtureModels
 {
   defineTypeNameAndDebug(isotropicMixture, 0);
   addToRunTimeSelectionTable(fluidThermoMixtureModel, isotropicMixture, dictionary);
 }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermoMixtureModels::isotropicMixture::isotropicMixture
(
  const word& name,
  const fvMesh& mesh,
  const word& phase1,
  const word& phase2
)
:
  fluidThermoMixtureModel(name, mesh, phase1, phase2),
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
  isViscDissipation_(readBool(dict_.lookup("enableViscousDissipation"))),  
  Cp1_(dict_.subDict(phase1).lookup("Cp")),
  Cp2_(dict_.subDict(phase2).lookup("Cp")),
  k1_0_(dimensionedScalar("k0",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.)),
  k1_1_(dimensionedScalar("k1",dimensionSet(1, 1, -3, -2, 0, 0, 0), 0.)),
  k1_2_(dimensionedScalar("k2",dimensionSet(1, 1, -3, -3, 0, 0, 0), 0.)),  
  k1_3_(dimensionedScalar("k3",dimensionSet(1, 1, -3, -4, 0, 0, 0), 0.)),
  k2_0_(dimensionedScalar("k0",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.)),
  k2_1_(dimensionedScalar("k1",dimensionSet(1, 1, -3, -2, 0, 0, 0), 0.)),
  k2_2_(dimensionedScalar("k2",dimensionSet(1, 1, -3, -3, 0, 0, 0), 0.)),  
  k2_3_(dimensionedScalar("k3",dimensionSet(1, 1, -3, -4, 0, 0, 0), 0.))
{ 
  // Fill-in kappaCoeffs
  scalar tmpV;
  dict_.subDict(phase1).subDict("kappaCoeffs").lookup("a0") >> tmpV; k1_0_.value() = tmpV;  
  dict_.subDict(phase1).subDict("kappaCoeffs").lookup("a1") >> tmpV; k1_1_.value() = tmpV;
  dict_.subDict(phase1).subDict("kappaCoeffs").lookup("a2") >> tmpV; k1_2_.value() = tmpV;
  dict_.subDict(phase1).subDict("kappaCoeffs").lookup("a3") >> tmpV; k1_3_.value() = tmpV;
  
  dict_.subDict(phase2).subDict("kappaCoeffs").lookup("a0") >> tmpV; k2_0_.value() = tmpV;  
  dict_.subDict(phase2).subDict("kappaCoeffs").lookup("a1") >> tmpV; k2_1_.value() = tmpV;
  dict_.subDict(phase2).subDict("kappaCoeffs").lookup("a2") >> tmpV; k2_2_.value() = tmpV;
  dict_.subDict(phase2).subDict("kappaCoeffs").lookup("a3") >> tmpV; k2_3_.value() = tmpV;
  
  // No need to initialize kappa here, as that is the first step before solving T (kappa
  // is only used in TEqn)     
}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluidThermoMixtureModels::isotropicMixture::correct
(
  const volVectorField& U,
  const surfaceScalarField& rhoPhi,
  const volSymmTensorField& tau,
  const volScalarField& alpha1,
  const volScalarField& rho,
  fv::options& fvOptions,
  int nNonOrtoC
)
{

int nNonOC(max(1,nNonOrtoC+1));

for (int jj = 0; jj<nNonOC; jj++)
{
   // Update kappa.
   // Note even if k1,k2,k3 == 0, alpha1 is always changing and kappa allways
   // needs to be updated 
  
   kappa_ = 
   alpha1*
   (      
      k1_0_ 
    + k1_1_*T_
    + k1_2_*T_*T_
    + k1_3_*T_*T_*T_
   )
   + (1.0-alpha1)*
   (      
      k2_0_ 
    + k2_1_*T_
    + k2_2_*T_*T_
    + k2_3_*T_*T_*T_
  );
 
  // TEqn
  volScalarField cpField(cp(alpha1)); 
  
  fvScalarMatrix TEqn
  (
      fvm::ddt(rho, T_)
    + fvm::div(rhoPhi, T_, "phi(rhoPhi,T)")
    - (1.0/cpField)*fvm::laplacian(kappa_, T_)
  );

  if (fvOptions.appliesToField(T_.name()))
  { 
    TEqn -= fvOptions(rho, T_);
  }
  
  if (isViscDissipation_)
    TEqn -= (tau && fvc::grad(U))/cpField;
  
  TEqn.relax();

  fvOptions.constrain(TEqn);

  TEqn.solve(); 
 
  fvOptions.correct(T_);
  
  //Info<< "Min/max T: " << min(T_).value() << ' ' << max(T_).value() << endl;  
}

}
 

// ************************************************************************* //
