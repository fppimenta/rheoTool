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

#include "Newtonian.H"
#include "addToRunTimeSelectionTable.H"
#include "Constant.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(Newtonian, 0);
    addToRunTimeSelectionTable(constitutiveEq, Newtonian, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::Newtonian::Newtonian
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("eta")),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor
        (
                "0",
                dimensionSet(1, -1, -2, 0, 0, 0, 0),
                pTraits<symmTensor>::zero  
        ),
        extrapolatedCalculatedFvPatchField<symmTensor>::typeName
    ),
    thermoEtaPtr_(thermoFunction::New("thermoEta", U.mesh(), dict)),   
    eta_
    (
        IOobject
        (
            "eta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            thermoEtaPtr_().type() == "Constant" ? IOobject::NO_WRITE : IOobject::AUTO_WRITE 
        ),
        U.mesh(),
        etaS_,
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    )    
{
  // Update eta for non-isothermal cases
  correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
Foam::tmp<Foam::fvVectorMatrix> Foam::constitutiveEqs::Newtonian::divTau
(
  const volVectorField& U
) const
{
  return
  (
    fvm::laplacian( etaS_/rho_, U, "laplacian(eta,U)")
  );
}

void Foam::constitutiveEqs::Newtonian::correct()
{
  // Update viscosity for non-isothermal flow (does nothing when
  // thermoEtaPtr_ is the isothermal function)
  thermoEtaPtr_->multiply(eta_,etaS_);
}


// ************************************************************************* //
