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

#include "Oldroyd_BSqrt.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs 
{
    defineTypeNameAndDebug(Oldroyd_BSqrt, 0);
    addToRunTimeSelectionTable(constitutiveEq, Oldroyd_BSqrt, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveEqs::Oldroyd_BSqrt::Oldroyd_BSqrt
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    constitutiveEq(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    b_
    (
        IOobject
        (
            "b" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    rho_(dict.lookup("rho")),
    etaS_(dict.lookup("etaS")),
    etaP_(dict.lookup("etaP")),
    lambda_(dict.lookup("lambda"))
{
 checkForStab(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveEqs::Oldroyd_BSqrt::correct()
{
    // Velocity gradient tensor
    volTensorField L = fvc::grad(U());
   
    volTensorField a(L*0.0);
    scalar t1(0.0), t2(0.0), t3(0.0), B1(0.0), B2(0.0), B3(0.0), w1(0.0), w2(0.0), w3(0.0), D(0.0);
    
    forAll(b_, cellI)
      {
 
       t1 = b_[cellI].yy()+b_[cellI].zz();
       t2 = b_[cellI].xx()+b_[cellI].zz(); 
       t3 = b_[cellI].xx()+b_[cellI].yy();

       B1 = b_[cellI].yz();  B2 = b_[cellI].xz();  B3 = b_[cellI].xy();

       w1 =  b_[cellI].xy()*L[cellI].xx() - b_[cellI].xx()*L[cellI].xy()
           + b_[cellI].yy()*L[cellI].yx() - b_[cellI].xy()*L[cellI].yy()
           + b_[cellI].yz()*L[cellI].zx() - b_[cellI].xz()*L[cellI].zy();

       w2 =  b_[cellI].xz()*L[cellI].xx() - b_[cellI].xx()*L[cellI].xz()
           + b_[cellI].zz()*L[cellI].zx() - b_[cellI].xz()*L[cellI].zz()
           + b_[cellI].yz()*L[cellI].yx() - b_[cellI].xy()*L[cellI].yz();

       w3 =  b_[cellI].xz()*L[cellI].xy() - b_[cellI].xy()*L[cellI].xz()
           + b_[cellI].yz()*L[cellI].yy() - b_[cellI].yy()*L[cellI].zy()
           + b_[cellI].zz()*L[cellI].zy() - b_[cellI].yz()*L[cellI].zz();

       D = t1*(t2*t3-B1*B1) - B2*(B2*t2+B1*B3) - B3*(B2*B1 + B3*t3);

       a[cellI].xy() = (  (t1*t2-B3*B3)*w1 - (B1*t1+B3*B2)*w2 + (B2*t2+B1*B3)*w3 )/D;
       a[cellI].xz() = ( -(B1*t1+B3*B2)*w1 + (t1*t3-B2*B2)*w2 - (B2*B1+B3*t3)*w3 )/D;
       a[cellI].yz() = (  (B2*t2+B1*B3)*w1 - (B2*B1+B3*t3)*w2 + (t2*t3-B1*B1)*w3 )/D;
 
       a[cellI].yx() = -a[cellI].xy();
       a[cellI].zx() = -a[cellI].xz();
       a[cellI].zy() = -a[cellI].yz();

    }     


    // Stress transport equation
    fvSymmTensorMatrix bEqn
    (
         fvm::ddt(b_)
       + fvm::div(phi(), b_) 
     ==
         symm((b_&L) + (a & b_))
       + (0.5/lambda_) * inv( b_.T() ) 
       - fvm::Sp(0.5/lambda_, b_)
    );
 
    bEqn.relax();
    bEqn.solve();

    // Convert from b to tau

    dimensionedSymmTensor Itensor
    ( 
        "Identity", 
        dimensionSet(0, 0, 0, 0, 0, 0, 0), 
        symmTensor::I 
    );

    tau_ = (etaP_/lambda_) * symm( (b_&b_) - Itensor);

    tau_.correctBoundaryConditions();
 
}


// ************************************************************************* //
