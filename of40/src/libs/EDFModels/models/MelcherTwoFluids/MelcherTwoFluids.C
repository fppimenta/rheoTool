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

#include "MelcherTwoFluids.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MelcherTwoFluids, 0);
    addToRunTimeSelectionTable(EDFEquation, MelcherTwoFluids, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MTFPhase::MTFPhase
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:  
    relPermi_(dict.lookup("relPerm")),
    sigmai_(dict.lookup("sigma")),
    phaseName_(name)
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MelcherTwoFluids::MelcherTwoFluids
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    EDFEquation(name, phi),
    rhoE_
    (
        IOobject
        (
            "rhoE" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    psi_
    (
        IOobject
        (
            "psi" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    relPerm_
    (
        IOobject
        (
            "relPerm" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar
        (
                "1",
                dimless,
                1.
        ),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    ),
    sigma_
    (
        IOobject
        (
            "sigma" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar
        (
                "1",
                dimensionSet(-1, -3, 3, 0, 0, 2, 0),
                1.
        ),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    ),
    extraE_(dict.lookupOrDefault<dimensionedVector>("extraEField", dimensionedVector("0", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero))),
    psiEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    rhoEEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("rhoEEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    maxIterPsi_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<scalar>("maxIter", 50)),
    maxIterRhoE_(phi.mesh().solutionDict().subDict("electricControls").subDict("rhoEEqn").lookupOrDefault<scalar>("maxIter", 50)),
    phases_(),
    nPhases_(2)
{
    PtrList<entry> specEntries(dict.lookup("phases"));
    
   //- Check if exactly 2 phases exist 
    if (specEntries.size()!=2)
     {
        FatalErrorIn("MelcherTwoFluids::MelcherTwoFluids")
                        << " Two phases must be specified in subDict 'phases' of dictionary 'electricProperties'. " 
                        << exit(FatalError);
     }
     
    phases_.setSize(nPhases_);
         
    forAll (phases_, phaseI)
    {    
        phases_.set
        (
            phaseI,
            new MTFPhase
            (
                specEntries[phaseI].keyword(),
                phi,
                specEntries[phaseI].dict()
            )
        );        
    }          
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::MelcherTwoFluids::Fe() const
{

    return
    (
         fvc::div(
                    relPerm_*epsilonK_*
                    (
                      fvc::grad(psi_)*fvc::grad(psi_) - 0.5*magSqr(fvc::grad(psi_))*Foam::tensor::I
                    ),
                   "div(MaxwellT)" 
                 ) 
    );
/*  

    return
    (
         -rhoE_*fvc::grad(psi_) - 0.5 * magSqr(fvc::grad(psi_)) * epsilonK_*fvc::grad( relPerm_ )
    );
*/
          
}

void Foam::MelcherTwoFluids::correct()
{

       correctAlphaProperties();
 
    //- Equation for the internal potential
       
       scalar res=GREAT; 
       scalar iter=0; 

       while (res > psiEqRes_ && iter < maxIterPsi_)
         { 

		fvScalarMatrix psiEqn
		(	  
		     fvm::laplacian(relPerm_*epsilonK_, psi_, "laplacian(eps,psi)")
		  == 
		    -rhoE_		           		    
		);
	        
	        psiEqn.relax();
		res=psiEqn.solve().initialResidual();

		iter++;
        } 
        
    //- Equation for the charge transport
   
       res=GREAT; 
       iter=0;

       while (res > rhoEEqRes_ && iter < maxIterRhoE_)
          { 
          
		fvScalarMatrix rhoEEqn
		(
		     fvm::ddt(rhoE_)
		   + fvm::div(phi(), rhoE_)
		   ==
		     fvc::laplacian(sigma_, psi_) 
		);
		
		rhoEEqn.relax();
		res=rhoEEqn.solve().initialResidual();
		
		iter++;
          } 

}

void Foam::MelcherTwoFluids::correctAlphaProperties()
{
    const volScalarField& alpha1(phi().mesh().lookupObject<volScalarField>(phases_[0].phaseName() ) );
    
    const volScalarField bAlpha1
    (
        min(max(alpha1, scalar(0)), scalar(1))
    );
     
    relPerm_ = bAlpha1*phases_[0].relPermi() + (1.-bAlpha1)*phases_[1].relPermi();
    //relPerm_ = 1./( bAlpha1/(phases_[0].relPermi() ) + (1.-bAlpha1)/(phases_[1].relPermi() ) );
    
    dimensionedScalar mSmall("mSmall", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 1e-20);
    
    sigma_ = 1./( bAlpha1/(phases_[0].sigmai() + mSmall) + (1.-bAlpha1)/(phases_[1].sigmai() + mSmall) ); 
  //  sigma_ = bAlpha1*phases_[0].sigmai() + (1.-bAlpha1)*phases_[1].sigmai(); 
   
}

// ************************************************************************* //
