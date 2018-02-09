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

#include "Ohmic.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Ohmic, 0);
    addToRunTimeSelectionTable(EDFEquation, Ohmic, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Ohmic::OhSpecie::OhSpecie
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:  
    zi_(dict.lookup("z")),
    Di_(dict.lookup("D")),
    namei_(name) 
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Ohmic::Ohmic
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    EDFEquation(name, phi),
    sigma_
    (
        IOobject
        (
            "sigma" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    phiE_
    (
        IOobject
        (
            "phiE" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    relPerm_(dict.lookup("relPerm")),
    Deff_("0", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.),
    extraE_(dict.lookupOrDefault<dimensionedVector>("extraEField", dimensionedVector("0", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero))),
    sigmaEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("sigmaEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    phiEEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    maxIterSigma_(phi.mesh().solutionDict().subDict("electricControls").subDict("sigmaEqn").lookupOrDefault<int>("maxIter", 50)),
    maxIterPhiE_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<int>("maxIter", 50)),
    species_(),
    nSpecies_(0)
{
    PtrList<entry> specEntries(dict.lookup("species"));
    nSpecies_ = specEntries.size();
    species_.setSize(nSpecies_);
      
    forAll (species_, specI)
    {
    
        species_.set
        (
            specI,
            new OhSpecie
            (
                specEntries[specI].keyword(),
                phi,
                specEntries[specI].dict()
            )
        );   
    }   
    
    // Compute Deff_ (assuming two species only)
    scalar z0 = mag(species_[0].zi().value());
    scalar z1 = mag(species_[1].zi().value());
       
    Deff_ =  2*( species_[0].Di()*species_[1].Di() )
           / ( species_[0].Di() + species_[1].Di() );  
           
     
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::Ohmic::Fe() const
{
    return
    (
         fvc::laplacian(epsilonK_*relPerm_, phiE_, "rhoE") * ( fvc::grad(phiE_) + extraE_) 
    );     
}

void Foam::Ohmic::correct()
{

       scalar res=GREAT; 
       scalar iter=0;  
  
   //- Equation for the conductivity 
       while (res > sigmaEqRes_ && iter < maxIterSigma_)
          { 
          
		fvScalarMatrix sigmaEqn
		(
		    fvm::ddt(sigma_)
		  + fvm::div(phi(), sigma_)
		  ==
		    fvm::laplacian(Deff_, sigma_, "laplacian(Deff,sigma)") 
		     
		);
		
		sigmaEqn.relax();
		res=sigmaEqn.solve().initialResidual();
		
		iter++;
          } 
  
   //- Equation for the current (rhoE)       
       res=GREAT;
       iter=0;  
   
       while (res > phiEEqRes_ && iter < maxIterPhiE_)
         { 

		fvScalarMatrix phiEEqn
		(	  
		     fvm::laplacian(sigma_, phiE_)  		    
		);
	        
	        phiEEqn.relax();
		res=phiEEqn.solve().initialResidual();

		iter++;
        } 

}

// ************************************************************************* //
