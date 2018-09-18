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

#include "NernstPlanck.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace EDFEquations
{
    defineTypeNameAndDebug(NernstPlanck, 0);
    addToRunTimeSelectionTable(EDFEquation, NernstPlanck, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquations::NernstPlanck::NPSpecie::NPSpecie
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:  
    ci_
    (
        IOobject
        (
            name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    zi_(dict.lookup("z")),
    Di_(dict.lookup("D")),
    namei_(name)  
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquations::NernstPlanck::NernstPlanck
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    EDFEquation(name, phi),
    solvePhiE_(checkForPhiE(name, phi)),
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
    phiE_
    (
        IOobject
        (
            "phiE" + name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::READ_IF_PRESENT,
            solvePhiE_ == false ? (IOobject::NO_WRITE) : (IOobject::AUTO_WRITE)
        ),
        phi.mesh(),
        dimensionedScalar
        (
                "zero",
                psi_.dimensions(),
                pTraits<scalar>::zero
        ),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
    ),
    relPerm_(dict.lookup("relPerm")),
    T_(dict.lookup("T")),
    extraE_(dict.lookupOrDefault<dimensionedVector>("extraEField", dimensionedVector("0", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero))),
    psiContrib_(dict.lookupOrDefault<bool>("psiContrib", true)),
    phiEEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    psiEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    ciEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("ciEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    maxIterPhiE_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<int>("maxIter", 50)),
    maxIterPsi_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<int>("maxIter", 50)),
    maxIterCi_(phi.mesh().solutionDict().subDict("electricControls").subDict("ciEqn").lookupOrDefault<int>("maxIter", 50)),
    nIterPNP_(phi.mesh().solutionDict().subDict("electricControls").lookupOrDefault<int>("nIterPNP", 2)),
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
            new NPSpecie
            (
                specEntries[specI].keyword(),
                phi,
                specEntries[specI].dict()
            )
        );  
    }     
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::EDFEquations::NernstPlanck::Fe() const
{    
    volScalarField rhoE( psi_ * dimensionedScalar("norm", epsilonK_.dimensions()/dimArea, 0.) );
    
    forAll (species_, i)
    {
      rhoE 
         += (
               species_[i].zi()*species_[i].ci()*FK_
            );
    }
    
    if (solvePhiE_)
     {
       if (psiContrib_)
        {
           return
           (
             -rhoE * ( fvc::grad(phiE_+psi_) - extraE_) 
           );
        }
       else
        {
           return
           (
             -rhoE * ( fvc::grad(phiE_) - extraE_) 
           ); 
        }   
     }
    else
     {
       if (psiContrib_)
        {
           return
           (
             -rhoE * ( fvc::grad(psi_) - extraE_) 
           );
        }
       else
        {
           return
           (
             -rhoE * (-extraE_) 
           ); 
        }  
     }      
}

void Foam::EDFEquations::NernstPlanck::correct()
{

// Electrokinetic coupling loop
 
for (int j=0; j<nIterPNP_; j++) 
 {
       Info << "PNP Coupling iteration: " << j << endl;
       
       scalar res=GREAT; 
       scalar iter=0;  
  
   //- Equation for the external potential (loop for the case
   //  of non-orthogonal grids)   
       if (solvePhiE_)
       {
         while (res > phiEEqRes_ && iter < maxIterPhiE_)
            { 
          
	   	fvScalarMatrix phiEEqn
		(
		    fvm::laplacian(phiE_)
		);
		
		phiEEqn.relax();
		res=phiEEqn.solve().initialResidual();
		
		iter++;
            }
       }
    
    //- Equation for the intrinsic potential
   
       res=GREAT;
       iter=0; 
   
       volScalarField souE(psi_ * dimensionedScalar("norm1",dimless/dimArea,0.));       

       forAll (species_, i)
         {
           souE 
             += (
                    -species_[i].zi()*species_[i].ci()*FK_
                   /(relPerm_*epsilonK_)
                );
         }

        while (res > psiEqRes_ && iter < maxIterPsi_)
         { 

		fvScalarMatrix psiEqn
		(	  
		     fvm::laplacian(psi_)
		  ==
		     souE		           		    
		);
	     
	        psiEqn.relax();
		res=psiEqn.solve().initialResidual();

		iter++;
        } 
  
        
    //- Nernst-Planck equation for each ionic specie  
    
      //  surfaceScalarField eMigFluxp( fvc::snGrad(phiE_+psi_) * phi().mesh().magSf() ); // Compute once, outside the loop 
        
        forAll (species_, i)
         {         
           res=GREAT;
           iter=0;
           
           volScalarField& ci = species_[i].ci();
                     
           dimensionedScalar cf(species_[i].Di() * eK_ * species_[i].zi() / (kbK_*T_));  
           
           while (res > ciEqRes_ && iter < maxIterCi_)
            { 
       
		 fvScalarMatrix ciEqn 
		  (
		   fvm::ddt(ci)
		+  fvm::div(phi(), ci, "div(phi,ci)")  
		  ==  
		    fvm::laplacian(species_[i].Di(), ci, "laplacian(D,ci)") 	 
              //    + fvm::div(eMigFluxp*cf, ci, "div(eMigFlux,ci)") 
                  + fvc::laplacian(ci*cf, phiE_+psi_, "laplacian(elecM)")	 
		  );      
		  	  
		  ciEqn.relax(phi().mesh().equationRelaxationFactor("ci"));
		  res=ciEqn.solve(phi().mesh().solver("ci")).initialResidual();

                  ci = Foam::max( dimensionedScalar("lowerLimit",ci.dimensions(), 0.), ci );
                  ci.correctBoundaryConditions();
                  
		  iter++;
           } 
                      
       }
       
       
 }       
      
 
}

// ************************************************************************* //
