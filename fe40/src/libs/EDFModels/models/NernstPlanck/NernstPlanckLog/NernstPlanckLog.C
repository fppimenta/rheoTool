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

#include "NernstPlanckLog.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NernstPlanckLog, 0);
    addToRunTimeSelectionTable(EDFEquation, NernstPlanckLog, dictionary);
}
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NPLSpecie::NPLSpecie
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:  
    logCi_
    (
        IOobject
        (
            name + "Log",
            phi.time().timeName(),
            phi.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phi.mesh()
    ),
    ci_
    (
        IOobject
        (
            name,
            phi.time().timeName(),
            phi.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Foam::exp(logCi_) * dimensionedScalar("1", dimensionSet(0, -3, 0, 0, 1, 0, 0), 1.0)
    ),
    zi_(dict.lookup("z")),
    Di_(dict.lookup("D")),
    namei_(name)  
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NernstPlanckLog::NernstPlanckLog
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
        zeroGradientFvPatchField<scalar>::typeName
    ),
    relPerm_(dict.lookup("relPerm")),
    T_(dict.lookup("T")),
    extraE_(dict.lookupOrDefault<dimensionedVector>("extraEField", dimensionedVector("0", dimensionSet(1, 1, -3, 0, 0, -1, 0), vector::zero))),
    psiContrib_(dict.lookupOrDefault<bool>("psiContrib", true)),
    phiEEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    psiEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    ciEqRes_(phi.mesh().solutionDict().subDict("electricControls").subDict("ciEqn").lookupOrDefault<scalar>("residuals", 1e-7)),
    maxIterPhiE_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<scalar>("maxIter", 50)),
    maxIterPsi_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<scalar>("maxIter", 50)),
    maxIterCi_(phi.mesh().solutionDict().subDict("electricControls").subDict("ciEqn").lookupOrDefault<scalar>("maxIter", 50)),
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
            new NPLSpecie
            (
                specEntries[specI].keyword(),
                phi,
                specEntries[specI].dict()
            )
        );  
    }     
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::NernstPlanckLog::Fe() const
{
    volScalarField rhoE( phiE_ * dimensionedScalar("norm", epsilonK_.dimensions()/dimArea, 0.) );
    
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

void Foam::NernstPlanckLog::correct()
{

for (int k=0; k<1; k++)
 {

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
           
    //- Equation for each ionic specie 
    
        dimensionedScalar icUnit("one", dimensionSet(0, 3, 0, 0, -1, 0, 0), 1.0);   
    
        forAll (species_, i)
         {
         
           res=GREAT;
           iter=0;
           
           volScalarField& ci(species_[i].ci());
           volScalarField& lci(species_[i].logCi());
           
           dimensionedScalar ca( species_[i].Di() * eK_ * species_[i].zi() / (kbK_*T_) ); 
           
           surfaceScalarField eMigFluxp( fvc::snGrad(phiE_+psi_) * phi().mesh().magSf() ); // Compute once, outside the loop                
           
           while (res > ciEqRes_ && iter < maxIterCi_)
            { 
                 
                   surfaceScalarField cif(Foam::exp(fvc::interpolate(lci))); 
           //      surfaceScalarField cif(fvc::interpolate(ci*icUnit));          
           

		 fvScalarMatrix ciLogEqn 
		  (
		    fvm::ddt(ci*icUnit, lci)
		 // + fvc::div(phi(), ci*icUnit, "div(phi,ciLog)")
		  ==  
		    fvm::laplacian(species_[i].Di()*cif, lci, "laplacian(D,ci)")   	  
		//  + fvc::laplacian(EMigD, phiE_+psi_, "laplacian(EMigD,phiE)")
		  + fvc::div((eMigFluxp*ca - phi()) * cif)
		    
		  );        
		
		  ciLogEqn.relax(phi().mesh().solutionDict().equationRelaxationFactor("ciLog"));
		  res=ciLogEqn.solve(phi().mesh().solutionDict().solver("ciLog")).initialResidual();

                  ci = Foam::exp(lci)/icUnit;
                  ci.correctBoundaryConditions();
                  
		  iter++;
           } 
       }
        
      Info << nl << endl;
      
 } // for del
}

// ************************************************************************* //
