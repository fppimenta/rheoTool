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
 
#include "NernstPlanckCoupled.H"
#include "addToRunTimeSelectionTable.H"

#include "zeroGradientFvPatchField.H"
#include "zeroIonicFluxFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
 namespace EDFEquations
 {
    defineTypeNameAndDebug(NernstPlanckCoupled, 0);
    addToRunTimeSelectionTable(EDFEquation, NernstPlanckCoupled, dictionary);
 }
} 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquations::NernstPlanckCoupled::NPSpecie::NPSpecie
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
{ 

 // If any zeroIonicFlux BC is set for ions, we need to tell them 
 // that the PNP is coupled, since the coeffs handling depends on that. 
 forAll(ci_.boundaryField(), pI)
   if (isType<zeroIonicFluxFvPatchScalarField>(ci_.boundaryField()[pI])) 
     const_cast<zeroIonicFluxFvPatchScalarField&>
     (
        refCast<const zeroIonicFluxFvPatchScalarField>
        ( 
          ci_.boundaryField()[pI]
        )
     ).setFlagCoupled(true);    
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::EDFEquations::NernstPlanckCoupled::NernstPlanckCoupled
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
 nSpecies_(0),
 cps_(NULL),
 solveAlone_(false)
{
 // Build species
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
 
 // Coupled solver
 solveAlone_ = !readBool(phi.mesh().solutionDict().subDict("coupledSolvers").subDict("ciPsi").lookup("solveWithUptau"));

 // Insert fields, either in the own system or in the one already existing for momentum
 if (solveAlone_)
 {
   cps_.reset
   (
     new coupledSolver("ciPsi", phi.mesh())
   );
 
   // Add the variables to the solver context
   // psi
   cps_->insertField(psi_);

   // ci's
   forAll(species_, si)
     cps_->insertField(species_[si].ci());
 }
 else
 {
   // Check if momentum is coupled
   bool solveUptauEnab_ = readBool(phi.mesh().solutionDict().subDict("coupledSolvers").subDict("Uptau").lookup("solveCoupledUp"));
   if (!solveUptauEnab_)
   {
     FatalErrorIn("Foam::EDFEquations::NernstPlanckCoupled::NernstPlanckCoupled\n")
     << "Cannot solve electric module coupled with momentum/continuity if "
     << "solveCoupledUp flag is set to false in fvSolution." 
     << abort(FatalError);
   }
   
   coupledSolver& cps = const_cast<coupledSolver&>
   (
     phi.mesh().lookupObject<coupledSolver>("Uptau")
   );
 
   // Add the variables to the solver context
   // psi
   cps.insertField(psi_);

   // ci's
   forAll(species_, si)
    cps.insertField(species_[si].ci()); 
 }
 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::EDFEquations::NernstPlanckCoupled::Fe() const
{    
 volScalarField rhoE
 (
    IOobject
    (
      "rhoE",
      psi_.time().constant(),
      psi_.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
    ), 
    psi_.mesh(),
    dimensionedScalar("0", psi_.dimensions()*epsilonK_.dimensions()/dimArea, 0.)
 );
    
 forAll (species_, i)
   rhoE+= species_[i].zi()*species_[i].ci()*FK_;
    
 if (solvePhiE_)
 {
   if (psiContrib_)
     return (solveAlone_ ? -rhoE * (fvc::grad(phiE_+psi_) - extraE_) : -rhoE * (fvc::grad(phiE_) - extraE_) );  
  
   else
     return ( -rhoE * (fvc::grad(phiE_) - extraE_)   );  
 }
 else
 {
   if (psiContrib_)
     return (solveAlone_ ? -rhoE * (fvc::grad(psi_) - extraE_) : rhoE * extraE_  );
    
   else
     return (  -rhoE * (-extraE_)  );
 }      
}
 
void Foam::EDFEquations::NernstPlanckCoupled::FeImplCoupled
(
  const dimensionedScalar rho
) 
const
{
 if (!solveAlone_ && psiContrib_)
 {
   coupledSolver& cps = psi_.mesh().lookupObjectRef<coupledSolver>("Uptau");
      
   volScalarField rhoE
   (
    IOobject
    (
      "rhoE",
      psi_.time().timeName(),
      psi_.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE,
      false
    ),
    psi_.mesh(),
    dimensionedScalar("0", psi_.dimensions()*epsilonK_.dimensions()/(dimArea*rho.dimensions()), 0.)
   );
    
   forAll(species_, i)
     rhoE += species_[i].zi()*species_[i].ci()*FK_/rho; // Division by rho in momentum equation
   
   cps.insertEquation
   (
     "U",
     psi_.name(),
     fvmb::grad(rhoE, psi_)  
   );
 }
}

void Foam::EDFEquations::NernstPlanckCoupled::correct()
{  
 if (solveAlone_)
  solveAlone();
 else
  solveWithUptau();
}

void Foam::EDFEquations::NernstPlanckCoupled::solveAlone()
{
 //- phiE is never coupled to other variables. It can be solved alone.  
 scalar res=GREAT; 
 scalar iter=0;  

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
 
 //- Solve the coupled system of psi and ci's 
 for (int j=0; j<nIterPNP_; j++) 
 {
   Info << "PNP Coupling iteration: " << j << endl;
  
   // Insert the equations in the system
  
   // psi_psi
   cps_->insertEquation
   (
    psi_.name(),
    psi_.name(),
    fvm::laplacian(psi_)
   ); 
 
   // Create alias of psi to ensure zero-flux on the boundaries 
   // where it should hold. psiCorr only differs from psi on
   // zero-flux boundaries, where it receives a zero-gradient.
   #include "createPsiCorr.H"

   forAll(species_, i)
   {
     volScalarField& ci = species_[i].ci();            
     dimensionedScalar cf(species_[i].Di() * eK_ * species_[i].zi() / (kbK_*T_));
  
     // Change psiCorr BC according to the BC's of ci
     #include "correctPsiCorr.H"

     // psi_ci
     cps_->insertEquation
     (
       psi_.name(),
       ci.name(),  
       fvm::Sp(species_[i].zi()*FK_/(relPerm_*epsilonK_), ci)	
     ); 
     
     // ci_ci (allow under-relaxation)
     fvScalarMatrix ciEqn
     (
       fvm::ddt(ci)
     + fvm::div(phi(), ci, "div(phi,ci)")  
     ==  
       fvm::laplacian(species_[i].Di(), ci, "laplacian(D,ci)") 
     + fvc::laplacian(ci*cf, phiE_, "laplacian(elecM)")	
     );
  
     if (phi().mesh().relaxEquation("ci"))
       ciEqn.relax(phi().mesh().equationRelaxationFactor("ci"));
     
     
     cps_->insertEquation
     (
       ci.name(),
       ci.name(),
       ciEqn
     );
  
     // ci_psi
     cps_->insertEquation
     (
       ci.name(),
       psi_.name(),
       fvm::laplacian(-ci*cf, psiCorr, "laplacian(elecM)")	 	 	
     );   
   }
 
   // Solve the coupled system
   cps_->solve();  
   
   // Ensure positive concentration 
   forAll(species_, i)
   {
    volScalarField& ci = species_[i].ci();   
    ci = Foam::max( dimensionedScalar("lowerLimit",ci.dimensions(), 0.), ci );
    ci.correctBoundaryConditions();
   }  
 } 
}

void Foam::EDFEquations::NernstPlanckCoupled::solveWithUptau()
{
 //- Note: phiE does contribute to momentum equation, but since 
 // it results from a Laplace equation, the field becomes stationary
 // upon convergence of the equation and its contribution to momentum
 // becomes a stationary source term. Thus, no need for implicit coupling.  
 scalar res=GREAT; 
 scalar iter=0;  

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
  
 // Lookup for the system
 coupledSolver& cps = phi().mesh().lookupObjectRef<coupledSolver>("Uptau");
 
 // Note: we cannot inner iterate the electric module only, as
 // we do in segregated solving.

 // Insert the equations in the system
  
 // psi_psi
 cps.insertEquation
 (
  psi_.name(),
  psi_.name(),
  fvm::laplacian(psi_)
 ); 
 
 // Create alias of psi to ensure zero-flux on the boundaries 
 // where it should hold. psiCorr only differs from psi on
 // zero-flux boundaries, where it receives a zero-gradient.
 #include "createPsiCorr.H"

 forAll(species_, i)
 {
   volScalarField& ci = species_[i].ci();            
   dimensionedScalar cf(species_[i].Di() * eK_ * species_[i].zi() / (kbK_*T_));
 
   // Change psiCorr BC according to the BC's of ci
   #include "correctPsiCorr.H"

   // psi_ci
   cps.insertEquation
   (
     psi_.name(),
     ci.name(),  
     fvm::Sp(species_[i].zi()*FK_/(relPerm_*epsilonK_), ci)	
   ); 
  
   // ci_ci
   cps.insertEquation
   (
     ci.name(),
     ci.name(),
     fvm::ddt(ci)
   + fvm::div(phi(), ci, "div(phi,ci)")  
   ==  
     fvm::laplacian(species_[i].Di(), ci, "laplacian(D,ci)") 
   + fvc::laplacian(ci*cf, phiE_, "laplacian(elecM)")	
   );
  
   // ci_psi
   cps.insertEquation
   (
     ci.name(),
     psi_.name(),
     fvm::laplacian(-ci*cf, psiCorr, "laplacian(elecM)")	 	 	
   );   
 }
 
 // NOTE: when all vars are solved coupled, we do not force [0,+Inf]
 // bound on ci. Negative ci could arise in harsh conditions.
}
 
// ************************************************************************* //
