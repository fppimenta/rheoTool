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

#include "DebyeHuckel.H"
#include "addToRunTimeSelectionTable.H"
 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DebyeHuckel, 0);
    addToRunTimeSelectionTable(EDFEquation, DebyeHuckel, dictionary);
} 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DebyeHuckel::DHSpecie::DHSpecie
(
    const word& name,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:  
    zi_(dict.lookup("z")),
    c0i_(dict.lookup("c0")),
    namei_(name) 
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DebyeHuckel::DebyeHuckel
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
    maxIterPhiE_(phi.mesh().solutionDict().subDict("electricControls").subDict("phiEEqn").lookupOrDefault<int>("maxIter", 50)),
    maxIterPsi_(phi.mesh().solutionDict().subDict("electricControls").subDict("psiEqn").lookupOrDefault<int>("maxIter", 50)),
    species_(),
    nSpecies_(0)
{
    PtrList<entry> specEntries(dict.lookup("species"));
    nSpecies_ = specEntries.size();
    species_.setSize(nSpecies_);
    
    Info << nl << endl;
    dimensionedScalar cum("0", dimensionSet(0, -3, 0, 0, 1, 0, 0), 0);   
    forAll (species_, specI)
    {
    
        species_.set
        (
            specI,
            new DHSpecie
            (
                specEntries[specI].keyword(),
                phi,
                specEntries[specI].dict()
            )
        );  
        
        cum += species_[specI].DebyeLengthP(relPerm_, T_); 
    }  
    
    dimensionedScalar DebL
     (
       sqrt
           (
              Foam::EDFEquation::epsilonK_* relPerm_ * Foam::EDFEquation::kbK_*T_     
            / (cum*Foam::EDFEquation::FK_*Foam::EDFEquation::eK_)
           )
     );
     
    Info << "Debye length: " << DebL.value() << " m." << nl << endl; 
     
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::DebyeHuckel::Fe() const
{
    
    volScalarField rhoE( psi_ * dimensionedScalar("norm", epsilonK_.dimensions()/dimArea, 0.) );
    
    forAll (species_, i)
    {
      rhoE 
         += (
               species_[i].zi()*species_[i].c0i()*FK_
             * ( 1. - eK_*species_[i].zi()*psi_/(kbK_*T_) )
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

void Foam::DebyeHuckel::correct()
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
   
       dimensionedScalar souE("SouE", psi_.dimensions()/dimArea, 0.);       
       dimensionedScalar souI("SouE", dimless/dimArea, 0.);
       
       forAll (species_, i)
         {
           souE 
             += (
                     species_[i].zi()*species_[i].c0i()*FK_
                   /(relPerm_*epsilonK_)
                );
                
           souI 
             += (
                   eK_*species_[i].zi()*species_[i].zi()*species_[i].c0i()*FK_
                   /(relPerm_*epsilonK_*kbK_*T_)
                );
         }

        while (res > psiEqRes_ && iter < maxIterPsi_)
         { 

		fvScalarMatrix psiEqn
		(	  
		     fvm::laplacian(psi_)
		  == 
		     fvm::Sp( souI, psi_)
		   - souE
		           		    
		);
	        
	        psiEqn.relax();
		res=psiEqn.solve().initialResidual();

		iter++;
        } 

}

// ************************************************************************* //
