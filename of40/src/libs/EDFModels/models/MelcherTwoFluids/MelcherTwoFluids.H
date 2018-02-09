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

#ifndef MelcherTwoFluids_H
#define MelcherTwoFluids_H

#include "EDFEquation.H"

namespace Foam
{

class MTFPhase
  {   
      
    //- Model parameters
    
      // Electric permittivity of the fluid  
        dimensionedScalar relPermi_;
      
      // Fluid conductivity   
        dimensionedScalar sigmai_;
        
      // Name of the phase
        word phaseName_;  
  
   // Private Member Functions

      //- Disallow default bitwise copy construct
        MTFPhase(const MTFPhase&);

      //- Disallow default bitwise assignment
        void operator=(const MTFPhase&);
        
      public:
      
       // Constructor
       
        MTFPhase
        (
            const word& name,
            const surfaceScalarField& phi,
            const dictionary& dict
        );
      
       // Destructor
        virtual ~MTFPhase()
        {}
     
     // Public Member Functions 
     
     // Const-access to relPermi   
      const dimensionedScalar relPermi() const
       {
           return relPermi_;
       } 
 
     // Const-access to sigmai  
      const dimensionedScalar sigmai() const
       {
           return sigmai_;
       }
       
     // Const-access to the name of phase 
      const word phaseName() const
       {
           return phaseName_;
       }
  };

/*---------------------------------------------------------------------------*\
                           Class MelcherTwoFluids Declaration
\*---------------------------------------------------------------------------*/

class MelcherTwoFluids
:
    public EDFEquation
{
    // Private data

       //- Externally applied electric potential
       volScalarField rhoE_;

       //- Internal electric potential
       volScalarField psi_;
       
       //- Alpha-dependent permittivity
       volScalarField relPerm_;
       
       //- Alpha-dependent conductivity
       volScalarField sigma_;

       //- Extra electric field (optional)
       dimensionedVector extraE_;

       //- Solution controls for psi and rhoE
          
            // Residuals treshold
            scalar psiEqRes_;
            scalar rhoEEqRes_;
            
            // Max number of iterations
            scalar maxIterPsi_;
            scalar maxIterRhoE_;
            
       //- List of phases
       PtrList<MTFPhase> phases_;
       
       //- Number of phases (limited to 2 by the VOF component)
       int nPhases_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MelcherTwoFluids(const MelcherTwoFluids&);

        //- Disallow default bitwise assignment
        void operator=(const MelcherTwoFluids&);
     
public:

    //- Runtime type information
    TypeName("MelcherTwoFluids");

    // Constructors

        //- Construct from components
        MelcherTwoFluids
        (
            const word& name,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~MelcherTwoFluids()
        {}


    // Member Functions

        //- Return the electric force for the momentum equation
         virtual tmp<volVectorField> Fe() const;

        //- Update the electric force
        virtual void correct();
        
        //- Update the permittivity and conductivity with alpha
        void correctAlphaProperties();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
