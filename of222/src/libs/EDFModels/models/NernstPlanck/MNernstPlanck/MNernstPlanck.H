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

#ifndef MNernstPlanck_H
#define MNernstPlanck_H

#include "EDFEquation.H"

namespace Foam
{

class MNPSpecie
  {   
      
    //- Model parameters
    
      // Molar concentration for specie i (mol/m3)
        volScalarField ci_;
    
      // Ionic valence (-)     
        dimensionedScalar zi_;
        
      // Diffusion coefficient (m2/s)     
        dimensionedScalar Di_;
        
     // Species name     
        word namei_;  
      
      // Effective size     
        dimensionedScalar ai_;
      
    // Private Member Functions

      //- Disallow default bitwise copy construct
        MNPSpecie(const MNPSpecie&);

      //- Disallow default bitwise assignment
        void operator=(const MNPSpecie&);
        
      public:
      
       // Constructor
       
        MNPSpecie
        (
            const word& name,
            const surfaceScalarField& phi,
            const dictionary& dict
        );
      
       // Destructor
        virtual ~MNPSpecie()
        {}
     
     // Public Member Functions 
     
     // Const access to zi   
      const dimensionedScalar zi() const
       {
           return zi_;
       } 
       
     // Const access to Di   
      const dimensionedScalar Di() const
       {
           return Di_;
       } 
       
     // Const access to namei   
      const word namei() const
       {
           return namei_;
       } 
       
     // Const access to a   
      const dimensionedScalar ai() const
       {
           return ai_;
       } 
 
     // Const access to ci   
      const volScalarField& ci() const
       {
           return ci_;
       }
       
      // Non-const access to ci   
      volScalarField& ci()
       {
           return ci_;
       }
  };

/*---------------------------------------------------------------------------*\
                           Class MNernstPlanck Declaration
\*---------------------------------------------------------------------------*/

class MNernstPlanck
:
    public EDFEquation
{
    // Private data
    
        //- True if phiE_ exists in the directory, false otherwise  
        bool solvePhiE_;
        
        //- Intrinsic electric potential
        volScalarField psi_;

        //- Externally applied electric potential
        volScalarField phiE_;        

        //- Model constants

            // Electrical permittivity of the fluid
            dimensionedScalar relPerm_;
            
            // Absolute temperature (K)  
            dimensionedScalar T_;

       //- Extra electric field (optional)
       dimensionedVector extraE_;
       
       //- Should the electric field from the intrinsic potential enter f_E
       bool psiContrib_;

       //- Solution controls for phiE and psi
          
            // Residuals treshold
            scalar phiEEqRes_;
            scalar psiEqRes_;
            scalar ciEqRes_;
            
            // Max number of iterations
            int maxIterPhiE_;
            int maxIterPsi_;
            int maxIterCi_;
            
       //- List of species
       PtrList<MNPSpecie> species_;
       
       //- Number of species
       int nSpecies_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        MNernstPlanck(const MNernstPlanck&);

        //- Disallow default bitwise assignment
        void operator=(const MNernstPlanck&);
         
public:

    //- Runtime type information
    TypeName("MNernstPlanck");

    // Constructors

        //- Construct from components
        MNernstPlanck
        (
            const word& name,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~MNernstPlanck()
        {}


    // Member Functions

        //- Return the electric force for the momentum equation
         virtual tmp<volVectorField> Fe() const;

        //- Update the electric-related component
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
