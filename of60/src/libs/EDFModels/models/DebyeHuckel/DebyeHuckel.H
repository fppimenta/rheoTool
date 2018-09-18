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
    
Class
    DebyeHuckel

SourceFiles
    DebyeHuckel.C 

\*---------------------------------------------------------------------------*/

#ifndef DebyeHuckel_H
#define DebyeHuckel_H

#include "EDFEquation.H"

namespace Foam
{
namespace EDFEquations
{

class DebyeHuckel
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
            
            // Max number of iterations
            int maxIterPhiE_;
            int maxIterPsi_;
       
       class DHSpecie
          {   
              
            //- Model parameters
            
              // Ionic valence (-)     
                dimensionedScalar zi_;
              
              // Bulk concentration (mol/m3)   
                dimensionedScalar c0i_;
                
              // Species name     
                word namei_;  
          
           // Private Member Functions

              //- Disallow default bitwise copy construct
                DHSpecie(const DHSpecie&);

              //- Disallow default bitwise assignment
                void operator=(const DHSpecie&);
                
              public:
              
               // Constructor
               
                DHSpecie
                (
                    const word& name,
                    const surfaceScalarField& phi,
                    const dictionary& dict
                );
              
               // Destructor
                ~DHSpecie() {}
             
             // Public Member Functions 
             
             // Const-access to zi   
              const dimensionedScalar zi() const
               {
                   return zi_;
               } 
         
             // Const-access to c0i   
              const dimensionedScalar c0i() const
               {
                   return c0i_;
               }
               
             // Const access to namei   
              const word namei() const
               {
                   return namei_;
               } 
               
             // Const-access to member func used in 
             // computing the Debye length scale  
              const dimensionedScalar DebyeLengthP
               (
                 const dimensionedScalar relPerm,
                 const dimensionedScalar T
               ) const
               {
                   
                   return 
                   (    
                      zi_*zi_*c0i_         
                   );
                             
               }
          };
       
       //- List of species
       PtrList<DHSpecie> species_;
       
       //- Number of species
       int nSpecies_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        DebyeHuckel(const DebyeHuckel&);

        //- Disallow default bitwise assignment
        void operator=(const DebyeHuckel&);
       
public:

    //- Runtime type information
    TypeName("DebyeHuckel");

    // Constructors

        //- Construct from components
        DebyeHuckel
        (
            const word& name,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~DebyeHuckel()
        {}


    // Member Functions

        //- Return the electric force for the momentum equation
         virtual tmp<volVectorField> Fe() const;

        //- Update the electric-related component
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace EDFEquations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
