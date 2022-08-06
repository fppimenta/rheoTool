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
    constitutiveEq

Description
    Abstract base class from which each specific GNF or viscoelastic model
    is derived.

    Based on the libviscoelasticTransportModels library written by Jovani
    L. Favero. 

\*---------------------------------------------------------------------------*/

#ifndef constitutiveEq_H
#define constitutiveEq_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "NamedEnum.H"
#include "runTimeSelectionTables.H"
#include "extrapolatedCalculatedFvPatchField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class constitutiveEq Declaration
\*---------------------------------------------------------------------------*/

class constitutiveEq
{
private:

    // Private data

        //- Name
        const word name_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Reference to face flux field
        const surfaceScalarField& phi_;
        

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constitutiveEq(const constitutiveEq&);

        //- Disallow default bitwise assignment
        void operator=(const constitutiveEq&);
        
        //- An intruder class to access protected members from derived classes
        friend class constitutiveEqProtectedIntr;
        
        //- Class for film casting
        friend class filmModel;
        
public:

    //- Stabilization options
    enum stabOptions
    {
      soNone,
      soBSD,
      soCoupling
    };
        
protected:

   // Protected Members
   
        //- Return B and Sigma matrices after decomposition of grad(U).T() for the LCM
        void decomposeGradU
        (
          const volTensorField& M,
          const volTensorField& eigVals, 
          const volTensorField& eigVecs,
          volTensorField& sigma, 
          volTensorField& B
        );

        //- Return eigenvectors/values of theta for the LCM
        void calcEig(const volSymmTensorField& theta, volTensorField& vals, volTensorField& vecs);
        
        //- Return the strain rate magnitude for GNF models
        inline tmp<volScalarField> strainRate(const volTensorField* gradU = nullptr)
        {
          if (gradU == nullptr)
          {
            return ( sqrt(2.0)*mag( symm(fvc::grad(U_)) ) );
          }
          else
          {
            return ( sqrt(2.0)*mag( symm(*gradU) ) );
          }
        }
        
        //- Check which stabilization method to use
        void checkForStab(const dictionary& dict);
        
        //- Check is solver is coupled
        void checkIfCoupledSolver
        (
          const dictionary& dict,
          volSymmTensorField& tau
        );
      
        
  // The following protected functions are intended to be used inside divTau() and divTauS() functions,
  // which should be general for almost all models (exception: multiMode). Any other use should be
  // carefully analyzed, with special attention to the mentioned exceptions.
        
      //- Return the viscosity field for GNF models (default behavior is a null field - for VE models)  
      virtual const volScalarField& eta() const
      {      
         return ((sqrt(2.0)*mag( symm(fvc::grad(U_)) ))*dimensionedScalar("zeroU", dimensionSet(1, -1, 0, 0, 0, 0, 0), 0));       
      }
       
      //- Return the solvent viscosity for VE models (default behavior is a null viscosity - for GNF models)  
      virtual const dimensionedScalar etaS() const
      {
         return (dimensionedScalar("zeroU", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0));
      }
       
      //- Return the polymeric viscosity for VE models (default behavior is a null viscosity - for GNF models)
      // Does not consider possible variation of etaP with shear rate (there is no need for it where is used)  
      virtual const dimensionedScalar etaP() const
      {
         return (dimensionedScalar("zeroU", dimensionSet(1, -1, -1, 0, 0, 0, 0), 0));
      }
      
      // Same as etaS, but takes into account temperature effects.
      // To be redefined in non-GNF non-isothermal models. When called for
      // models not thermo-ready (hasThermo == false) it returns etaS. Tmp fix (>v5)
      // to get tauTotal() working in rheoTestFoam with both thermo and non-thermo
      // models.
      virtual tmp<volScalarField> etaSThermo() const
      {         
        tmp<volScalarField> eta
        (
          new volScalarField
          (
            IOobject
            (
              "etaSThermo",
              U().instance(),
              U().db()
            ),
            U().mesh(),
            etaS(),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
          )
        );
        return eta;
      }
      
      // Same as etaP, but takes into account temperature effects.
      // To be redefined in non-GNF non-isothermal models
      virtual tmp<volScalarField> etaPThermo() const
      {
         
        tmp<volScalarField> eta
        (
          new volScalarField
          (
            IOobject
            (
              "etaPThermo",
              U().instance(),
              U().db()
            ),
            U().mesh(),
            etaP(),
            extrapolatedCalculatedFvPatchField<scalar>::typeName
          )
        );        
        return eta;
      }
      
      // Is the model prepared to work in non-isothermal conditions
      virtual bool hasThermo() const
      {
         return false;
      }
       
      static const NamedEnum<stabOptions, 3> stabOptionNames_;
      stabOptions stabOption_;
      
      //- Should the constitutive equation be solved coupled with momentum
      // and continuity?
      bool solveCoupled_;
      
      // Hot fix functions for bug in gcc > 9.0
      tmp<volTensorField> innerP
      (
        const volTensorField& t1,
        const volTensorField& t2,
        bool isFirstT
      ) const;
      
      tmp<volTensorField> innerP
      (
        const volTensorField& t1,
        const tmp<volTensorField>& t2,
        bool isFirstT
      ) const;

public:

    //- Runtime type information
    TypeName("constitutiveEq");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            constitutiveEq,
            dictionary,
            (
                const word& name,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const dictionary& dict
            ),
            (name, U, phi, dict)
        );


    // Constructors

        //- Construct from components
        constitutiveEq
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Selectors

        //- Return a reference to the selected viscoelastic law
        static autoPtr<constitutiveEq> New
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~constitutiveEq()
        {}

    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return the velocity field
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the face flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }
                
        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const = 0;
        
        //- Return the total stress tensor (solvent + polymeric)
        //  using alpha1 as weight
        virtual tmp<volSymmTensorField> tauTotal() const;
        
        //- Return the density
        virtual const dimensionedScalar rho() const = 0;
         
        //- Return the coupling term for the momentum equation (since it is virtual non-pure,
        // it does not need to be redefined in derived classes, although it can if any model needs
        // a different implementation)
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U) const;
        
        //- Return the solvent extra-stress contribution + stabilizing terms for phase alpha
        // (since it is virtual non-pure, it does not need to be redefined in derived classes, 
        // although it can if any model needs a different implementation)
        virtual tmp<fvVectorMatrix> divTauS(const volVectorField& U, const volScalarField& alpha) const;
        
        //- Same as divTau but takes into account temperature effects (to be used in non-isothermal solvers)
        virtual tmp<fvVectorMatrix> divTauThermo(const volVectorField& U) const;
        
        //- Same as divTauS but takes into account temperature effects (to be used in non-isothermal solvers)
        virtual tmp<fvVectorMatrix> divTauSThermo(const volVectorField& U, const volScalarField& alpha) const;
        
        //- Adds implicitly the extra-stress polymer contribution in momentum equation 
        virtual void divTauImplCoupled() const;
        
        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const = 0;

        //- Correct the viscoelastic stress: alpha is the color function in two phase-flows
        virtual void correct
        (
          const volScalarField* alpha = nullptr,
          const volTensorField* gradU = nullptr
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
