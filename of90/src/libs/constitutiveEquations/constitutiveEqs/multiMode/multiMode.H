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
    multiMode

SourceFiles
    multiMode.C

\*---------------------------------------------------------------------------*/

#ifndef multiMode_H
#define multiMode_H

#include "constitutiveEq.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveEqs
{ 

/*---------------------------------------------------------------------------*\
                           Class multiMode Declaration
\*---------------------------------------------------------------------------*/

class multiMode
:
    public constitutiveEq
{
    // Private data

        //- Transported viscoelastic stress
        mutable volSymmTensorField tau_;

        //- List of models
        PtrList<constitutiveEq> models_;
        
        //- Solvent viscosity (read from dict; can not
        // use etaS() from models_).
        dimensionedScalar etaS_; 

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMode(const multiMode&);

        //- Disallow default bitwise assignment
        void operator=(const multiMode&);
        

protected:

    // There is no need to provide implementations for etaS,
    // etaP and coupledUtau() because functions divTau() and
    // divTauS() are directly implemented.
    
    // Note v3.0: etaS is needed for tauTotal()
    
    //- Return the solvent viscosity. In a multimode model,
    //  it is the sum of etaS from all modes.
    virtual const dimensionedScalar etaS() const
    {
      return etaS_;
    }
          
    // Same as etaS, but takes into account temperature effects.
    virtual tmp<volScalarField> etaSThermo() const;
    
    // Same as etaP, but takes into account temperature effects.
    // Note: to be used ONLY for stabilization purposes
    virtual tmp<volScalarField> etaPThermo() const;
    
    // Is the model prepared to work in non-isothermal conditions
    virtual bool hasThermo() const; 

public:

    //- Runtime type information
    TypeName("multiMode");

    // Constructors

        //- Construct from components
        multiMode
        (
            const word& name,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMode()
        {}


    // Member Functions

        //- Return the viscoelastic stress tensor
        virtual tmp<volSymmTensorField> tau() const;
        
        //- Return the density
        virtual const dimensionedScalar rho() const;

        //- Return the coupling term for the momentum equation
        virtual tmp<fvVectorMatrix> divTau(const volVectorField& U) const;
        
        //- Return the solvent extra-stress contribution + stabilizing terms for phase alpha
        virtual tmp<fvVectorMatrix> divTauS(const volVectorField& U, const volScalarField& alpha) const;
        
        //- Same as divTau but takes into account temperature effects (to be used in non-isothermal solvers)
        virtual tmp<fvVectorMatrix> divTauThermo(const volVectorField& U) const;
        
        //- Same as divTauS but takes into account temperature effects (to be used in non-isothermal solvers)
        virtual tmp<fvVectorMatrix> divTauSThermo(const volVectorField& U, const volScalarField& alpha) const;
        
        //- Return true if GNF (non-elastic)
        virtual bool isGNF() const
        {
          return false;
        };

        //- Correct the viscoelastic stress: alpha is the color function in two phase-flows
        virtual void correct
        (
          const volScalarField* alpha = nullptr,
          const volTensorField* gradU = nullptr
        );
        
        //- Adds implicitly the extra-stress polymer contribution in momentum equation 
        virtual void divTauImplCoupled() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constitutiveEqs 

/*---------------------------------------------------------------------------*\
                           Class constitutiveEqProtectedIntr Declaration
\*---------------------------------------------------------------------------*/

class constitutiveEqProtectedIntr
{
    
private:
  constitutiveEqProtectedIntr() {}
  
  // Only multiMode is allowed to create this class
  friend class constitutiveEqs::multiMode;

public:
  tmp<volScalarField> etaSThermo
  (
    const constitutiveEq* constitutiveEqPtr
  ) const
  {
    return constitutiveEqPtr->etaSThermo();
  }
  
  tmp<volScalarField> etaPThermo
  (
    const constitutiveEq* constitutiveEqPtr
  ) const
  {
    return constitutiveEqPtr->etaPThermo();
  }
  
  bool hasThermo
  (
    const constitutiveEq* constitutiveEqPtr
  ) const
  {
    return constitutiveEqPtr->hasThermo();
  }
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
