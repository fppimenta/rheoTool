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
    filmModel

Description
    Base class which can be directly called by applications to use any of 
    the models derived from class constitutivEq.

    Based on the libviscoelasticTransportModels library written by Jovani
    L. Favero. 

\*---------------------------------------------------------------------------*/

#ifndef filmModel_H
#define filmModel_H


#include "constitutiveEq.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class filmModel Declaration
\*---------------------------------------------------------------------------*/

class filmModel
:
    public IOdictionary
{
private:

    // Private data
    
        //- Reference to velocity field
        const volVectorField& U_;

        //- Non-isothermal?
        bool isThermo_;
        
        //- Heat transfer coefficient
        dimensionedScalar ht_;
        
        //- Heat capacity
        dimensionedScalar Cp_;
        
        //- Air temperature
        dimensionedScalar Tair_;

        //- Constitutive Eq
        autoPtr<constitutiveEq> eqPtr_;
        
        //- Temperature
        autoPtr<volScalarField> TPtr_;
        
        //- Height function
        volScalarField h_;
        
        //- Pressure field
        volScalarField p_;
        
        //- Do we need phi to be absolute before updating points on free-surface?
        bool absFluxNeeded_;
        
        //- Apply zero traction on outlet 
        bool isZeroTractionBC_;
        
        //- Enable postProcessing ?
        bool ppEnabled_;
        
        //- Dir name to output post-processing data
        fileName ppDir_;
        
        //- Sampling frequency (unit: nº of time-steps)
        label ppFreq_;
        
        //- Counter helper for post-processing
        label ppCnt_;
        
        //- File to where the film width at the outlet is written
        autoPtr<OFstream> outW_;
        
        //- File to where the film thickness along the outlet is written
        autoPtr<OFstream> outH_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        filmModel(const filmModel&);

        //- Disallow default bitwise assignment
        void operator=(const filmModel&);


public:

    //- Runtime type information
    TypeName("filmModel");


    // Constructors

        
        //- Construct from components
        filmModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

    // Destructor

        virtual ~filmModel()
        {}


    // Member Functions
    
        //- Retrieve flux balance accounting for thickness (except on free-surface)
        void checkFlux
        (
          const volVectorField& U,
          const surfaceScalarField& phiRel
        ) const;
    
        //- Update the free-surface boundary condition, so that the displacement gets defined
        // in the next time mesh.update() is called. 
        void updateFreeSurface
        (
          bool isFinalIter,
          const volVectorField& U,
          const surfaceScalarField& phiAbs
        ) const;
         
      
        //- Retrieve the thickness field
        const volScalarField& h () const
        {
          return h_; 
        }
        
        //- Retrieve the velocity field
        const volVectorField& U () const
        {
          return U_; 
        }
        
        //- Retrieve absFluxNeeded_
        const bool& absoluteFluxNeeded() const
        {
          return absFluxNeeded_;
        }
        
        //- Solve the equation for h
        void updateHeight(const surfaceScalarField& phi);
    
        //- Create the pressure field
        void updatePressure(const volTensorField& U);
        
        //- Computes dw/dz from du/dx and dv/dy, based on continuity 
        void correctGradU
        (
          volTensorField& gradU
        ) const;

        //- Return the viscoelastic stress tensor
        tmp<volSymmTensorField> tau() const;
        
        //- Return the stress tensor
        tmp<volSymmTensorField> tauTotal
        (
          const volVectorField& U
        ) const;
        
        //- Return the density
        const dimensionedScalar rho() const;

        //- Return the coupling term for the momentum equation
        tmp<fvVectorMatrix> divTau(const volVectorField& U);
         
        //- Is GNF or VE
        bool isGNF() const;
        
        //- Force the zero flux condition across all free-surface boundaries
        void forcedZeroFluxFreeSurface
        (
          bool force,
          surfaceScalarField& phiRel
        ) const;
       
        //- Correct the viscoelastic stress
        void correctStresses
        (
          const volVectorField& U
        );

        //- Correct temperature
        void correctThermo
        (
          const volVectorField& U,
          const surfaceScalarField& phi
        );
        
        //- Post-processing
        void postProcess
        (
          const volVectorField& U
        );

        //- Read constitutiveProperties dictionary
        virtual bool read();       
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
