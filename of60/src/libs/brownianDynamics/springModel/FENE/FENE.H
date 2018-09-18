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
    FENE
    
Description
     
    This class is part of rheoTool.

SourceFiles
    FENE.C

\*---------------------------------------------------------------------------*/

#ifndef FENE_H
#define FENE_H

#include "springModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace springModels
{

/*---------------------------------------------------------------------------*\
                           Class FENE Declaration
\*---------------------------------------------------------------------------*/

class FENE
:
    public springModel
{
    // Private data
    
         //- Maximum fractional lenght of a spring
         scalar MAXFRACL_;
     
    // Private Member Functions

        //- Disallow default bitwise copy construct
        FENE(const FENE&);

        //- Disallow default bitwise assignment
        void operator=(const FENE&);

protected:

 // Protected member functions  
    
    //- Spring force (call per molecule)
    virtual tmp<vectorField> fSpringI(vectorField&, label, bool, bool);
        
    //- Function f of Newton-Raphson method
    virtual void fSIM
    (
      label mi,
      label cmpi, 
      const Field<scalar>& xStar, 
      const Field<vector>& x,
      Field<scalar>& fm
    );

    //- Jacobian of f in Newton-Raphson method
    virtual void jacobianSIM
    (
      label mi,
      label cmpi,
      const Field<vector>& x,
      scalarSquareMatrix& J
    ); 

public:


    //- Runtime type information
    TypeName("FENE");


    // Constructors

        //- Construct from components
        FENE
        (
            const dictionary& dict,
            const volVectorField& U,
            sPCloudInterface& sPCI
        );

    // Destructor

        virtual ~FENE()
        {}


    // Member Functions
   	
	//- Return true if the model limits the length of springs
	virtual bool isLimitedSpring()
	{
	   return true;
	}        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace springModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
