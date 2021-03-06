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
    BoomerAMG

SourceFiles
    BoomerAMG.C

\*---------------------------------------------------------------------------*/

#ifndef BoomerAMG_H
#define BoomerAMG_H

#include "HypreIJPreconditioner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJPreconditioners
{ 

/*---------------------------------------------------------------------------*\
                           Class BoomerAMG Declaration
\*---------------------------------------------------------------------------*/

class BoomerAMG
:
    public HypreIJPreconditioner
{
    // Private data

      
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BoomerAMG(const BoomerAMG&);

        //- Disallow default bitwise assignment
        void operator=(const BoomerAMG&);
        
protected:

       
 
public:

    //- Runtime type information
    TypeName("BoomerAMG");

    // Constructors

        //- Construct from components
        BoomerAMG
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~BoomerAMG()
        {}


    // Member Functions
    
        //- Set user-defined options controlling the preconditioner
        virtual void initialize(HYPRE_Solver&) const;
        
        //- Return precond functor to setup
        virtual HYPRE_PtrToSolverFcn setupPtr() const;               
        
        //- Return precond functor to solve
        virtual HYPRE_PtrToSolverFcn solvePtr() const;
        
        //- Destroy
        virtual void destroy(HYPRE_Solver&) const;
        
        //- Return the name of the precond
        virtual const word precondName() const
        {
          return "AMG";
        }
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HypreIJPreconditioners 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
