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
    ParaSails

SourceFiles
    ParaSails.C

\*---------------------------------------------------------------------------*/

#ifndef ParaSails_H
#define ParaSails_H

#include "HypreIJPreconditioner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HypreIJPreconditioners
{ 

/*---------------------------------------------------------------------------*\
                           Class ParaSails Declaration
\*---------------------------------------------------------------------------*/

class ParaSails
:
    public HypreIJPreconditioner
{
    // Private data

      
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        ParaSails(const ParaSails&);

        //- Disallow default bitwise assignment
        void operator=(const ParaSails&);
        
protected:

       
 
public:

    //- Runtime type information
    TypeName("ParaSails"); 

    // Constructors

        //- Construct from components
        ParaSails
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~ParaSails()
        {}


    // Member Functions
    
        //- Correct
        virtual void initialize(HYPRE_Solver&) const;
                
        virtual HYPRE_PtrToSolverFcn setupPtr() const;               
        
        virtual HYPRE_PtrToSolverFcn solvePtr() const;
        
        virtual void destroy(HYPRE_Solver&) const;
        
        virtual const word precondName() const
        {
          return "ParaSails";
        }
        
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HypreIJPreconditioners 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
