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
    HypreIJPreconditioner
    
Decription
    Base class. This class only implements methods, which it applies to
    the objects received.

\*---------------------------------------------------------------------------*/

#ifndef HypreIJPreconditioner_H
#define HypreIJPreconditioner_H


#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HypreIJPreconditioner Declaration
\*---------------------------------------------------------------------------*/

class HypreIJPreconditioner
{
private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        HypreIJPreconditioner(const HypreIJPreconditioner&);

        //- Disallow default bitwise assignment
        void operator=(const HypreIJPreconditioner&);
        
protected:

   // Protected Members
   
     // Solver dictionary
     dictionary dict_;
      

public:

    //- Runtime type information
    TypeName("HypreIJPreconditioner");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            HypreIJPreconditioner,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from components
        HypreIJPreconditioner
        (
            const dictionary& dict
        );


    // Selectors

        //- Return a reference  
        static autoPtr<HypreIJPreconditioner> New
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~HypreIJPreconditioner()
        {}


    // Member Functions

        //- Set user-defined options controlling the preconditioner
        virtual void initialize(HYPRE_Solver&) const = 0;
        
        //- Return precond functor to setup
        virtual HYPRE_PtrToSolverFcn setupPtr() const = 0;               
        
        //- Return precond functor to solve
        virtual HYPRE_PtrToSolverFcn solvePtr() const = 0;
        
        //- Destroy
        virtual void destroy(HYPRE_Solver&) const = 0;
        
        //- Return the name of the precond
        virtual const word precondName() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
