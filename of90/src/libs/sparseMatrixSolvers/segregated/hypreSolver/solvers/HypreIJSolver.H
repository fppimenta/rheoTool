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
    HypreIJSolver
    
Decription
    Base class. This class only implements methods, which it applies to
    the objects received.

\*---------------------------------------------------------------------------*/

#ifndef HypreIJSolver_H
#define HypreIJSolver_H


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
                        Class HypreIJSolver Declaration
\*---------------------------------------------------------------------------*/

class HypreIJSolver
{
private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        HypreIJSolver(const HypreIJSolver&);

        //- Disallow default bitwise assignment
        void operator=(const HypreIJSolver&);
        
protected:

   // Protected Members
   
        // Solver dictionary
        dictionary dict_;
         
        // Number of iterations to convergence 
        int nIters_;
        
        // Final residual using hypre's definition
        scalar residual_;   

public:

    //- Runtime type information
    TypeName("HypreIJSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            HypreIJSolver,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from components
        HypreIJSolver
        (
            const dictionary& dict
        );


    // Selectors

        //- Return a reference  
        static autoPtr<HypreIJSolver> New
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~HypreIJSolver()
        {}


    // Member Functions

        //- Set user-defined options controlling the solver 
        virtual void initialize(HYPRE_Solver&) const = 0;
        
        //- Setup and compute preconditioner if it applies
        virtual void setup
        (
          HYPRE_Solver&, 
          HYPRE_Solver&,
          HYPRE_ParCSRMatrix&,
          HYPRE_ParVector&,
          HYPRE_ParVector&,
          HYPRE_PtrToSolverFcn,
          HYPRE_PtrToSolverFcn,
          bool
        ) const = 0;
               
        //- Solve
        virtual void solve
        (
          HYPRE_Solver&, 
          HYPRE_ParCSRMatrix&,
          HYPRE_ParVector&,
          HYPRE_ParVector&
        ) = 0;
        
        //- Destroy
        virtual void destroy(HYPRE_Solver&) const = 0;
        
        //- Return the name of the solver
        virtual const word solverName() const = 0;
        
        //- Return the number of iterations
        virtual int nIterations() const
        {
          return nIters_;
        }
        
        //- Return hypre's residual
        virtual scalar residual() const
        {
          return residual_;
        }
        
        //- Print solver statistics to terminal 
        virtual void printInfo
        (
          word tName,
          word pcName,
          scalar initResidual,
          scalar finalResidual
        ) const;     
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
