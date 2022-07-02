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
    EigenIterDirSolver

\*---------------------------------------------------------------------------*/

#ifndef EigenIterDirSolver_H
#define EigenIterDirSolver_H


#define EIGEN_STOR_ORDER Eigen::RowMajor


#include "IOdictionary.H"
#include "Switch.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class EigenIterDirSolver Declaration
\*---------------------------------------------------------------------------*/

class EigenIterDirSolver
{
private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        EigenIterDirSolver(const EigenIterDirSolver&);

        //- Disallow default bitwise assignment
        void operator=(const EigenIterDirSolver&);
        
protected:

   // Protected Members
   
        // Solver dictionary
        dictionary dict_;
         
        // Number of iterations to convergence 
        int nIters_;
        
        // Final residual using eigen's definition
        scalar residual_;  

public:

    //- Runtime type information
    TypeName("EigenIterDirSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            EigenIterDirSolver,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        //- Construct from components
        EigenIterDirSolver
        (
            const dictionary& dict
        );


    // Selectors

        //- Return a reference  
        static autoPtr<EigenIterDirSolver> New
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~EigenIterDirSolver()
        {}


    // Member Functions

        //- Set user-defined options controlling the preconditioner
        virtual void initialize( ) = 0;
        
        //- Compute preconditioner
        virtual void setup
        (
          Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>& matrix,
          bool isSymmetric,
          bool reusePattern
        ) = 0;
               
        //- Solve 
        virtual void solve
        (
          const Eigen::VectorXd& source,
          Eigen::VectorXd& x 
        ) = 0;
        
        //- Return solver name        
        virtual word solverName() const = 0;
        
        //- Is the solver iterative or direct ?
        virtual bool isIterative() const = 0;
        
        //- Return number of iterations
        virtual int nIterations() const
        {
          return nIters_;
        }
        
        //- Return Eigen's residual
        virtual scalar residual() const
        {
          return residual_;
        }
        
        //- Print solver statistics to terminal 
        virtual void printInfo
        (
          word tName,
          scalar initResidual,
          scalar finalResidual
        ) const;     
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
