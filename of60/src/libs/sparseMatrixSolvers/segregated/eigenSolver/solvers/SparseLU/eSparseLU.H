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
    SparseLU

SourceFiles
    SparseLU.C

\*---------------------------------------------------------------------------*/

#ifndef SparseLU_H
#define SparseLU_H

#include "EigenIterDirSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Eigen::SparseLU<
         Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>,   
         Eigen::COLAMDOrdering<int>
        >  sparseLU;
         
namespace Foam
{
namespace EigenIterDirSolvers
{ 

/*---------------------------------------------------------------------------*\
                           Class SparseLU Declaration
\*---------------------------------------------------------------------------*/

class SparseLU
:
    public EigenIterDirSolver
{
    // Private data

        // Ptr to a sparseLU solver
        autoPtr< sparseLU > sparseLU_;
          
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        SparseLU(const SparseLU&);

        //- Disallow default bitwise assignment
        void operator=(const SparseLU&);
        
protected:

       
 
public:

    //- Runtime type information
    TypeName("SparseLU");

    // Constructors

        //- Construct from components
        SparseLU
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~SparseLU()
        {}


    // Member Functions
    
        //- Set user-defined options controlling the preconditioner
        virtual void initialize();
        
        //- Compute preconditioner
        virtual void setup
        (
           Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>& matrix,
           bool isSymmetric,
           bool reusePattern
        );               
        
        //- Solve
        virtual void solve
        (
           const Eigen::VectorXd& source,
           Eigen::VectorXd& x 
        );
        
        //- Return solver name   
        virtual word solverName() const
        {          
          return "SparseLU";
        }
        
        //- Is the solver iterative or direct ?
        virtual bool isIterative() const
        {
          return false;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace EigenIterDirSolvers 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
