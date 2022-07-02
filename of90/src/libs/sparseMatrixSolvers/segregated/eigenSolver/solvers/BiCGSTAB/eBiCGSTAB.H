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
    BiCGSTAB

SourceFiles
    BiCGSTAB.C

\*---------------------------------------------------------------------------*/

#ifndef BiCGSTAB_H
#define BiCGSTAB_H

#include "EigenIterDirSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Eigen::BiCGSTAB<
         Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>,
         Eigen::IncompleteLUT<double>
        >  ILUBiCGSTAB;

typedef Eigen::BiCGSTAB<
         Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>,
         Eigen::DiagonalPreconditioner<double> 
        >  DiagBiCGSTAB;

namespace Foam
{
namespace EigenIterDirSolvers
{ 

/*---------------------------------------------------------------------------*\
                           Class BiCGSTAB Declaration
\*---------------------------------------------------------------------------*/

class BiCGSTAB
:
    public EigenIterDirSolver
{
    // Private data

        // Ptr to an ILU::BiCGSTAB precond-solver
        autoPtr< ILUBiCGSTAB > solverILU_;
        
        // Ptr to a Diag::BiCGSTAB precond-solver
        autoPtr< DiagBiCGSTAB > solverDiag_;
       
        // Name of preconditioner to be used
        word precondType_;
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BiCGSTAB(const BiCGSTAB&);

        //- Disallow default bitwise assignment
        void operator=(const BiCGSTAB&);
        
protected:

       
 
public:

    //- Runtime type information
    TypeName("BiCGSTAB");

    // Constructors

        //- Construct from components
        BiCGSTAB
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~BiCGSTAB()
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
         if (precondType_ == "ILUT")
          return "BiCGSTAB:ILUT";
         else
          return "BiCGSTAB:Diag";
        }
        
        //- Is the solver iterative or direct ?
        virtual bool isIterative() const
        {
          return true;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace EigenIterDirSolvers 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
