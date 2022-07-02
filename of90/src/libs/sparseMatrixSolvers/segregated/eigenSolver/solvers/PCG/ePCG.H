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
    ConjugateGradient

SourceFiles
    ConjugateGradient.C

\*---------------------------------------------------------------------------*/

#ifndef ConjugateGradient_H
#define ConjugateGradient_H

#include "EigenIterDirSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Eigen::ConjugateGradient<
           Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>, 
           Eigen::Lower|Eigen::Upper,
           Eigen::IncompleteCholesky<double, Eigen::Lower, Eigen::AMDOrdering<int> >
        >  ICCPCG; 
        
typedef Eigen::ConjugateGradient<
           Eigen::SparseMatrix<double, EIGEN_STOR_ORDER>,
           Eigen::Lower|Eigen::Upper,
           Eigen::DiagonalPreconditioner<double> 
        >  DiagPCG;
                

namespace Foam
{
namespace EigenIterDirSolvers
{ 

/*---------------------------------------------------------------------------*\
                           Class ConjugateGradient Declaration
\*---------------------------------------------------------------------------*/

class ConjugateGradient
:
    public EigenIterDirSolver
{
    // Private data
        
        // Ptr to an ICC::CG precond-solver
        autoPtr<ICCPCG> solverICC_;
        
        // Ptr to a Diag::CG precond-solver
        autoPtr<DiagPCG> solverDiag_;
       
        // Name of preconditioner to be used
        word precondType_;
 
    // Private Member Functions

        //- Disallow default bitwise copy construct
        ConjugateGradient(const ConjugateGradient&);

        //- Disallow default bitwise assignment
        void operator=(const ConjugateGradient&);
        
protected:

       
 
public:

    //- Runtime type information
    TypeName("PCG");

    // Constructors

        //- Construct from components
        ConjugateGradient
        (
            const dictionary& dict
        );


    // Destructor

        virtual ~ConjugateGradient()
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
         if (precondType_ == "ICC")
          return "ConjugateGradient:ICC";
         else
          return "ConjugateGradient:Diag";
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
