/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::eigenSolver
    
Description
    Interface to iterative/direct solvers of Eigen. 
 
\*---------------------------------------------------------------------------*/

#ifndef eigenSolver_H
#define eigenSolver_H

#include "sparseSolver.H"

#include "EigenIterDirSolver.H"
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Eigen::SparseMatrix<double, EIGEN_STOR_ORDER> spmat; 
typedef Eigen::Triplet<double> trip;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                       Class eigenSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class eigenSolver
:
    public sparseSolver<Type>
{
    // Private data
       
       // Create and save A, b, x and precond? Those will be re-used and kept
       // on the heap during all the simulation. There is of course memory 
       // overhead, but it allows reusing preconditioners.  
       Switch saveSystem_;
       
       // List of matrix of coefs (1 per valid component of Type)
       PtrList< spmat > listA_;
       
       // List of RHS (actually, there is only 1 element since it can be reused
       // for all components)
       PtrList< VectorXd > listb_;
       
       // List of solutions (actually, there is only 1 element since it can be reused
       // for all components)
       PtrList< VectorXd > listx_; 
       
       // List of sparse solvers (1 per valid component of Type)  
       List< autoPtr<EigenIterDirSolver> > listSolvers_;
       
       // Counter for the times solve() is called  
       int times_;
       
       // Frequency to update the preconditioner (based on times_ counter) 
       int updatePrecondFreq_;
       
       // Should the matrix of coeffs be updated all times solve() is called?
       bool updateA_;
       
       // Number of valid components of Type
       int nValidCmp_;
      
       // Sum of all coeffs of A to check for its variations over time
       // and warn the user in case they change and updateA_ is disabled 
       scalar sumofA_;
       
       // The time index upon call to ctor. It is needed to subtract it to the current
       // time index in order to get indices zeroed from call to ctor.
       label initTimeIndex;
       
       // Flag to detect auto mode to regulate preconditioner update. Only possible for
       // svaeSystem = true and updateA = true.
       bool autoPrecond;
       
       // CPU time for a solve call that includes updating the preconditioner. Only used
       // in auto mode.
       scalar solveTime1;
         
       // CPU time for a solve call immediately after updating the preconditioner. Only used
       // in auto mode.
       scalar solveTime2;
       
       // Time-step index before which matrices and solvers are always updated, in case
       // saveSystem_ is enabled
       static int nEvalInit_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        eigenSolver(const eigenSolver&);

        //- Disallow default bitwise assignment
        void operator=(const eigenSolver&);
        
        //- Check if the field being solved has boundary types that cannot be handled by the interface
        // Ideally this function should not exist.
        void checkLimitations
        (
          const GeometricField<Type, fvPatchField, volMesh>&
        );
        
        //- Solve the system reusing elements (A and preconditioner)
        void solveReuse
        (
          fvMatrix<Type>&
        );
        
        //- Solve the system creating all elements from scratch
        void solveNoReuse
        (
          fvMatrix<Type>&
        );
        
        //- Assemble A, b and x in eigen format. Values for A and b
        // are transfered from fvMatrix and x is set to 0.
        void assembleEigenAbx
        (
          spmat&,
          VectorXd&,
          VectorXd&,
          fvMatrix<Type>&,
          const GeometricField<Type, fvPatchField, volMesh>&,
          int
        );
        
        //- Assemble b and x in eigen format. b is built from 
        //  fvMatrix values, whereas x is set to 0.
        void assembleEigenBx
        (
          VectorXd&,
          VectorXd&,
          fvMatrix<Type>&,
          const GeometricField<Type, fvPatchField, volMesh>&,
          int
        );
        
        //- Transfers eigen solution vector (x) to the respective
        //  component of T
        void transferEigenSolution
        (
          VectorXd&,
          GeometricField<Type, fvPatchField, volMesh>&,
          int
        );
        
        //- Compute sum of A and check if it changes between calls to solve(). 
        // Check is performed once, controlled by nEvalInit_. This is how we
        // detect if the coeffs of A are changing and if updateA_ needs to be
        // enabled. Criteria not 100 % safe.
        virtual void checkMatrixSum
        (
          const scalarField& rowSum,
          const word name,
          const int tindex,
          const int vcmpt  
        );

public:

    //- Runtime type information
    TypeName("eigenSolver");


    // Constructors

        //- Construct from mesh
        eigenSolver
        (
            const GeometricField<Type, fvPatchField, volMesh>& T,
            const fvMesh& mesh,
            const dictionary& fvSolution
        );
        
    // Destructor

        virtual ~eigenSolver()
        {
           
        }
 
    // Member Functions
    
        // Does the sparse solver relies on any external package? 
        virtual bool isExternalSolver() const
        {
          return true;
        }
        
        // Interface to applications
        virtual void solve
        (
            fvMatrix<Type>&
        );
};

  
} // End namespace Foam

#ifdef NoRepository
    #include "eigenSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
