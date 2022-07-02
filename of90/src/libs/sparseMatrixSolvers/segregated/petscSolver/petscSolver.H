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
    Foam::petscSolver 
    
Description
    Interface to the KSP context of Petsc.   
 

\*---------------------------------------------------------------------------*/

#ifndef petscSolver_H
#define petscSolver_H

#include "sparseSolver.H"
#include <petscksp.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                       Class petscSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class petscSolver
:
    public sparseSolver<Type>
{
    // Private data
       
       // Create and save A, b, x and precond? Those will be re-used and kept
       // on the heap during all the simulation. There is of course memory 
       // overhead, but it allows reusing preconditioners.  
       Switch saveSystem_;
       
       // List of matrix of coefs (1 per valid component of Type)
       PtrList< Mat > listA_; 
       
       // List of RHS (actually, there is only 1 element since it can be reused
       // for all components)
       PtrList< Vec > listb_;
       
       // List of solutions (actually, there is only 1 element since it can be reused
       // for all components)
       PtrList< Vec > listx_;
       
       // List of linear solvers context (1 per valid component of Type)  
       PtrList<KSP> listSolvers_;
       
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
       
       // Time-step index before which matrices and solvers are always updated, in case
       // saveSystem_ is enabled
       static int nEvalInit_;
       
       // Petsc error variable 
       PetscErrorCode ierr;
       
       // Prefix for all Petsc Mat, Vec and KSP options. The prefix is the name of the
       // variable being solved. Used to distinguish between vars in a unique database. 
       std::string prefix_;
       
       // The time index upon call to ctor. It is needed to subtract it to the current
       // time index in order to get indices zeroed from call to ctor.
       label initTimeIndex;
       
       // Flag to detect auto mode to regulate preconditioner update. Only possible for
       // svaeSystem = true and updateA = true.
       bool autoPrecond;
       
       // True if there are cyclicAMI patches in the mesh. 
       bool isThereCyclicAMI_;
       
       // The ID corresponding to the given mesh region
       label meshID_;
       
       // CPU time for a solve call that includes updating the preconditioner. Only used
       // in auto mode.
       scalar solveTime1;
         
       // CPU time for a solve call immediately after updating the preconditioner. Only used
       // in auto mode.
       scalar solveTime2;
        
       // Array used for allocating A. Contains the number of expected diagonal entries (proc universe) 
       // in the matrix, per row. The array size is nLocalCells * nVarsToSolve. Computed
       // once for static mesh or repeatedly if topo changes.
       labelList maxInProcBlocks_;
       
       // Array used for allocating of A. Contains the number of expected diagonal entries (proc universe) 
       // in the matrix, per row. The array size is the local size of T (number of local rows). Computed
       // once for static mesh or repeatedly if topo changes.
       labelList maxInProcFaces_;
       
       // Array used for allocating of A. Contains the number of expected off-diagonal entries (proc universe) 
       // in the matrix, per row. The array size is the local size of T (number of local rows). It is
       // empty for single proc run. Computed once for static mesh or repeatedly if topo changes.
       labelList maxOutProcFaces_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        petscSolver(const petscSolver&);

        //- Disallow default bitwise assignment
        void operator=(const petscSolver&);
        
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
        
        //- Build the arrays used to pre-allocate A
        void computeAllocationPetsc
        (
          fvMatrix<Type>&,
          const GeometricField<Type, fvPatchField, volMesh>&
        );
        
        //- Assemble A, b and x in petsc format. Values for A and b
        // are transfered from fvMatrix and x is set to 0.
        void assemblePetscAbx
        (
          Mat&,
          Vec&,
          Vec&,
          fvMatrix<Type>&,
          const GeometricField<Type, fvPatchField, volMesh>&,
          int,
          int
        );  
        
        //- Assemble b and x in petsc format. b is built from 
        //  fvMatrix values, whereas x is set to 0.
        void assemblePetscbx
        (
          Vec&,
          Vec&,
          fvMatrix<Type>&,
          const GeometricField<Type, fvPatchField, volMesh>&,
          int,
          int
        );  
        
        //- Transfers petsc solution vector (x) to the respective
        //  component of T
        void transferPetscSolution
        (
          Vec&,
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
    TypeName("petscSolver");


    // Constructors

        //- Construct from mesh
        petscSolver
        (
            const GeometricField<Type, fvPatchField, volMesh>& T,
            const fvMesh& mesh,
            const dictionary& fvSolution
        );
        
    // Destructor

        virtual ~petscSolver()
        {
          // Need to destroy explicitly the variables allocated
          // dynamically if saveSystem_ was enabled   
          if (saveSystem_)
          {
             forAll(listSolvers_, i)
             {
               KSPDestroy(&listSolvers_[i]);             
             }
             
             VecDestroy(&listb_[0]);
             VecDestroy(&listx_[0]);
             
             forAll(listA_, i)
             {
               MatDestroy(&listA_[i]); 
             }
          }
          
          this->counterPetsc_--;
          
          // Finalize petsc for the last object of this type
          if (this->counterPetsc_==0)
           ierr = PetscFinalize();
          
          // Since we start MPI, we are also responsible to finalize it 
          if (this->counterHypre_ + this->counterPetsc_ == 0 && !Pstream::parRun())
            MPI_Finalize();       
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
    #include "petscSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
