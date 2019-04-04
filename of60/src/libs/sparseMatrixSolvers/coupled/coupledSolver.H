/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::coupledSolver
    
Description
    Solver for coupled systems using Petsc.  
 
\*---------------------------------------------------------------------------*/

#ifndef coupledSolver_H
#define coupledSolver_H

#include "fieldTypes.H"
#include "volFields.H"
#include "fvMatrix.H"

#include <petscksp.h>
 
#include "IPstream.H"
#include "OPstream.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "cyclicFvPatch.H"

#include "sparseSolverBase.H"
#include "LMatrix.H"
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
// For LMatrix
using namespace Foam::fvmb; 

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                           Class coupledSolver Declaration
\*---------------------------------------------------------------------------*/

class coupledSolver
:
   public IOList<label>,
   protected sparseSolverBase
{
    // Private data
        
        // Reference to mesh
        const fvMesh& mesh_;
        
        // List containing pointers to all scalar/vector/tensor/symmTensor/sphericalTensor
        // variables being solved for
        List<volScalarField*> varScalarList;
        List<volVectorField*> varVectorList;
        List<volTensorField*> varTensorList;
        List<volSymmTensorField*> varSymmTensorList;
        List<volSphericalTensorField*> varSphericalTensorList;
        
        // List with the location and type of the variables solved for. For example, if
        // p is the 2nd field in varScalarList and the 5th field inserted by insertField(),
        // then varLoc[5] = {1, ftscalar}.   
        List<labelList>  varLoc;
        
        // Concatenation of the variable names
        word    varNames;
        
        // List with names of the variables being solved. The order is the same as
        // in varLoc and corresponds to the order of calls to insertField().
        List<word>   varNamesList;
        
        // Same as varNamesList, but accounts for valid components name 
        List<word> varNamesCmpList;
        
        // varLoc only contains one entry per variable. Some of the variables are multicomponent
        // like vectors and one need to know the location of each of its components in the concatened
        // solution vector. firstCmpList contains the position of the first component of each field,
        // taking into acount only the solvable components (ie, in a 2D case eg, a vector has only 2 solvable
        // components). Example (3D): if varLoc = { {0, ftscalar}, {0, ftvector}, {0, fttensor}, {1, ftscalar}}, then
        // firstCmpList = {0, 1, 4, 13, 14}. Last element is the number of blocks in the matrix.
        // Order is the same as in varLoc.
        labelList    firstCmpList;
        
        // Local number of cells in process
        int nLocalCells;
        
        // Total number of cells (reduce(nLocalCells))  
        int nGlobalCells;
        
        // False until the first equation is inserted. After it becomes true,
        // no more fields can be inserted
        bool isSet;
        
        // Becomes true upon entering createSystem(). It is used to signal
        // that the A, b, x have been created and sized.
        bool isSysSized;
        
        // Becomes true upon entering createSystem(). Points the need of creating
        // sizing x. It can only be reverted to false if the mesh topo is changing.
        bool resetX;
        
        // Counter for the times solve() is called    
        int times_;
         
        // Field types
        enum {ftscalar, ftvector, ftsymmTensor, fttensor, ftsphericalTensor};
        
        // Petsc error        
        PetscErrorCode ierr;
        
        // Matrix of coefficients
        Mat A;
        
        // Matrix of coefficients in the previous time-step (only used to check
        // if the matrix is changing over time, in case saveSystem is enabled) 
        Mat A0;
        
        // Source vector
        Vec b;
        
        // Solution vector (all fields being solved are concatened sequentially in
        // in the order of varLoc)
        Vec x;
        
        // Solver context
        KSP ksp;
         
        // Sum of matrix coefficients (used to check if matrix is changing)
        PetscScalar Asum;
        
        // Should we save A, b and ksp for all the simulation time ? x is 
        // always saved in order to be possible to compute initial residuals.  
        bool saveSystem_;
       
        // Name of the coupled solver (the one specified in fvSolution)
        word name_;
       
       // Prefix for all Petsc Mat, Vec and KSP options. The prefix is the name of the
       // coupled solver. Used to distinguish between solvers in a unique database. 
       std::string prefix_;
       
       // Frequency to update the preconditioner (based on times_ counter) 
       int updatePrecondFreq_;
       
       // Should the matrix of coeffs be updated all times solve() is called?
       bool updateA_;
       
       // Enables a robust test to detected if A coeffs are changing between time-steps.
       bool isRobustSumCheck;
       
       // This flag is set to true once we test for the steadyness of the A coeffs. This
       // test is done once (the first time initTimeFlag becomes false) 
       bool sumCheckDone_;
       
       // This flag becomes false once the time index becomes equal to nEvalInit_. It has
       // been created because in a inner iteration loop time().timeIndex() is constant,
       // and does no allow to know the first of such iterations. On the other hand, this 
       // flag is switched the first time time().timeIndex() == nEvalInit_ is reached and
       // thus allows to detect first iteration of inner loop.
       // Only relevant if saveSystem_ = true && updateA_ = false .
       bool initTimeFlag;
       
       // The time index upon call to ctor. It is needed to subtract it to the current
       // time index in order to get indices zeroed from call to ctor.
       label initTimeIndex;
        
       // Time-step index before which matrices and solvers are always updated, in case
       // saveSystem_ is enabled
       static int nEvalInit_;
       
       // Flag to detect auto mode to regulate preconditioner update. Only possible for
       // svaeSystem = true and updateA = true.
       bool autoPrecond;
       
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
       
       // Array used for allocating A. Contains the number of expected off-diagonal entries (proc universe) 
       // in the matrix, per row. The array size is nLocalCells * nVarsToSolve. It is
       // empty for single proc run. Computed once for static mesh or repeatedly if topo changes.
       labelList maxOutProcBlocks_;
       
      
    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledSolver(const coupledSolver&);

        //- Disallow default bitwise assignment
        void operator=(const coupledSolver&);
        
        //- Create, size and allocate A, b, x and ksp
        void createSystem();
        
        //- Build the arrays used to pre-allocate A
        void computeAllocationPetsc
        (
          int nloc,
          int nglb
        );
  
        //- Assemble A and b in petsc format. Values for A and b
        //  are transfered from either a fvMatrix or a Lmatrix.
        template<class eqType>
        void assemblePetscAb
        (
          Mat& A,
          Vec& b,
          Vec& x,
          eqType& eqn,
          int cmpI,
          int rowBias,
          int colBias
        );
        
        //- Assemble b in petsc format. Values are transfered from
        //  either a fvMatrix or a Lmatrix.
        template<class eqType>
        void assemblePetscb
        (
          Vec& b,
          eqType& eqn,
          int cmpI,
          int rowBias
        );
        
        //- Transfers petsc solution vector (x) to the respective
        //  volField component. 
        template<class Type>
        void transferPetscSolution
        (
          Vec&,
          GeometricField<Type, fvPatchField, volMesh>&,
          int,
          int
        );
        
        //- This function is called for each field Type and grabs from
        //  x all the fields of that Type
        template<class Type>
        void getSolution();
        
        //- Computes residuals using openfoam's definition for segregated solvers.
        //  The formula needs improved because equations that are inserted may present 
        //  different scales, and less dominant ones will not express their residuals. 
        //  It also checks for A coeffs variation if needed.
        void getResiduals 
        ( 
          const Mat&,
          const Vec&,
          const Vec&,
          scalarList&
        );
        
        //- This is where Petsc solver is called
        void solvePetsc();
        
        //- Set all coeffs of A to zero
        void resetA();
        
        //- Set all elements of b to zero
        void resetb();
        
        //- Overloaded function returning a reference to the varTypeList
        //  which has the Type given as argument to the function. 
        inline List<Foam::volScalarField*>& varTypeList(scalar s)
        {
         return varScalarList; 
        }
        
        inline List<Foam::volVectorField*>& varTypeList(vector v)
        {
         return varVectorList; 
        }
        
        inline List<Foam::volSymmTensorField*>& varTypeList(symmTensor st)
        {
         return varSymmTensorList; 
        }
        
        inline List<Foam::volTensorField*>& varTypeList(tensor t)
        {
         return varTensorList; 
        }
        
        inline List<Foam::volSphericalTensorField*>& varTypeList(sphericalTensor spht)
        {
         return varSphericalTensorList; 
        }
        
        //- Overloaded function returning the enum Type corresponding
        //  to the Type of the argument passed to the function
        inline int ftType(scalar s) { return ftscalar; }
        inline int ftType(vector s) { return ftvector; }
        inline int ftType(symmTensor s) { return ftsymmTensor; }
        inline int ftType(tensor s) { return fttensor; }
        inline int ftType(sphericalTensor s) { return ftsphericalTensor; }
        
        //- Write A and b for debug 
        void printSystem(Mat& A, Vec& b);
       
public:

     //- Runtime type information
        TypeName("coupledSolver");
    
    // Constructors

        //- Construct given mesh
        coupledSolver
        (
          const word& name,
          const fvMesh& mesh
        );
        
    // Destructor
        virtual ~coupledSolver();


    // Member Functions
 
        // Insert fields
        // Note: insertField() should be called for all fields of the matrix
        // before insertEquation() is called for the first time.
        template<class Type>
        void insertField(GeometricField<Type, fvPatchField, volMesh>&);

        
        // Insert equations from fvMatrix
        // rowField is the name of variable for which we compute the equation
        // and colField is the name of the variable which is contributing to 
        // that equation.         
        template<class Type>
        void insertEquation
        (
          word rowField,
          word colField,
          fvMatrix<Type>& matrix
        );
        
        // Insert equations from tmp<fvMatrix>
        template<class Type>
        void insertEquation
        (
          word rowField,
          word colField,
          const tmp<fvMatrix<Type>>& matrix
        );
        
        // Insert equations from LMatrix
        template<class Type>
        void insertEquation
        (
          word rowField,
          word colField,
          LMatrix<Type>& matrix
        );
        
        // Insert equations from tmp<LMatrix>
        template<class Type>
        void insertEquation
        (
          word rowField,
          word colField,
          const tmp<LMatrix<Type>>& matrix
        );
        
        // Public interface to compute the solution
        void solve();
             
        // Access to mesh
        const fvMesh& mesh()
        {
          return mesh_;
        }
                 
};
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// All templated functions
#include "coupledSolverTemplates.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeTypeLMatrix(Type)                                                          \
template void Foam::coupledSolver::assemblePetscAb<Foam::fvmb::LMatrix<Foam::Type> >   \
(                                                                                      \
  Mat& A,                                                                              \
  Vec& b,                                                                              \
  Vec& x,                                                                              \
  Foam::fvmb::LMatrix<Foam::Type>& eqn,                                                \
  int cmpI,                                                                            \
  int rowBias,                                                                         \
  int colBias                                                                          \
);                                                                                     \
template void Foam::coupledSolver::assemblePetscb<Foam::fvmb::LMatrix<Foam::Type> >    \
(                                                                                      \
  Vec& b,                                                                              \
  Foam::fvmb::LMatrix<Foam::Type>& eqn,                                                \
  int cmpI,                                                                            \
  int rowBias                                                                          \
);                                                                                     \
                                                                                     
        
#define makeTypefvMatrix(Type)                                                         \
template void Foam::coupledSolver::assemblePetscAb<Foam::fvMatrix<Foam::Type> >        \
(                                                                                      \
  Mat& A,                                                                              \
  Vec& b,                                                                              \
  Vec& x,                                                                              \
  Foam::fvMatrix<Foam::Type>& eqn,                                                     \
  int cmpI,                                                                            \
  int rowBias,                                                                         \
  int colBias                                                                          \
);                                                                                     \
template void Foam::coupledSolver::assemblePetscb<Foam::fvMatrix<Foam::Type> >         \
(                                                                                      \
  Vec& b,                                                                              \
  Foam::fvMatrix<Foam::Type>& eqn,                                                     \
  int cmpI,                                                                            \
  int rowBias                                                                          \
);  

                        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
