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
#include "cyclicAMIFvPatch.H"   
#include "cyclicFvPatchField.H" 

#include "regionCoupledAMIFvPatch.H" 
#include "regionCoupledBaseFvPatch.H"
#include "coupledTFvPatchScalarField.H"
 
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
    
        // Struct with mesh-related info
        struct meshListS
        {
          const fvMesh* mesh;      // ptr to the mesh
          label ID;                // corresponding mesh ID in sharedData
          label nValidCmp;         // valid components solved for this mesh 
  
          meshListS()
          :
          mesh(NULL),
          ID(-1),
          nValidCmp(0)
          { }
          
        };
        
        // Struct with var-related info
        struct varInfoS
        {
          // Data
          label localID;         // index of the field in var<Type>List
          label typeV;           // label corresponding to the type of the field
          label nValidCmp;       // number of solvable components for that field 
          label nCells;          // global size of the field (reduceOp) 
          label meshID;          // index of the corresponding mesh in meshList
          word  name;            // name of the field 
          label firstCmp;        // the position (index) of the first component of the field in the solution vector x
          label firstElem;       // the global index of the first element of the first solvable component of the field in vector x        
        };
        
        // List with mesh-related info. One entry per mesh.    
        List<meshListS> meshList;
        
        // List containing pointers to all scalar/vector/tensor/symmTensor/sphericalTensor
        // variables being solved for
        List<volScalarField*> varScalarList;
        List<volVectorField*> varVectorList;
        List<volTensorField*> varTensorList;
        List<volSymmTensorField*> varSymmTensorList;
        List<volSphericalTensorField*> varSphericalTensorList;
   
        // List with field-related info. One entry per field.     
        List<varInfoS> varInfo;
        
        // Concatenation of the variable names
        word    varNames;
         
        // Same as varInfoS.name, but accounts for valid components name 
        List<word> varNamesCmpList;
               
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
        autoPtr<Mat> A0;
        
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
       
        // Name of the coupled solver (the one specified in fvSolution, different
        // from the one registered in the database, which contains reference to 
        // the mesh).
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
       
       // True if there are cyclicAMI patches in the mesh. 
       bool isThereCyclicAMI_;
       
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
        
        //- Check if the field being solved has boundary types that cannot be handled by the interface
        // Ideally this function should not exist.
        void checkLimitations
        (
          const fvMesh&
        );
        
        //- Create, size and allocate A, b, x and ksp
        void createSystem();
        
        //- Build the arrays used to pre-allocate A
        void computeAllocationPetsc
        (
          int nloc,
          int nglb
        );
       
        // Return the position in varInfo of the field having
        // name == fName
        label findField
        (
          const word& fName,
          const word& meshName  
        ) const
        {
          label id(-1);
          forAll(varInfo, i)
          {
            // Check both field name and mesh name. This allows having fields with same name
            // but owned by different meshes.
            if (varInfo[i].name == fName && meshList[varInfo[i].meshID].mesh->name() == meshName)
             id = i;
          }
          
          return id;
        }
  
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
          int colBias,
          int rowVarID   
        );
        
        //- Assemble b in petsc format. Values are transfered from
        //  either a fvMatrix or a Lmatrix.
        template<class eqType>
        void assemblePetscb
        (
          Vec& b,
          eqType& eqn,
          int cmpI,
          int rowBias,
          int rowVarID   
        );
        
        //- Transfers petsc solution vector (x) to the respective
        //  volField component. 
        template<class Type>
        void transferPetscSolution
        (
          Vec&,
          GeometricField<Type, fvPatchField, volMesh>&,
          int,
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

        //- Construct given time
        coupledSolver
        (
           const word& solverName,
           const word& meshName,
           const Time& time,
           const dictionary& dict
        );
        
    // Destructor
        virtual ~coupledSolver();


    // Member Functions
    
    
        // Insert mesh
        // Note: all meshes of a coupled system must be inserted before the fields
        void insertMesh(const fvMesh&);
 
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
        
        // Interface to compute the solution
        void solve();                 
};
 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// All templated functions
#ifdef NoRepository
    #include "coupledSolverTemplates.C"
#endif
 
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
  int colBias,                                                                         \
  int rowVarID                                                                         \
);                                                                                     \
template void Foam::coupledSolver::assemblePetscb<Foam::fvmb::LMatrix<Foam::Type> >    \
(                                                                                      \
  Vec& b,                                                                              \
  Foam::fvmb::LMatrix<Foam::Type>& eqn,                                                \
  int cmpI,                                                                            \
  int rowBias,                                                                         \
  int rowVarID                                                                         \
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
  int colBias,                                                                         \
  int rowVarID                                                                         \
);                                                                                     \
template void Foam::coupledSolver::assemblePetscb<Foam::fvMatrix<Foam::Type> >         \
(                                                                                      \
  Vec& b,                                                                              \
  Foam::fvMatrix<Foam::Type>& eqn,                                                     \
  int cmpI,                                                                            \
  int rowBias,                                                                         \
  int rowVarID                                                                         \
);  

                        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
