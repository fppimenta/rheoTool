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
    Foam::sparseSolver

Description
    Abstract base class for sparse solvers. Provides interfaces for external
    packages of linear solvers, as well as to the default OpenFOAM solvers. 

SourceFiles
    sparseSolver.C

\*---------------------------------------------------------------------------*/

#ifndef sparseSolver_H
#define sparseSolver_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "fvMatrix.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "cyclicFvPatch.H"
#include "sparseSolverBase.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sparseSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class sparseSolver :
   protected sparseSolverBase
{      
    // Private Member Functions

        //- Disallow copy construct
        sparseSolver(const sparseSolver&);

        //- Disallow default bitwise assignment
        void operator=(const sparseSolver&);
        
protected:

    // Protected data

        // Reference to mesh
        const fvMesh& mesh_;

        // Solution dictionary
        const dictionary& solDict_; 
        
        //- Add boundary contribution to source 
        void addBoundarySource
        (
           Field<Type>& source,
           fvMatrix<Type>& matrix,
           const GeometricField<Type, fvPatchField, volMesh>& T,
           const bool couples=true
        ) const;
        
        //- Add boundary contribution to matrix (diagonal) 
        void addBoundaryDiag
        (
           scalarField& diag,
           const direction cmpt,
           fvMatrix<Type>& matrix
        ) const;
        
        //- Compute residuals using OpenFOAM's definition (scaling). We also
        // use OpenFOAM machinery to perform the A*x product. 
        scalar getFoamResiduals
        (
           const GeometricField<Type, fvPatchField, volMesh>& T,
           fvMatrix<Type>& matrix,
           const scalarField& sourceCmpt,
           const scalarField& psiCmpt, 
           const scalarField& saveDiag,
           const FieldField<Field, scalar>& bouCoeffsCmpt,
           const lduInterfaceFieldPtrsList& interfaces,
           const int nEvalInit,
           const bool saveSystem,
           const direction cmpt,
           const int vcmpt        
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
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            sparseSolver,
            dictionary,
            (
                const GeometricField<Type, fvPatchField, volMesh>& T,
                const fvMesh& mesh,
                const dictionary& fvSolution
            ),
            (T, mesh, fvSolution)
        );


    // Constructors

        //- Construct from components
        sparseSolver
        (
            const GeometricField<Type, fvPatchField, volMesh>& T,
            const fvMesh& mesh,
            const dictionary& fvSolution
        )
        :
        mesh_(mesh),
        solDict_(fvSolution.subDict("solvers"))
        {}    
        
        
    // Selectors

        //- Return a pointer to a new sparseSolver created on freestore
        static autoPtr<sparseSolver> New
        (
            const GeometricField<Type, fvPatchField, volMesh>& T,
            const fvMesh& mesh,
            const dictionary& fvSolution
        );


    //- Destructor
    virtual ~sparseSolver();


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }
        
        //- Does the sparse solver relies on any external package? 
        virtual bool isExternalSolver() const = 0;
        
        //- Interface to applications
        virtual void solve
        (
            fvMatrix<Type>&
        ) = 0;
        
        //-- Reimplement upon need
        virtual void solve
        (
          fvMatrix<Type>&,
          const dictionary&
        );
        
        virtual void solve
        (
          const tmp<fvMatrix<Type>>&
        );
        
        //-- Reimplement upon need
        virtual void solve
        (
          const tmp<fvMatrix<Type>>&,
          const dictionary&
        );
            
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables


#define makeSparseTypeSolver(SS, Type)                                          \
  defineNamedTemplateTypeNameAndDebug(Foam::SS<Foam::Type>, 0);                \
    namespace Foam                                                             \
    {                                                                          \
                                                                               \
            sparseSolver<Type>::adddictionaryConstructorToTable<SS<Type>>       \
                add##SS##Type##ConstructorToTable_;                            \
                                                                               \
    }


#define makeSparseSolver(SS)                                                    \
                                                                               \
makeSparseTypeSolver(SS, scalar)                                                \
makeSparseTypeSolver(SS, vector)                                                \
makeSparseTypeSolver(SS, tensor)                                                \
makeSparseTypeSolver(SS, symmTensor)                                            \
makeSparseTypeSolver(SS, sphericalTensor)                                                

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
#ifdef NoRepository
    #include "sparseSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
