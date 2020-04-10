template<class Type>
void Foam::coupledSolver::insertField  
(
  GeometricField<Type, fvPatchField, volMesh>& field
)
{
  if (isSet)
  {
    FatalErrorIn
    (
        "Foam::coupledSolver::addField"
    )   << nl << nl << "Cannot insert fields after inserting equation. "
        << exit(FatalError);
  }

  typename pTraits<Type>::labelType validComponents
  (
    field.mesh().template validComponents<Type>()
  );
 
  int i = 0;
  for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
  {
    if (component(validComponents, cmpt) == -1) continue; 
    
    varNamesCmpList.append(word(field.name()+pTraits<Type>::componentNames[cmpt]));
    i++;
  }
  
  // Element of type Type
  Type type(pTraits<Type>::zero); 
 
  varNames += field.name() + ", ";
   
  varTypeList(type).append(&field);
  varLoc.append({varTypeList(type).size()-1, ftType(type)});
  varNamesList.append(field.name());
  
  // Set the position for the first component of the NEXT variable to be added.
  // A Type variable may add up to N blocks because it has
  // up to N solvable components
  firstCmpList.append(firstCmpList.last() + i);
}

template<class Type>
void Foam::coupledSolver::insertEquation
(
  word rowField,
  word colField,
  fvMatrix<Type>& matrix
)
{
 createSystem();
 
 typename pTraits<Type>::labelType validComponents
 (
    mesh().template validComponents<Type>()
 );
 
 // The index retrieved is the position in varLoc where the field lies
 // because both varNamesList and varLoc are built sequentially, in the
 // same way
 int bRow = findIndex(varNamesList, rowField);
 int bCol = findIndex(varNamesList, colField);
 
 int i = 0;
 for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
 {
   if (component(validComponents, cmpt) == -1) continue; 
   
   int rowBias = (firstCmpList[bRow]+i)*nGlobalCells;
   int colBias = (firstCmpList[bCol]+i)*nGlobalCells;
   
   if (initTimeFlag || updateA_ || !saveSystem_)
   {
     assemblePetscAb
     (
       A,
       b,
       x,
       matrix,
       cmpt,  
       rowBias,
       colBias
     );
   }
   else
   {
     assemblePetscb
     (
       b,
       matrix,
       cmpt,  
       rowBias
     );
   }
   
   i++;
 } 
 
 isSet = true;
}

template<class Type>
void Foam::coupledSolver::insertEquation
(
  word rowField,
  word colField,
  const tmp<fvMatrix<Type>>& tmatrix
)
{
  insertEquation(rowField, colField, const_cast<fvMatrix<Type>&>(tmatrix()));
  tmatrix.clear();
}
 
 
template<class Type>
void Foam::coupledSolver::insertEquation
(
  word rowField,
  word colField,
  LMatrix<Type>& matrix
)
{
 createSystem();
  
 // The index retrieved is the position in varLoc where the field lies
 // because both varNamesList and varLoc are built sequentially, in the
 // same way
 int bRow = findIndex(varNamesList, rowField);
 int bCol = findIndex(varNamesList, colField);
 
 int fRow = firstCmpList[bRow];
 int fCol = firstCmpList[bCol];
 
 List<labelList>& rel = matrix.rowColList();
 
 forAll(rel, i)
 {   
   int rowBias = (fRow + rel[i][0])*nGlobalCells;
   int colBias = (fCol + rel[i][1])*nGlobalCells;
   direction cmpt = rel[i][2];
   
   if (initTimeFlag || updateA_ || !saveSystem_)
   {
     assemblePetscAb
     (
       A,
       b,
       x,
       matrix,
       cmpt,  
       rowBias,
       colBias
     );
   }
   else
   {
     assemblePetscb
     (
       b,
       matrix,
       cmpt,  
       rowBias
     );
   }
 } 
 
 isSet = true;
}


template<class Type>
void Foam::coupledSolver::insertEquation
(
  word rowField,
  word colField,
  const tmp<LMatrix<Type>>& tLMatrix
)
{
  insertEquation(rowField, colField, const_cast<LMatrix<Type>&>(tLMatrix()));
  tLMatrix.clear();
}

template<class Type>
void Foam::coupledSolver::getSolution()  
{
 Type type(pTraits<Type>::zero); 
   
 List<GeometricField<Type, fvPatchField, volMesh>*>& varList(varTypeList(type));

 forAll(varList, i)
 {
   int g = findIndex(varLoc, {i, ftType(type)}); 
   int fCmp = firstCmpList[g];
   
   typename pTraits<Type>::labelType validComponents
   (
     mesh().template validComponents<Type>()
   );
   
   // varList[i]->correctBoundaryConditions();
   
   int j = 0;
   for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
   {
     if (component(validComponents, cmpt) == -1) continue; 
     
     int rowBias = (fCmp+j)*nGlobalCells;
     transferPetscSolution(x, *varList[i], cmpt, rowBias);  
     
     j++;
   }
   
   varList[i]->correctBoundaryConditions();
   
   // Info << min(*varList[i])<< endl;  
 }
}

template<class Type>
void Foam::coupledSolver::transferPetscSolution
(
  Vec& x,
  Foam::GeometricField<Type, fvPatchField, volMesh>& T,
  int cmpI,
  int rowBias
)
{
 int ilower = this->sharedData.ilower;
 int n = T.size();
 
 // If run in parallel, we cannot VecGetValues() from processors different
 // from the one we are. Therefore we need to scatter and gather the values
 // we need to our processor. 
 if (Pstream::parRun())
 {
   Vec xloc;  
   VecScatter scatter;  
   IS from, to;  
   PetscScalar *values;

   scalarField tt(n, 0.);
   std::vector<int> idx_from(n);
   std::vector<int> idx_to(n);
   for (int i = 0; i < n; i++)
   {
    idx_from[i] = ilower + i + rowBias;
    idx_to[i] = i;
   }

   VecCreateSeq(PETSC_COMM_SELF,n,&xloc);

   ISCreateGeneral(PETSC_COMM_SELF,n,&idx_from[0],PETSC_COPY_VALUES,&from);
   ISCreateGeneral(PETSC_COMM_SELF,n,&idx_to[0],PETSC_COPY_VALUES,&to);

   VecScatterCreate(x,from,xloc,to,&scatter);
   VecScatterBegin(scatter,x,xloc,INSERT_VALUES,SCATTER_FORWARD);
   VecScatterEnd(scatter,x,xloc,INSERT_VALUES,SCATTER_FORWARD);

   VecGetArray(xloc,&values);
   for (int i = 0; i < n; i++)
    tt[i] = values[i];

   T.primitiveFieldRef().replace(cmpI, tt);

   ISDestroy(&from);
   ISDestroy(&to);
   VecScatterDestroy(&scatter);
   VecDestroy(&xloc);
 }
 else
 { 
   std::vector<int> rows(n);
        
   for (int i = 0; i < n; i++)
    rows[i] = ilower + i + rowBias;
    
   scalarField tt(n, 0.);
   VecGetValues(x, n, &rows[0], tt.begin());  
   
   T.primitiveFieldRef().replace(cmpI, tt);
 }
} 

template<class eqType>
void Foam::coupledSolver::assemblePetscAb
(
  Mat& A,
  Vec& b,
  Vec& x,
  eqType& eqn,
  int cmpI,
  int rowBias,
  int colBias
)
{
 int ilower = this->sharedData.ilower;
 
 // Start filling the matrix/vector
 
 //- Off diagonal elements   
 const lduAddressing& lduA = eqn.lduAddr();
 int col; int row;
 int nFaces = lduA.lowerAddr().size();
 
 // Symmetric matrices only have upper(). No need to force creation of lower().
 if (eqn.symmetric())
 {
   scalarField upper(eqn.upper().component(cmpI));
   for (register label face=0; face<nFaces; face++)
   {    
     // ij Off-diagonal
     row = ilower + lduA.upperAddr()[face] + rowBias;
     col = ilower + lduA.lowerAddr()[face] + colBias;
     ierr = MatSetValues(A,1,&row,1,&col,&upper[face],ADD_VALUES);CHKERRV(ierr);
     
     // ji Off-diagonal
     row = ilower + lduA.lowerAddr()[face] + rowBias;
     col = ilower + lduA.upperAddr()[face] + colBias;
     ierr = MatSetValues(A,1,&row,1,&col,&upper[face],ADD_VALUES);CHKERRV(ierr);   
   }
 }
 else
 {
   scalarField lower(eqn.lower().component(cmpI));
   scalarField upper(eqn.upper().component(cmpI)); 
   for (register label face=0; face<nFaces; face++)
   {   
     // ij Off-diagonal
     row = ilower + lduA.upperAddr()[face] + rowBias;
     col = ilower + lduA.lowerAddr()[face] + colBias;
     ierr = MatSetValues(A,1,&row,1,&col,&lower[face],ADD_VALUES);CHKERRV(ierr);
     
     // ji Off-diagonal
     row = ilower + lduA.lowerAddr()[face] + rowBias;
     col = ilower + lduA.upperAddr()[face] + colBias;
     ierr = MatSetValues(A,1,&row,1,&col,&upper[face],ADD_VALUES);CHKERRV(ierr);   
   }
 }
 
 // Diagonal and source  
 scalarField source(eqn.source().component(cmpI));  
 scalarField diag(eqn.diag().component(cmpI));
 int n = source.size(); 
 for (int cellI=0; cellI<n; cellI++) 
 {  
   // Diagonal
   row = ilower + cellI + rowBias; 
   col = ilower + cellI + colBias; 
   ierr = MatSetValues(A,1,&row,1,&col,&diag[cellI],ADD_VALUES);CHKERRV(ierr);
   
   // Source vector   
   ierr = VecSetValues(b,1,&row,&source[cellI],ADD_VALUES);CHKERRV(ierr);
 }
 
 const fvBoundaryMesh& bMesh(mesh().boundary());
 
 //- Contribution from BCs
 forAll(bMesh, patchI)
 {
   scalarField bC = eqn.boundaryCoeffs()[patchI].component(cmpI);  
   scalarField iC = eqn.internalCoeffs()[patchI].component(cmpI);  
   const labelUList& addr = lduA.patchAddr(patchI);      
    
   // Non-coupled
   if (!bMesh[patchI].coupled())
   {
     forAll(addr, facei)
     {
       // Matrix of coefs - Diagonal
       row = ilower + addr[facei] + rowBias;  
       col = ilower + addr[facei] + colBias;  
       ierr = MatSetValues(A,1,&row,1,&col,&iC[facei],ADD_VALUES);CHKERRV(ierr);
       
       // Source vector
       ierr = VecSetValues(b,1,&row,&bC[facei],ADD_VALUES);CHKERRV(ierr);       
     }
   }
   // Processor or processorCyclic.
   // Note: from exchange ops point of view both patch types are equivalent.
   // Only difference is when building the sharedData struct, where the names
   // used to find neig patch in neig proc are different.
   else if (   isType<processorFvPatch>(bMesh[patchI])
            || isType<processorCyclicFvPatch>(bMesh[patchI])
           )
   { 
     // Matrix of coefs - Diagonal
     forAll(addr, facei)
     {      
       row = ilower + addr[facei] + rowBias;
       col = ilower + addr[facei] + colBias;   
       ierr = MatSetValues(A,1,&row,1,&col,&iC[facei],ADD_VALUES);CHKERRV(ierr);       
     } 
     
     // Matrix of coefs - off-diagonal (row -> this processor; col -> other processors)
     // Unique face contribution to off-diag (set instead of add). If use add, then the
     // inter-processor off-diag coefs need to be reset to 0 before.
     forAll(this->sharedData.procInfo, pI)
     {
       if (this->sharedData.procInfo[pI][2] == patchI)
       {
         forAll(bC, facei)
         {        
           double v = -bC[facei];
           int row = this->sharedData.fCo[pI][facei] + rowBias;  
           int col = this->sharedData.fCn[pI][facei] + colBias; 
           ierr = MatSetValues
           (
             A,
             1,
             &row,
             1,
             &col,
             &v,
             ADD_VALUES
           ); 
           CHKERRV(ierr);     
         } 
       }
     }  
   }
   // Cyclic. 
   // Note: cyclic patches (if not empty after decomposition) always contain
   // both own and neig cells in the same processor. Therefore, no data exchange 
   // is needed, and biasing with ilower is enough to account for global indexing.  
   else if (isType<cyclicFvPatch>(bMesh[patchI]))
   { 
     // Row
     const labelList& owfC = mesh().boundaryMesh()[patchI].faceCells();
     
     const fvPatch& cyclicPatch = bMesh[patchI];
    
     const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
        (
            cyclicPatch
        ).neighbFvPatch();
      
     // Col   
     const labelList& nbFC = nbrPatch.patch().faceCells();
     
     forAll(owfC, facei)
     {                  
       // Matrix of coefs - Diagonal  
       row = ilower + owfC[facei] + rowBias;  
       col = ilower + owfC[facei] + colBias; 
       ierr = MatSetValues(A,1,&row,1,&col,&iC[facei],ADD_VALUES);CHKERRV(ierr);  
       
       // Matrix of coefs - off-diagonal
       col = ilower + nbFC[facei] + colBias; 
       double v = -bC[facei];
       ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);       
     } 
       
   }
    
 }  
} 

template<class eqType>
void Foam::coupledSolver::assemblePetscb
(
  Vec& b,
  eqType& eqn,
  int cmpI,
  int rowBias
)
{
 int ilower = this->sharedData.ilower;
 
 // Start filling the vector
 
 const lduAddressing& lduA = eqn.lduAddr();
 int row;
 
 // Diagonal and source  
 scalarField source(eqn.source().component(cmpI));  
 int n = source.size(); 
 for (int cellI=0; cellI<n; cellI++) 
 {  
   // Diagonal
   row = ilower + cellI + rowBias; 
      
   // Source vector   
   ierr = VecSetValues(b,1,&row,&source[cellI],ADD_VALUES);CHKERRV(ierr);
 }
 
 const fvBoundaryMesh& bMesh(mesh().boundary());
 
 //- Contribution from BCs
 forAll(bMesh, patchI)
 {  
   // Non-coupled
   if (!bMesh[patchI].coupled())
   {
     scalarField bC = eqn.boundaryCoeffs()[patchI].component(cmpI);   
     const labelUList& addr = lduA.patchAddr(patchI); 
     
     forAll(addr, facei)
     {
       row = ilower + addr[facei] + rowBias;  
              
       // Source vector
       ierr = VecSetValues(b,1,&row,&bC[facei],ADD_VALUES);CHKERRV(ierr);       
     }
   }    
 }  
} 
 
