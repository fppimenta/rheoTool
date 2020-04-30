//- Note: each BC here must 'continue' the loop once it is entered to ensure it does
// not contribute twice to the matrix of coeffs

// Special BC for temperature-temperature coupling 
if (!bMesh[patchI].coupled())
{
  // First check if the patch is of the specified type and the if var is scalar (temperature is scalar)
  if (isType<regionCoupledAMIFvPatch>(pfvPatch) && varInfo[rowVarID].typeV == ftscalar)
  {     
    // Now check if field BC is coupledAMI                
    if ( isType<coupledTFvPatchScalarField>( varScalarList[varInfo[rowVarID].localID]->boundaryField()[patchI] ) )
    { 
     // camipp and neicampi always lay in the same processor, but one of them can be empty
     // if the other half of a given AMI was sent to another processor. In pratice, nothing
     // happens bellow if camipp is zero-sized.
     const regionCoupledBaseFvPatch& camipp = refCast<const regionCoupledAMIFvPatch>(pfvPatch);     
     const regionCoupledBaseFvPatch& neicamipp = refCast<const regionCoupledAMIFvPatch>(pfvPatch).neighbPatch();   
     
     const labelList& ownFC = camipp.faceCells();   
     const labelList& neiFC = neicamipp.faceCells();    
      
     const coupledTFvPatchScalarField& camiPatchField = refCast<const coupledTFvPatchScalarField>
     (
      varScalarList[varInfo[rowVarID].localID]->boundaryField()[patchI]
     );  
     
     // The ID of T from the neigb region
     int bCol = findField(camiPatchField.nbrFieldName(), camiPatchField.nbrMeshName());
 
     // The column bias
     int fRow = varInfo[bCol].firstElem;
 
     // The neiFC gets ilower from the proc/mesh where neicamipp is. This is ensured after
     // the call to distribute(neiFCproc).
     labelList neiFCproc(neiFC+this->sharedData[meshList[varInfo[bCol].meshID].ID].ilower); 
      
     // AMIs are shared by at least 2 patches, but the AMI interpolator is only accessible
     // from one of them, which is considered the "owner" patch.
     if (camipp.owner())
     {      
         // Distribute faceCells of neighbour (tgt) patch if parallel. We get
         // the tgt faceCells transfered to the proc where the src patch (camipp) is.
         if (camipp.AMI().singlePatchProc() == -1) 
          camipp.AMI().tgtMap().distribute(neiFCproc); 
         
         const labelListList& srcAd = camipp.AMI().srcAddress();
         const scalarListList& srcW = camipp.AMI().srcWeights();
         forAll(srcAd, k)
         {             
           // Matrix of coefs - Diagonal  
           row = ilower + ownFC[k] + rowBias; 
           col = ilower + ownFC[k] + colBias;  
           ierr = MatSetValues(A,1,&row,1,&col,&iC[k],ADD_VALUES);CHKERRV(ierr); 
         
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (camipp.AMI().applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (camipp.AMI().srcWeightsSum()[k] < camipp.AMI().lowWeightCorrection())
             {
               col = ilower + ownFC[k] + colBias;
               double v = -bC[k];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr); 
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(srcAd[k], kk)
               {                    
                col = neiFCproc[srcAd[k][kk]] + fRow;  
                double v = -bC[k]*srcW[k][kk];
                ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);         
               }
             }                                           
           }
           else
           {
              // Matrix of coefs - off-diagonal
              forAll(srcAd[k], kk)
              {                        
               col = neiFCproc[srcAd[k][kk]] + fRow; 
               double v = -bC[k]*srcW[k][kk];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);           
              }   
           }
         }
     }
     else
     {        
         // Distribute faceCells of neighbour (src) patch if parallel. We get
         // the src faceCells transfered to the proc where the tgt patch (camipp) is
         if (neicamipp.AMI().singlePatchProc() == -1)
           neicamipp.AMI().srcMap().distribute(neiFCproc);  
       
         const labelListList& tgtAd = neicamipp.AMI().tgtAddress();
         const scalarListList& tgtW = neicamipp.AMI().tgtWeights();
         forAll(tgtAd, k)
         {                         
           // Matrix of coefs - Diagonal 
           row = ilower + ownFC[k] + rowBias; 
           col = ilower + ownFC[k] + colBias;     
           ierr = MatSetValues(A,1,&row,1,&col,&iC[k],ADD_VALUES);CHKERRV(ierr);            
                     
           // If the applyLowWeightCorrection option is enabled, the zero-gradient
           // condition must be applied when the weightSum is less than a treshold.  
           if (neicamipp.AMI().applyLowWeightCorrection())
           {
             // Apply implicit zero-gradient
             if (neicamipp.AMI().tgtWeightsSum()[k] < neicamipp.AMI().lowWeightCorrection())
             {
               col = ilower + ownFC[k] + colBias;
               double v = -bC[k];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr); 
             }
             // Distribute weighted coefficients 
             else
             {
               // Matrix of coefs - off-diagonal
               forAll(tgtAd[k], kk)
               {                
                col = neiFCproc[tgtAd[k][kk]] + fRow; 
                double v = -bC[k]*tgtW[k][kk];
                ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);          
               }  
             }                                           
           }
           else
           {
              // Matrix of coefs - off-diagonal
              forAll(tgtAd[k], kk)
              {                
               col = neiFCproc[tgtAd[k][kk]] + fRow; 
               double v = -bC[k]*tgtW[k][kk];
               ierr = MatSetValues(A,1,&row,1,&col,&v,ADD_VALUES);CHKERRV(ierr);
              }   
           }              
         }       
    }
    
    // If the BC is of the specified type, then proceed to next BC. Not adding 'continue' 
    // would make us enter again in case 'non-coupled'.
    continue;
   }   
 }  
}
