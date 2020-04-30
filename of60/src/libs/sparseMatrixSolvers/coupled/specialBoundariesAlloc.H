//- Note: each BC here must 'continue' the loop once it is entered to ensure it does
// not contribute twice to the matrix of coeffs

// Special BC for temperature-temperature coupling 
if (!bMesh[patchI].coupled())
{
  // First check if the patch is of the specified type 
  // Note: here we consider that all fvPatchFields in this patch will couple to other field (over-sized)  
  if (isType<regionCoupledAMIFvPatch>(bMesh[patchI]))
  {           
     const fvPatch& pfvPatch = bMesh[patchI];  
        
     const regionCoupledBaseFvPatch& camipp = refCast<const regionCoupledAMIFvPatch>(pfvPatch);     
     const regionCoupledBaseFvPatch& neicamipp = refCast<const regionCoupledAMIFvPatch>(pfvPatch).neighbPatch();   
     
     const labelList& ownFC = camipp.faceCells();   
          
     if (camipp.owner())
     {                       
       const labelListList& srcAd = camipp.AMI().srcAddress();
 
       // Each patch of the AMI interface in a different processor
       if (camipp.AMI().singlePatchProc() == -1) 
       {            
         forAll(srcAd, facei)
         {
           forAll(srcAd[facei], kk)
           {
             if (camipp.AMI().applyLowWeightCorrection())
             {            
               if (camipp.AMI().srcWeightsSum()[facei] < camipp.AMI().lowWeightCorrection())
               {
                 maxInProcFaces[ownFC[facei]] += 1;
               }
               else
               {
                 maxOutProcFaces[ownFC[facei]] += 1;                 
               }
             }
             else
             {
               maxOutProcFaces[ownFC[facei]] += 1;                 
             }
           }            
         }      
       }
       // Both patches of the AMI interface in same processor
       else
       {
         forAll(srcAd, facei)
         {
           forAll(srcAd[facei], kk)
           {
             maxInProcFaces[ownFC[facei]] += 1;
           }            
         }
       }      
     }  
     else
     {      
       const labelListList& tgtAd = neicamipp.AMI().tgtAddress();
      
       // Each patch of the AMI interface in a different processor
       if (neicamipp.AMI().singlePatchProc() == -1) 
       {             
         forAll(tgtAd, facei)
         {
           forAll(tgtAd[facei], kk)
           {
             if (neicamipp.AMI().applyLowWeightCorrection())
             {           
               if (neicamipp.AMI().tgtWeightsSum()[facei] < neicamipp.AMI().lowWeightCorrection())
               {
                 maxInProcFaces[ownFC[facei]] += 1;
               }
               else
               {                 
                 maxOutProcFaces[ownFC[facei]] += 1;
               }
             }
             else
             {                 
               maxOutProcFaces[ownFC[facei]] += 1;
             }
           }            
         }      
       }
       // Both patches of the AMI interface in same processor
       else
       {
         forAll(tgtAd, facei)
         {
           forAll(tgtAd[facei], kk)
           {
             maxInProcFaces[ownFC[facei]] += 1;
           }           
         }
       }      
     }   
       
  // If the BC is of the specified type, then proceed to next BC. Not adding 'continue' 
  // would make us enter again in case 'non-coupled'.
  continue;
    
 }  
}
