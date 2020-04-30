// Special BC for temperature-temperature coupling: nothing in this BC 
// is explicit 

// Special BC for temperature-temperature coupling 
if (!bMesh[patchI].coupled())
{
  // First check if the patch is of the specified type and the if var is scalar (temperature is scalar)
  if (isType<regionCoupledAMIFvPatch>(bMesh[patchI]) && varInfo[rowVarID].typeV == ftscalar)
  {     
    // Now check if field BC is coupledAMI                
    if ( isType<coupledTFvPatchScalarField>( varScalarList[varInfo[rowVarID].localID]->boundaryField()[patchI] ) )
    { 
      continue;
    }
  }
}
