for (int jj = 0; jj<nInIter; jj++) 
{      
    Info << "Inner iteration:  " << jj << nl << endl; 
    
    if (solveCoupled)
    {            
      // Add/Solve electric equations
      elecM.correct();
           
      if (!frozenFlow)
      {          
        // When we enter here we know that there is no gravity,
        // because this has been checked earlier. phig = 0 and p = p_rgh
      
        #include "pUEqn.H" 
        
        // Add/solve constitutive equation 
        constEq.correct();      
        
        // Solve all coupled
        cpsFluid[i].solve(); 
        
        phi = fvc::flux(U) + pRC - fvc::snGrad(p_rgh)*fvc::interpolate(rAU)*mesh.magSf(); 
             
        #include "continuityErrs.H"
        
        fvOptions.correct(U);
        
        p == p_rgh;

        if (p_rgh.needReference())
        {
           p += dimensionedScalar
           (
             "p",
              p.dimensions(),
              pRefValue - getRefCellValue(p, pRefCell)
           );
           p_rgh = p;
        }            
      }        
    } 
    else
    {
      if (!frozenFlow)
      {
        #include "UEqn.H"
        #include "pEqn.H"
                
        // ---- Solve constitutive equation ----	
        constEq.correct();
      }
      
      // ---- Update electric terms ----
      elecM.correct();
    } 
    
    // ---- Update thermo ---
    // Can be solved coupled/segregated independent of p-U
    thermo.correct(U,phi,constEq.tauTotal(),cpsT,fvOptions,pimple.nNonOrthCorr());  
    
    // --- Passive Scalar transport
    if (sPS)
    {
      #include "CEqn.H"
    }        
 
    postProc.update();
} 
