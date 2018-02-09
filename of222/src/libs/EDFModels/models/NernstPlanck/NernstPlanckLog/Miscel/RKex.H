  dimensionedScalar dt = phi().mesh().time().deltaT();
  
  volScalarField y0(lci.oldTime());
  volScalarField y00(lci.oldTime().oldTime());
  
  surfaceScalarField cif0(Foam::exp(fvc::interpolate(y0)));
  surfaceScalarField cif00(Foam::exp(fvc::interpolate(y00)));
  
  // k1
  
  volScalarField 
  k1 = dt*
        (
           fvc::laplacian(species_[i].Di()*cif0, y0, "laplacian(D,ci)")   	  
         + fvc::laplacian(ca*cif0, phiE_.oldTime()+psi_.oldTime(), "laplacian(EMigD,phiE)")
        )/Foam::exp(y0); 
   
  // km1      
  volScalarField 
  km1 = dt * 
        (
           fvc::laplacian(species_[i].Di()*cif00, y00, "laplacian(D,ci)")   	  
         + fvc::laplacian(ca*cif00, phiE_.oldTime().oldTime()+psi_.oldTime().oldTime(), "laplacian(EMigD,phiE)")
        )/Foam::exp(y00); 
  
  // k2       
        
  volScalarField 
  y2 = y0 + 5.*k1/12; 
  
  surfaceScalarField cif2(Foam::exp(fvc::interpolate(y2))); 
  
  volScalarField 
  k2 = dt * 
        (
           fvc::laplacian(species_[i].Di()*cif2, y2, "laplacian(D,ci)")   	  
         + fvc::laplacian(ca*cif2, phiE_.oldTime()+psi_.oldTime(), "laplacian(EMigD,phiE)")
        )/Foam::exp(y2);  
        
  // km2
  
  volScalarField 
  ym2 = y00 + 5.*km1/12; 
  
  surfaceScalarField cifm2(Foam::exp(fvc::interpolate(ym2))); 
  
  volScalarField 
  km2 = dt * 
        (
           fvc::laplacian(species_[i].Di()*cifm2, ym2, "laplacian(D,ci)")   	  
         + fvc::laplacian(ca*cifm2, phiE_.oldTime().oldTime()+psi_.oldTime().oldTime(), "laplacian(EMigD,phiE)")
        )/Foam::exp(ym2);  
        
  // yn1
  
  lci = y0 + (k1 + km1 + 2.*k2 - 2.*km2)*0.5;
           
