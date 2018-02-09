// Write

word status("Exceed_Niter");
if (relError < 1e-8)
{
 status= "Converged";
}

res << gammaEpsilonDotL_[idx] << tab 
    << extraStress[0].component(0) << tab 
    << extraStress[0].component(1) << tab 
    << extraStress[0].component(2) << tab 
    << extraStress[0].component(3) << tab 
    << extraStress[0].component(4) << tab 
    << extraStress[0].component(5) << tab 
    << status << tab
    << relError
    << endl;  

// Restart 

idx++;
cnt = 0;

if ( idx == gammaEpsilonDotL_.size() ) // Does the last value has been reached?
  {
     Info << nl << "************************" << 
     nl << nl << "Your ramp is finished !" << 
     nl << nl << "************************" << nl << endl;

     return 0;
  }

runTime.setDeltaT((lambdaMin+lambdaMax)*0.5/gammaEpsilonDotL_[idx]); // Good approach for some cases

gradU = gradU_*gammaEpsilonDotL_[idx];

forAll(U.boundaryField(), patchi)
{

   if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(1.,0.,0.) ) >= SMALL )
    {	
	fixedValueFvPatchVectorField& Upatch = refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchi]);
    	Upatch == vectorField(Upatch.size(), vector(gradU.xx()*.5*dx,gradU.xy()*.5*dx,gradU.xz()*.5*dx));
    }
    else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(1.,0.,0.) ) <= -SMALL )
    {
        fixedValueFvPatchVectorField& Upatch = refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchi]);
	Upatch == vectorField(Upatch.size(), -vector(gradU.xx()*.5*dx,gradU.xy()*.5*dx,gradU.xz()*.5*dx));
    }
    else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,1.,0.) ) >= SMALL )
    {
        fixedValueFvPatchVectorField& Upatch = refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchi]);
	Upatch == vectorField(Upatch.size(), vector(gradU.yx()*.5*dy,gradU.yy()*.5*dy,gradU.yz()*.5*dy));
    }
    else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,1.,0.) ) <= -SMALL )
    {
        fixedValueFvPatchVectorField& Upatch = refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchi]);
	Upatch == vectorField(Upatch.size(), -vector(gradU.yx()*.5*dy,gradU.yy()*.5*dy,gradU.yz()*.5*dy));
    }
    else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,0.,1.) ) >= SMALL )
    {
        fixedValueFvPatchVectorField& Upatch = refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchi]);
	Upatch == vectorField(Upatch.size(), vector(gradU.zx()*.5*dz,gradU.zy()*.5*dz,gradU.zz()*.5*dz));
    }
    else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,0.,1.) ) <= -SMALL )
    {
        fixedValueFvPatchVectorField& Upatch = refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchi]);
	Upatch == vectorField(Upatch.size(), -vector(gradU.zx()*.5*dz,gradU.zy()*.5*dz,gradU.zz()*.5*dz));
    }
	
}



