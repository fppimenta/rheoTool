    Info<< "Reading constitutiveProperties\n" << endl;
    IOdictionary constitutiveProperties
    (
	IOobject
	(
	"constitutiveProperties",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE,
	false
	)
    );

    tensor gradU_(constitutiveProperties.subDict("rheoTestFoamParameters").lookup("gradU"));
    scalarField gammaEpsilonDotL_(constitutiveProperties.subDict("rheoTestFoamParameters").lookup("gammaEpsilonDotL"));
   
    bool ramp_(constitutiveProperties.subDict("rheoTestFoamParameters").lookupOrDefault<Switch>("ramp", false));

    int idx = 0; scalar cnt = 0;

    tensor oldExtraStress(Foam::tensor::zero);

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

     
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
 
tensor gradU = gradU_*gammaEpsilonDotL_[0];

// Change U BC in order to have the desired grad(U)

 //- First get the face distancing, in order to ensure a correct gradient, even if the user changed the mesh by mistake

vector xp(vector::zero), xm(vector::zero), yp(vector::zero), ym(vector::zero), zp(vector::zero), zm(vector::zero);

forAll(U.boundaryField(), patchi)
  {

   if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(1.,0.,0.) ) >= SMALL )
    {	
	    xp = mesh.boundaryMesh()[patchi].faceCentres()[0];
    }
   else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(1.,0.,0.) ) <= -SMALL )
    {
        xm = mesh.boundaryMesh()[patchi].faceCentres()[0];
    }
   else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,1.,0.) ) >= SMALL )
    {
        yp = mesh.boundaryMesh()[patchi].faceCentres()[0];
    }
   else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,1.,0.) ) <= -SMALL )
    {
        ym = mesh.boundaryMesh()[patchi].faceCentres()[0];
    }
   else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,0.,1.) ) >= SMALL )
    {
        zp = mesh.boundaryMesh()[patchi].faceCentres()[0];
    }
   else if( ( mesh.boundaryMesh()[patchi].faceAreas()[0] & vector(0.,0.,1.) ) <= -SMALL )
    {
        zm = mesh.boundaryMesh()[patchi].faceCentres()[0];
    }
	
  }
  
scalar dx = (xp-xm).x();
scalar dy = (yp-ym).y();
scalar dz = (zp-zm).z();
 
  //- Now change the u BC to get the correct gradient

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


//*********************    

#include "createPhi.H"

phi *= 0.; // Homogeneous flow: grad(tau) = 0;

int    nInIter = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<int>("nInIter", 0);

// Create constitutive equation
constitutiveModel constEq(U, phi);

bool isGNF(constEq.isGNF()); 

volSymmTensorField extraStress
(
 IOobject
 (
    "extraStress",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
 ),
 constEq.tauTotal()*1. 
);

// Create file (assume single-processor; won't make sense multiprocessor one 1 cell!)

fileName oFile(runTime.path()/"Report");
OFstream res(oFile);

// Write file header
if (ramp_)
 {
   res << "******* Ramp mode: variable ε̇|γ̇ **********" << 
   nl << "ε̇|γ̇"<< tab <<"extStressXX" << tab << "extStressXY" << tab << "extStressXZ" << tab <<"extStressYY" << tab << "extStressYZ" << tab << "extStressZZ" << tab << "Status" << tab << "Relative_Error" << endl;
 }
else
 {
   res << "******* Transient mode: fixed ε̇|γ̇ = " << gammaEpsilonDotL_[0] << " **********" << 
   nl << "t"<< tab <<"extStressXX" << tab << "extStressXY" << tab << "extStressXZ" << tab <<"extStressYY" << tab << "extStressYZ" << tab << "extStressZZ" << endl;
 }

// Set lambdaMax for VE models (distinguish between multi-mode and single-mode cases)

word CM = constitutiveProperties.subDict("parameters").lookup("type");
scalar lambdaMax(0.0);
scalar lambdaMin(1e20);

if (!isGNF) // doesn't care about ramp, although lambdaMax/Min have no sense if ramp is false
{
   if(CM=="multiMode")
    {
      PtrList<entry> modelEntries(constitutiveProperties.subDict("parameters").lookup("models"));

      forAll(modelEntries, modelI)
       {
         dimensionedScalar lambdaI(modelEntries[modelI].dict().lookup("lambda"));
         lambdaMax = Foam::max(lambdaMax, lambdaI.value());
         lambdaMin = Foam::min(lambdaMin, lambdaI.value());    
       }
    }
   else
    {
      dimensionedScalar lambdatmp(constitutiveProperties.subDict("parameters").lookup("lambda"));
         
      lambdaMax = lambdatmp.value();
      lambdaMin = lambdaMax;
    }
}


// Automatic selection of time-step in ramp mode

if (ramp_ && !isGNF) 
   {
       nInIter = 1; // Doesn't make sense when seeking for steady values
       runTime.setDeltaT((lambdaMin+lambdaMax)*0.5/gammaEpsilonDotL_[0]);
       runTime.setEndTime(1e8); // Just to avoid an early exit
   }
else if(ramp_ && isGNF) 
   {
       nInIter = 1; // Doesn't make sense when seeking for steady values
       runTime.setDeltaT(lambdaMin/20);
       lambdaMax = 1; // Any value, since it will exit from the error
       lambdaMin = 1; // Any value, since it will exit from the error
       runTime.setDeltaT(1); // Any value, since no equation is solved
       runTime.setEndTime(lambdaMax*1e8); // Just to avoid an early exit
   }

