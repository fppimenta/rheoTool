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
    
    #include "createPhi.H"
 
    // Create the electric model    
    EDFModel elecM(phi);
    
    volVectorField E
    (
      IOobject
      (
 	"E",
 	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
        IOobject::AUTO_WRITE
     ),
     mesh,
     dimensionedVector
      (
          "zero",
          dimensionSet(1, 1, -3, 0, 0, -1, 0),
          vector::zero
      )
    );
 
    
   
    
