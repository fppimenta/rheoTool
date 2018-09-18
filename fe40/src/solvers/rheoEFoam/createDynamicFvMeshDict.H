IOobject dictHeader
 (
       "dynamicMeshDict",
       runTime.constant(),
       runTime,
       IOobject::MUST_READ,
       IOobject::NO_WRITE,
       false
 );
 
// If dict does not exist, we assume the mesh is static
// and we create/write a dummy dict 
if (!dictHeader.headerOk())
 {    
   IOobject dictHeader2
   (
       "dynamicMeshDict",
       runTime.constant(),
       runTime,
       IOobject::NO_READ,
       IOobject::NO_WRITE,
       false
   );
        
   IOdictionary dictDyM(dictHeader2);
   dictDyM.add("dynamicFvMesh", "staticFvMesh");
   
   dictDyM.writeObject(IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED);  
 }

