// Read molecule properties
   
   IOdictionary molcProperties
    (
        IOobject
        (
            "moleculesControls",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

   bool   writeContFields_ = molcProperties.subDict("externalFlow").lookupOrDefault<Switch>("writeFields", true);
    
   scalar execTimeLagrang(0.);
    
   sPCloudInterface molecules(U);
