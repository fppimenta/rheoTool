{
#   include "continuityErrs.H"

    wordList pcorrTypes
    (
        pd.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<pd.boundaryField().size(); i++)
    {
        if (pd.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", pd.dimensions(), 0.0),
        pcorrTypes
    );

    dimensionedScalar rUAf("(1|A(U))", dimTime/rho.dimensions(), 1.0);

    adjustPhi(phi, U, pcorr);

    mesh.schemesDict().setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rUAf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pdRefCell, pdRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }

#   include "continuityErrs.H"
}
