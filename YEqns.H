{
    UPtrList<volScalarField>& Y1 = phase1.YActiveRef();
    UPtrList<volScalarField>& Y2 = phase2.YActiveRef();
    FlotationSoursePart FS(flotation.flotationMassTransfer());
    scalar alphaThreshold = 1e-3;
    if (!phase1.pure())
    {
        forAll(Y1, i)
        {
            fvScalarMatrix Y1iEqn
            (
                fvm::ddt(alpha1, rho1, Y1[i])
                + fvm::div(phase1.alphaRhoPhi(), Y1[i], "div(alphaRhoPhi." + Y1[i].name() + ",Yi)")
             ==
               fvOptions(alpha1, rho1, Y1[i])
            );
            if (Y1[i].name() == "particle.gas" && delay.value() < runTime.value() )
            {
                Y1iEqn -= FS.Q_phase1;
            }
            Y1iEqn.relax();
            Y1iEqn.solve("Yi");
            if (Y1[i].name() == "particle.gas" && delay.value() < runTime.value() )
            {
                forAll(alpha1, cellI)
                {
                    if ((alpha1[cellI] < alphaThreshold) || (Y1[i][cellI] < 0))
                    {
                        Y1[i][cellI] = 0;
                    }
                }
                forAll(Y1[i].boundaryField(), patchI)
                {
                    forAll(Y1[i].boundaryField()[patchI], faceI)
                    {
                        label cellI = Y1[i].boundaryField()[patchI].patch().faceCells()[faceI];
                        if ((alpha1[cellI] < alphaThreshold) || (Y1[i].boundaryFieldRef()[patchI][faceI] < 0))
                        {
                            Y1[i].boundaryFieldRef()[patchI][faceI] = 0.0;
                        }
                    }
                }
                particleMass_inAir = alpha1 * rho1 * Y1[i];

                Info << average(Y1[i]).value() << " " <<min(Y1[i]).value() << " " << max(Y1[i]).value() << nl;
            }
        }
    }
    if (!phase2.pure())
    {
        forAll(Y2, i)
        {
            fvScalarMatrix Y2iEqn
            (
                fvm::ddt(alpha2, rho2, Y2[i])
                + fvm::div(phase2.alphaRhoPhi(), Y2[i], "div(alphaRhoPhi." + Y2[i].name() + ",Yi)")
             ==
               fvOptions(alpha2, rho2, Y2[i])
            );
            if (Y2[i].name() == "particle.liquid" && delay.value() < runTime.value() )
            {
                Y2iEqn -= FS.Q_phase2;
            }
            Y2iEqn.relax();
            Y2iEqn.solve("Yi");
            if (Y2[i].name() == "particle.liquid" && delay.value() < runTime.value() )
            {
                forAll(alpha2, cellI)
                {
                    if (alpha2[cellI] < alphaThreshold)
                    {
                        Y2[i][cellI] = 0;
                    }
                }
                forAll(Y2[i].boundaryField(), patchI)
                {
                    forAll(Y2[i].boundaryField()[patchI], faceI)
                    {
                        label cellI = Y2[i].boundaryField()[patchI].patch().faceCells()[faceI];
                        if (alpha2[cellI] < alphaThreshold)
                        {
                            Y2[i].boundaryFieldRef()[patchI][faceI] = 0.0;
                        }
                    }
                }
            }
        }
    }
    fluid.correct();
}
