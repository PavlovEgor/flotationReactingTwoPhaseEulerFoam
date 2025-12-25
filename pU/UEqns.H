Info<< "Constructing momentum equations" << endl;

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVolume/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVolume/dimTime);

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    // volScalarField deltaM = flotation.gasSourse();
    // volVectorField momentumFlotationSource = deltaM * (U1 - U2);

    {
        U1Eqn =
        (
            phase1.UEqn()
         ==
           *momentumTransfer[phase1.name()]
        //   + momentumFlotationSource
          + fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        fvOptions.constrain(U1Eqn);
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            phase2.UEqn()
         ==
           *momentumTransfer[phase2.name()]
        //   - momentumFlotationSource
          + fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        fvOptions.constrain(U2Eqn);
        fvOptions.correct(U2);
    }
}
