{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());
    
    UPtrList<volScalarField>& Y1 = phase1.YActiveRef();
    UPtrList<volScalarField>& Y2 = phase2.YActiveRef();

    FlotationSoursePart FS(flotation.flotationMassTransfer());
    // volScalarField Q_gas = 1e-3 * flotation.gasSourseKa();

    // int indGas, indLiq, indParGas, indParLiq;
    
    if (!phase1.pure())
    {

        forAll(Y1, i)
        {
            fvScalarMatrix Y1iEqn
            (
                phase1.YiEqn(Y1[i])
             ==
                *massTransfer[Y1[i].name()]
               + fvOptions(alpha1, rho1, Y1[i])
            );

            if (Y1[i].name() == "particle.gas" && delay.value() < runTime.value() ){
                Y1iEqn = (  fvm::ddt(alpha1, rho1, Y1[i])
                            + fvm::div(phase1.alphaRhoPhi(), Y1[i], "div(alphaRhoPhi.particle,Yi)")
                            ==
                            *massTransfer[Y1[i].name()]
                            + fvOptions(alpha1, rho1, Y1[i])
                            + FS.Q_phase1
                        );
                // Y1iEqn -= FS.Q_phase1;
                // Y1iEqn -= Q_gas;
                // Info << "Q_gas : " << average(Q_gas).value() << endl;
            }

            Y1iEqn.relax();
            Y1iEqn.solve("Yi");

            if (Y1[i].name() == "particle.gas" && delay.value() < runTime.value() ){
                forAll(alpha1, cellI){
                    if (alpha1[cellI] < 1e-3)
                        Y1[i][cellI] = 0;
                }

                particleMass_inAir = alpha1 * rho1 * Y1[i];
            }

        }

    }
    if (!phase2.pure())
    {

        forAll(Y2, i)
        {
            fvScalarMatrix Y2iEqn
            (
                phase2.YiEqn(Y2[i])
             ==
               *massTransfer[Y2[i].name()]
              + fvOptions(alpha2, rho2, Y2[i])
            );


            if (Y2[i].name() == "particle.liquid" && delay.value() < runTime.value() ){
                Y2iEqn -= FS.Q_phase2;
                // Y2iEqn += Q_gas;
            }

            Y2iEqn.relax();
            Y2iEqn.solve("Yi");

            if (Y2[i].name() == "particle.liquid" && delay.value() < runTime.value() ){
                forAll(alpha2, cellI){
                    if (alpha2[cellI] < 1e-3)
                        Y2[i][cellI] = 0;
                }
            }
        }
    }
}
