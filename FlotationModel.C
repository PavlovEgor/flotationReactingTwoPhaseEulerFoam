#include "FlotationModel.H"


FlotationModel::FlotationModel(const fvMesh& mesh, twoPhaseSystem& fluid, dimensionedScalar& delay)
    : mesh_(mesh), 
    fluid_(fluid),
    rhoP("rhoP", dimDensity, 0.0), 
    d_p("d_p", dimLength, 0.0),
    d_b("d_b", dimLength, 0.0),
    nu_F("nu_F", dimViscosity, 0.0),
    sigma("sigma", dimForce / dimLength, 0.0)
{


    Info << "Reading flotationProperties" << endl;
    IOdictionary flotationProperties
    (
        IOobject
        (
            "flotationProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    delay = dimensionedScalar("delay", dimTime, readScalar(flotationProperties.lookup("delay")));

    // Particle Properties
    rhoP     = dimensionedScalar("rhoP", dimDensity, readScalar(flotationProperties.lookup("rho0"))); // density
    d_p      = dimensionedScalar("d_p", dimLength, readScalar(flotationProperties.lookup("d_p"))); // diameter 

    // Liquid-Gas Properties
    d_b =
    dimensionedScalar(
        "d_b",
        dimLength,
        readScalar(flotationProperties.lookup("d_b"))
    ); // diameter of bubble
    nu_F =
    dimensionedScalar(
        "nu_F",
        dimViscosity,
        readScalar(flotationProperties.lookup("nu_F"))
    );
    sigma =
    dimensionedScalar(
        "sigma",
        dimForce / dimLength, // Размерность поверхностного натяжения: [N/m] или [kg/s²]
        readScalar(flotationProperties.lookup("sigma"))
    );

    // Particle-Liquid-Gas Properties

    theta_c_deg = readScalar(flotationProperties.lookup("theta_c")); // contact angle in degree
    theta_c_rad = M_PI * theta_c_deg / 180.0; // contact angle in rad


    N_max = (((d_b / d_p) * (d_b / d_p)) / 2).value();

    name_of_particles = flotationProperties.getOrDefault<word>("name_of_particles", "particle");

    //////////////////////////////////////////////////////
    //                     Attachment                  //

    const dictionary& attachmentModelDict =
        flotationProperties.subDict("attachmentModel");

    attachmentModel = attachmentModelDict.getOrDefault<word>("type", "YoonLuttrell");
    inductionTimeModel = attachmentModelDict.getOrDefault<word>("inductionTimeModel", "KohSchwarz");

    if(inductionTimeModel == "Dai"){
        A = readScalar(attachmentModelDict.lookup("A"));
        B = readScalar(attachmentModelDict.lookup("B"));
    }

    //////////////////////////////////////////////////////
    //                     Collision                    //

    const dictionary& collisionModelDict =
        flotationProperties.subDict("collisionModel");

    collisionFrequencyModel = collisionModelDict.getOrDefault<word>("collisionFrequencyModel", "SaffmanTurner");
    collisionEfficiencyModel = collisionModelDict.getOrDefault<word>("collisionEfficiencyModel", "Sutherland");

    //////////////////////////////////////////////////////
    //                     Stability                    //

    const dictionary& stabilityModelDict =
        flotationProperties.subDict("stabilityModel");

    stabilityEfficiencyModel = stabilityModelDict.getOrDefault<word>("stabilityEfficiencyModel", "Schulze");
    detachmentFrequencyModel = stabilityModelDict.getOrDefault<word>("detachmentFrequencyModel", "MikaFuerstenau");

    C_1 = readScalar(stabilityModelDict.lookup("C_1"));
    A_s = readScalar(stabilityModelDict.lookup("A_s"));

    //////////////////////////////////////////////////////

}

FlotationModel::~FlotationModel(){}

volScalarField FlotationModel::collisionEfficiency(){
    volScalarField Ec(
                IOobject
                    (
                        "collisionEfficiency",
                        mesh_.time().constant(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ), 
                mesh_, 
                dimensionedScalar("ONE", dimless, 1.0));

    if (collisionEfficiencyModel == "Sutherland"){

        Ec = 3 * (d_p / d_b);
        // return Ec;


    } else if (collisionEfficiencyModel == "Gaudin"){

        Ec = (3.0 / 2.0) * (d_p / d_b) * (d_p / d_b);

    } else if (collisionEfficiencyModel == "YoonLuttrell"){

        volScalarField Re_b = Foam::mag(fluid_.phase1().URef()) * d_b / nu_F;
        Ec = ((3.0 / 2.0) + (4.0 * Foam::pow(Re_b, 0.72) / 15.0)) * (d_p / d_b) * (d_p / d_b);

    } else if (collisionEfficiencyModel == "WeberPaddock"){

        volScalarField Re_b = Foam::mag(fluid_.phase1().URef()) * d_b / nu_F;
        Ec = Foam::sin(theta_c_rad) * (1 + d_p / d_b) * (1 + d_p / d_b) * (Foam::mag(fluid_.phase2().URef()) / Foam::mag(fluid_.phase1().URef()));

    } else {
        Info << "Unknown Collision model: " << collisionEfficiencyModel << nl << exit(FatalError);
    }

    return Ec;

}


volScalarField FlotationModel::attachmentEfficiency(){

    dimensionedScalar t_i("t_i", dimTime, 0.0); 

    if (inductionTimeModel == "KohSchwarz"){
        t_i = dimensionedScalar("t_i", dimTime, (75.0 / theta_c_deg) * Foam::pow(d_p.value(), 0.6)); 
    } else if(inductionTimeModel == "Dai"){
        t_i = dimensionedScalar("t_i", dimTime, A * Foam::pow(d_p.value(), B)); 
    }  else{
        Info << "Unknown Induction time model: " << inductionTimeModel << nl << exit(FatalError);
    }

    if (attachmentModel == "YoonLuttrell"){

        volScalarField theta_a = 2 * Foam::atan(Foam::exp(-3 * Foam::mag(fluid_.phase1().URef()) * t_i / (d_b * ((d_b / d_p) + 1))));
        return Foam::pow(Foam::sin(theta_a), 2);

    }else if (attachmentModel == "DobbyFinch"){

        volScalarField theta_a = 2 * Foam::atan(Foam::exp(-t_i * (Foam::mag(fluid_.phase1().URef()) + Foam::mag(fluid_.phase2().URef()) * (2.0 + (Foam::pow(d_b / (d_p + d_b), 3)) / (d_p + d_b)))));
        volScalarField sinTheta_c = theta_a;
        volScalarField Re_b = Foam::mag(fluid_.phase1().URef()) * d_b / nu_F;
        forAll(Re_b, i){
            if (Re_b[i] > 20.0 ){
                sinTheta_c = Foam::sin(78.1 - 7.37 * Foam::log(Re_b));
            } else if (Re_b[i] > 1.0) {
                sinTheta_c = Foam::sin(85.5 - 12.49 * Foam::log(Re_b));
            } else {
                sinTheta_c = Foam::sin(85.0 - 2.5 * Foam::log(Re_b));
            }
        }

        return Foam::pow(Foam::sin(theta_a) / sinTheta_c, 2);

    }else {
        Info << "Unknown Attachment model: " << attachmentModel << nl << exit(FatalError);
    }
}

volScalarField FlotationModel::stabilityEfficiency(const volScalarField& epsilon){
    if (stabilityEfficiencyModel == "Schulze"){

        dimensionedScalar g("g", dimAcceleration, 9.81);

        volScalarField down = (d_p * d_p) * ((rhoP - fluid_.phase2().rho()) * g + rhoP * (1.9 * Foam::pow(epsilon, 2.0/3.0) * Foam::pow((d_p + d_b)/2, -1.0/3.0)))
        + (3.0 / 2.0) * d_p * (Foam::sin(M_PI - theta_c_rad/2) * Foam::sin(M_PI - theta_c_rad/2))
        * ((4 * sigma / d_b) - d_b * fluid_.phase2().rho() * g);
        volScalarField oneDivBo_prime = Foam::mag(6 * sigma * Foam::sin(M_PI + theta_c_rad/2) * Foam::sin(M_PI + theta_c_rad/2)) / down;

        return 1 - Foam::exp(A_s * (1 - (oneDivBo_prime * pos(oneDivBo_prime - 1) + 1 * (1 - pos(oneDivBo_prime - 1)))));

    } else{
        Info << "Unknown stabilityEfficiencyModel: " << stabilityEfficiencyModel << nl << exit(FatalError);
    }
}

volScalarField FlotationModel::K_attachment(const volScalarField& epsilon, const volScalarField& Es){

    volScalarField Z
    (IOobject
        (
            "Z",
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ), 
    mesh_, 
    dimensionedScalar("ZERO_freq", dimensionSet(0, 3, -1, 0, 0, 0, 0), 0.0));

    if (collisionFrequencyModel == "SaffmanTurner"){
        Z = Foam::sqrt((8.0 * M_PI) / (15.0 * nu_F)) * (((d_p + d_b) / 2) * ((d_p + d_b) / 2) * ((d_p + d_b) / 2)) * Foam::sqrt(epsilon);
    } else if (collisionFrequencyModel == "SchubertBischofberger"){
        Z = 5.0 * (d_p * d_p / 4.0) * 
            (0.4 *Foam::pow(epsilon, 4.0/9.0) * Foam::pow(d_b, 7.0/9.0) / Foam::pow(nu_F, 1.0/3.0)) * 
            ((fluid_.phase2().rho() - fluid_.phase1().rho()) / fluid_.phase2().rho());
    } else{
        Info << "Unknown collisionFrequencyModel: " << collisionFrequencyModel << nl << exit(FatalError);
    }

    volScalarField k1 = Z * collisionEfficiency() * 
                            attachmentEfficiency() *
                            Es;

    volScalarField beta = (fluid_.phase1().rho()/ fluid_.phase2().rho()) * 
                            (fluid_.phase1().Y(name_of_particles) / N_max) * 
                            (d_b * d_b * d_b / (d_p * d_p * d_p));

    // Info << "beta : " << average(beta).value() << endl;
    // Info << "fluid_.phase1().Y(name_of_particles) : " << average(fluid_.phase1().Y(name_of_particles)).value() << endl;


    return k1 * (1 - beta) * fluid_.phase2() / (( M_PI / 6.0) * (d_b * d_b * d_b));
}

volScalarField FlotationModel::K_detachment(const volScalarField& epsilon, const volScalarField& Es){

    volScalarField Z_prime
    (IOobject
        (
            "Z",
            mesh_.time().constant(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ), 
    mesh_, 
    dimensionedScalar("ZERO_freq", dimless/dimTime, 0.0));

    if (detachmentFrequencyModel == "MikaFuerstenau"){
        Z_prime = Foam::sqrt(C_1) * Foam::pow(d_p + d_b, -2.0/3.0) * Foam::cbrt(epsilon);
    } else{
        Info << "Unknown detachmentFrequencyModel: " << detachmentFrequencyModel << nl << exit(FatalError);
    }

    return Z_prime * (1 - Es) / N_max;
}

FlotationSoursePart FlotationModel::flotationMassTransfer(){

    // const volScalarField& epsilon = mesh_.lookupObject<volScalarField>("epsilon.liquid");

    const volScalarField& k = mesh_.lookupObject<volScalarField>("k.liquid");
    const volScalarField& omega = mesh_.lookupObject<volScalarField>("omega.liquid");
    volScalarField epsilon = k * omega * scalar(0.09);

    const volScalarField& Es = stabilityEfficiency(epsilon);

    volScalarField Ka = K_attachment(epsilon, Es);
    volScalarField Kd = K_detachment(epsilon, Es);

    volScalarField alphaRho1 = fluid_.phase1() * fluid_.phase1().rho();
    volScalarField alphaRho2 = fluid_.phase2() * fluid_.phase2().rho();


    // fvScalarMatrix Q_gas = fvm::Sp(-Kd * fluid_.phase1() * fluid_.phase1().rho(), fluid_.phase1().Y(name_of_particles)) + 
    //                                 Ka * fluid_.phase2() * fluid_.phase2().rho() * fluid_.phase2().Y(name_of_particles);

    // fvScalarMatrix Q_liquid = Kd * fluid_.phase1() * fluid_.phase1().rho() * fluid_.phase1().Y(name_of_particles) - 
    //                   fvm::Sp(Ka * fluid_.phase2() * fluid_.phase2().rho(), fluid_.phase2().Y(name_of_particles));

    // fvScalarMatrix Q_gas = fvm::Sp(dimensionedScalar("ZERO", dimDensity/dimTime, 0.0), fluid_.phase1().Y(name_of_particles)) -
    //                                 Kd * fluid_.phase1() * fluid_.phase1().rho() * fluid_.phase1().Y(name_of_particles) + 
    //                                 Ka * fluid_.phase2() * fluid_.phase2().rho() * fluid_.phase2().Y(name_of_particles);

    // fvScalarMatrix Q_liquid = fvm::Sp(dimensionedScalar("ZERO", dimDensity/dimTime, 0.0), fluid_.phase2().Y(name_of_particles)) +
    //                                 Kd * fluid_.phase1() * fluid_.phase1().rho() * fluid_.phase1().Y(name_of_particles)- 
    //                                 Ka * fluid_.phase2() * fluid_.phase2().rho() * fluid_.phase2().Y(name_of_particles);

    fvScalarMatrix Q_gas = fvm::Sp(-(Kd + Ka) * alphaRho1, fluid_.phase1().Y(name_of_particles)) +
                            Ka * (alphaRho1 * fluid_.phase1().Y(name_of_particles) + alphaRho2 * fluid_.phase2().Y(name_of_particles));

    fvScalarMatrix Q_liquid = fvm::Sp(-(Kd + Ka) * alphaRho2, fluid_.phase2().Y(name_of_particles)) +
                            Kd * (alphaRho1 * fluid_.phase1().Y(name_of_particles) + alphaRho2 * fluid_.phase2().Y(name_of_particles));

    // Info << "KaarY2 : " << average(Ka * fluid_.phase2() * fluid_.phase2().rho() * fluid_.phase2().Y(name_of_particles)).value() << endl;                                
    // Info << "Ka : " << average(Ka).value() << nl << "Kd : " << average(Kd).value() << endl;

    return FlotationSoursePart(1e-3 * Q_gas, 1e-3 * Q_liquid);

}

volScalarField FlotationModel::gasSourseKa(){

    const volScalarField& epsilon = mesh_.lookupObject<volScalarField>("epsilon.liquid");
    const volScalarField& Es = stabilityEfficiency(epsilon);

    volScalarField Ka = K_attachment(epsilon, Es);

    return Ka * fluid_.phase2() * fluid_.phase2().rho() * fluid_.phase2().Y(name_of_particles);

}
