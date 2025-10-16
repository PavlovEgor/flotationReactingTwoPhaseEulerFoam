#ifndef FM
#define FM

#include "fvCFD.H"
#include "volFields.H"
#include "twoPhaseSystem.H"
#include "IOdictionary.H"

struct FlotationSoursePart
{
    fvScalarMatrix Q_phase1;
    fvScalarMatrix Q_phase2;

    FlotationSoursePart(fvScalarMatrix _Q_phase1, fvScalarMatrix _Q_phase2) :
        Q_phase1(_Q_phase1), Q_phase2(_Q_phase2) {}
    
    FlotationSoursePart(const FlotationSoursePart& FS) :
        Q_phase1(FS.Q_phase1), Q_phase2(FS.Q_phase2) {}

};

class FlotationModel
{
private:

    const fvMesh& mesh_;
    twoPhaseSystem& fluid_; 

    dimensionedScalar rhoP;
    dimensionedScalar d_p;
    dimensionedScalar d_b;
    dimensionedScalar nu_F;
    dimensionedScalar sigma;

    word name_of_particles;

    word attachmentModel;
    word inductionTimeModel;
    scalar A;
    scalar B;

    word collisionFrequencyModel;
    word collisionEfficiencyModel;

    word stabilityEfficiencyModel;
    word detachmentFrequencyModel;
    scalar C_1;
    scalar A_s;

    scalar theta_c_deg;
    scalar theta_c_rad;
    scalar N_max;


public:
    FlotationModel(const fvMesh& mesh, twoPhaseSystem& fluid, dimensionedScalar& delay);

    ~FlotationModel();

    volScalarField collisionEfficiency();

    volScalarField attachmentEfficiency();

    volScalarField stabilityEfficiency(const volScalarField& epsilon);

    volScalarField K_attachment(const volScalarField& epsilon, const volScalarField& Es);

    volScalarField K_detachment(const volScalarField& epsilon, const volScalarField& Es);

    FlotationSoursePart flotationMassTransfer();

    volScalarField gasSourseKa();
};


#endif
