

const volScalarField& my_epsilon = mesh.lookupObject<volScalarField>("epsilon.liquid");

volScalarField theta_a = 2 * Foam::atan(Foam::exp(-3 * Foam::mag(U1) * t_i / (d_b * ((d_b / d_p) + 1))));
volScalarField E_a = Foam::pow(Foam::sin(theta_a), 2);

volScalarField down = (d_p * d_p) * ((rhoP - phase2.rho()) * Foam::mag(g) + rhoP * (1.9 * Foam::pow(my_epsilon, 2.0/3.0) * Foam::pow((d_p + d_b)/2, -1.0/3.0)))
    + (3.0 / 2.0) * d_p * (Foam::sin(M_PI - theta_c/2) * Foam::sin(M_PI - theta_c/2))
    * ((4 * sigma / d_b) - d_b * phase2.rho() * Foam::mag(g));

// volScalarField down = (d_p * d_p) * ((rhoP - phase2.rho()) * Foam::mag(g) + rhoP * (1.9 * Foam::pow(my_epsilon, 2.0/3.0) * Foam::pow((d_p + d_b)/2, -1.0/3.0)))
//     + (3.0 / 2.0) * d_p * (Foam::sin(theta_c/2) * Foam::sin(theta_c/2))
//     * ((4 * sigma / d_b) - d_b * phase2.rho() * Foam::mag(g));

oneDivBo_prime = Foam::mag(6 * sigma * sin(M_PI + theta_c/2) * sin(M_PI + theta_c/2)) / down;

E_s = 1 - Foam::exp(A_s * (1 - (oneDivBo_prime * pos(oneDivBo_prime - 1) + 1 * (1 - pos(oneDivBo_prime - 1)))));

// dimensionedScalar E_s = 1 - Foam::exp(A_s * (1 - 2)) ;




volScalarField k1 = Z_0 * Foam::sqrt(my_epsilon) * E_c * E_a * (1 - E_s);
volScalarField k2 = Z_prime_0 * Foam::cbrt(my_epsilon) * E_s;
volScalarField beta = (rho1/ rho2) * (phase1.Y("particle") / Nmax) * (d_b * d_b * d_b / (d_p * d_p * d_p));


KA = k1 * (1 - beta) * alpha1 / (( M_PI / 6.0) * (d_b * d_b * d_b));
KD = k2 / Nmax;

Info << "my_epsilon: " << average(E_s) << nl << endl;
// Info << "E_a: " << E_a << nl << endl;



