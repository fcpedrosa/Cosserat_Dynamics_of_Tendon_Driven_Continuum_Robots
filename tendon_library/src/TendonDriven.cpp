#include "TendonDriven.hpp"

// overloaded constructor
template <std::size_t N, std::size_t numTendons>
TendonDriven<N, numTendons>::TendonDriven(double E, double G, double rad, double mass, double L, double alpha, double T, double dt, double tendonOffset, double tendonCompliance, double basePlateToMotor, double dampBendTwist, double dampShearExt, double dragCoeff) : m_E(E), m_G(G), m_rad(rad), m_mass(mass), m_L(L), m_alpha(alpha), m_T(T), m_dt(dt)
{
    this->m_g = { 0.0, 0.0, -9.81 };
	this->m_A = M_PI * pow(this->m_rad, 2); // cross-sectional area
    constexpr double I = M_PI_4 * pow(this->rad, 4);
    this->m_rho = this->m_mass / (this->m_L * this->m_A);
    blaze::diagonal(this->m_Kse) = {this->m_G * this->M_A, this->m_G * this->m_A, this->m_E * this->m_A};
    blaze::diagonal(this->m_Kbt) = { this->m_E * I, this->m_E * I, this->m_G * (I + I) };
	blaze::diagonal(this->m_J) = { I, I, 2 * I };
	this->m_BVPAccuracy = 1e-6;
	this->m_solType = cosseratEquations::STATIC_SOLUTION;
	this->m_method = rootFindingMethod::MODIFIED_NEWTON_RAPHSON;
	this->m_tau = 0.0;
	this->m_tendonCompliance = blaze::uniform(numTendons, tendonCompliance);
	blaze::diagonal(this->m_Bse) = blaze::uniform(3UL, dampShearExt);
	blaze::diagonal(this->m_Bbt) = blaze::uniform(3UL, dampBendTwist);
	blaze::diagonal(this->m_C) = blaze::uniform(3UL, dragCoeff);
	blaze::StaticVector<double, 3UL> v = { tendonOffset, 0.0, 0.0 };
	this->m_tendonOffset.fill(v);
	blaze::StaticVector<double, numTendons> angles;
	angles = blaze::linspace<blaze::rowVector>( numTendons, 0.0, 2 * M_PI );
	size_t idx = 0UL;
    for (auto d : angles)
    {
        this->m_r[idx] = { cos(angles[idx]), sin(angles[idx]), 0.0 };
        idx++; 
    }
	this->m_time = std::make_unique<timeManager>(this->m_alpha, this->m_dt);
	this->m_ODE = std::make_unique<ODESystem>();
	this->m_PDE = std::make_unique<PDESystem>();
}

// copy constructor
template <std::size_t N, std::size_t numTendons>
TendonDriven<N, numTendons>::TendonDriven(const TendonDriven &rhs) : m_E(rhs.m_E), m_G(rhs.m_G), m_rad(rhs.m_rad), m_mass(rhs.m_mass), m_L(rhs.m_L), m_alpha(rhs.m_alpha), m_T(rhs.m_T), m_dt(rhs.m_dt), m_tendonOffset(rhs.m_tendonOffset), m_tendonCompliance(rhs.m_tendonCompliance), m_basePlateToMotor(rhs.m_basePlateToMotor), m_A(rhs.m_A), m_rho(rhs.m_rho), m_Kse(rhs.m_Kse), m_Kbt(rhs.m_Kbt), m_J(rhs.m_J), m_BVPAccuracy(rhs.m_BVPAccuracy), m_solType(rhs.m_solType), m_method(rhs.m_method), m_tau(rhs.m_tau), m_Bse(rhs.m_Bse), m_Bbt(rhs.m_Bbt), m_C(rhs.m_C), m_r(rhs.m_r), m_g(rhs.m_g) 
{
	this->m_time = std::make_unique<timeManager>(this->m_alpha, this->m_dt);
	this->m_ODE = std::make_unique<ODESystem>();
	this->m_PDE = std::make_unique<PDESystem>();
}

// move constructor
template <std::size_t N, std::size_t numTendons>
TendonDriven<N, numTendons>::TendonDriven(TendonDriven &&rhs) noexcept
{
    // handling self-assignment
    if (this != &rhs)
    {
        this->m_E = rhs.m_E;
		this->m_G = rhs.m_G;
		this->m_rad = rhs.m_rad;
		this->m_mass = rhs.m_mass;
		this->m_L = rhs.m_L;
		this->m_alpha = rhs.m_alpha;
		this->m_T = rhs.m_T;
		this->m_dt = rhs.m_dt;
		this->m_tendonOffset = rhs.m_tendonOffset;
		this->m_tendonSlack = std::move(rhs.m_tendonSlack);
		this->m_tendonCompliance = std::move(rhs.m_tendonCompliance);
		this->m_basePlateToMotor = rhs.m_basePlateToMotor;
		this->m_A = rhs.m_A;
		this->m_rho = rhs.m_rho;
		this->m_Kse = std::move(rhs.m_Kse);
		this->m_Kbt = std::move(rhs.m_Kbt);
		this->m_J = std::move(rhs.m_J);
		this->m_BVPAccuracy = rhs.m_BVPAccuracy;
		this->m_solType = rhs.m_solType;
		this->m_method = rhs.m_method;
		this->m_tau = std::move(rhs.m_tau);
		this->m_Bse = std::move(rhs.m_Bse);
		this->m_Bbt = std::move(rhs.m_Bbt);
		this->m_C = std::move(rhs.m_C);
		this->m_r = std::move(rhs.m_r);
		this->m_g = std::move(rhs.m_g);
		this->m_y = std::move(rhs.m_y);
		this->m_z = std::move(rhs.m_z);
		this->m_zh = std::move(rhs.m_zh);
		this->m_time = std::move(rhs.m_time);
		this->m_ODE = std::move(rhs.m_ODE);
		this->m_PDE = std::move(rhs.m_PDE);
    }
}

// CTR destructor
template <std::size_t N, std::size_t numTendons>
TendonDriven<N, numTendons>::~TendonDriven()
{
	// nothing to be done
}

// copy assignment operator
template <std::size_t N, std::size_t numTendons>
TendonDriven<N, numTendons> &TendonDriven<N, numTendons>::operator=(const TendonDriven<N, numTendons> &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_E = rhs.m_E;
		this->m_G = rhs.m_G;
		this->m_rad = rhs.m_rad;
		this->m_mass = rhs.m_mass;
		this->m_L = rhs.m_L;
		this->m_alpha = rhs.m_alpha;
		this->m_T = rhs.m_T;
		this->m_dt = rhs.m_dt;
		this->m_tendonOffset = rhs.m_tendonOffset;
		this->m_tendonSlack = rhs.m_tendonSlack;
		this->m_tendonCompliance = rhs.m_tendonCompliance;
		this->m_basePlateToMotor = rhs.m_basePlateToMotor;
		this->m_A = rhs.m_A;
		this->m_rho = rhs.m_rho;
		this->m_Kse = rhs.m_Kse;
		this->m_Kbt = rhs.m_Kbt;
		this->m_J = rhs.m_J;
		this->m_BVPAccuracy = rhs.m_BVPAccuracy;
		this->m_solType = rhs.m_solType;
		this->m_method = rhs.m_method;
		this->m_tau = rhs.m_tau;
		this->m_Bse = rhs.m_Bse;
		this->m_Bbt = rhs.m_Bbt;
		this->m_C = rhs.m_C;
		this->m_r = rhs.m_r;
		this->m_g = rhs.m_g;
		this->m_y = rhs.m_y;
		this->m_z = rhs.m_z;
		this->m_zh = rhs.m_zh;
		this->m_time = rhs.m_time;
		this->m_ODE = rhs.m_ODE;
		this->m_PDE = rhs.m_PDE;
	}

	return *this;
}

// move assignment operator
template <std::size_t N, std::size_t numTendons>
TendonDriven<N, numTendons> &TendonDriven<N, numTendons>::operator=(TendonDriven<N, numTendons> &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_E = rhs.m_E;
		this->m_G = rhs.m_G;
		this->m_rad = rhs.m_rad;
		this->m_mass = rhs.m_mass;
		this->m_L = rhs.m_L;
		this->m_alpha = rhs.m_alpha;
		this->m_T = rhs.m_T;
		this->m_dt = rhs.m_dt;
		this->m_tendonOffset = rhs.m_tendonOffset;
		this->m_tendonSlack = std::move(rhs.m_tendonSlack);
		this->m_tendonCompliance = std::move(rhs.m_tendonCompliance);
		this->m_basePlateToMotor = rhs.m_basePlateToMotor;
		this->m_A = rhs.m_A;
		this->m_rho = rhs.m_rho;
		this->m_Kse = std::move(rhs.m_Kse);
		this->m_Kbt = std::move(rhs.m_Kbt);
		this->m_J = std::move(rhs.m_J);
		this->m_BVPAccuracy = rhs.m_BVPAccuracy;
		this->m_solType = rhs.m_solType;
		this->m_method = rhs.m_method;
		this->m_tau = std::move(rhs.m_tau);
		this->m_Bse = std::move(rhs.m_Bse);
		this->m_Bbt = std::move(rhs.m_Bbt);
		this->m_C = std::move(rhs.m_C);
		this->m_r = std::move(rhs.m_r);
		this->m_g = std::move(rhs.m_g);
		this->m_y = std::move(rhs.m_y);
		this->m_z = std::move(rhs.m_z);
		this->m_zh = std::move(rhs.m_zh);
		this->m_time = std::move(rhs.m_time);
		this->m_ODE = std::move(rhs.m_ODE);
		this->m_PDE = std::move(rhs.m_PDE);
	}

	return *this;
}

// returns the residue (violation of boundary condition) for a particular initial guess
template <std::size_t N, std::size_t numTendons>
blaze::StaticVector<double, 6UL + numTendons> TendonDriven<N, numTendons>::residueFunction(const blaze::StaticVector<double, 6UL + numTendons> &initGuess)
{	
	// initGuess := [internal force n(0), internal moment m(0), tendon_lengths]
    blaze::StaticVector<double, 3UL> v0, u0, e3;
    e3 = {0.0, 0.0, 1.0};
    v0 = e3 + blaze::inv(this->m_Kse) * blaze::subvector<0UL, 3UL>(initGuess); // v = v* + K^(-1)*R^T*n
    u0 = blaze::inv(this->m_Kbt) * *blaze::subvector<3UL, 3UL>(initGuess);     // u = u* + K^(-1)*R^T*m

    // initial conditions for numerical integration
    blaze::StaticVector<double, 19UL + numTendons> y0, dyds;
    blaze::subvector<0UL, 19UL>(y0) = {
        0.0, 0.0, 0.0,          // initial position
        1.0, 0.0, 0.0, 0.0,     // initial orientation (the 1-0-0-0 quaternion yields the identity matrix)
        v[0UL], v[1UL], v[2UL], // initial strain
        u[0UL], u[1UL], u[2UL], // initial curvature
        0.0, 0.0, 0.0,          // initial linear velocity q_s
        0.0, 0.0, 0.0           // initial angular velocity w_s
    }                           
    blaze::subvector<20UL, numTendons>(y0) = blaze::uniform(numTendons, this->m_basePlateToMotor); // initial tendon length
    blaze::column<0UL>(this->m_y) = y0;

    // initGuess := [internal force n(0), internal moment m(0), tendon_lengths]
    this->m_tau = blaze::max(blaze::subvector<6UL, numTendons>(initGuess), 0.0);          // tendon-displacement control => tau >= 0.0
    this->m_tendonSlack = -blaze::min(blaze::subvector<6UL, numTendons>(initGuess), 0.0); // tension and slack are mutually exclusive => slack < 0

    // Starting the numerical integration of the model equations -- Must define if Static or Dynamic solutions are to be used
    // numerical integration of the model equations (handles static or dynamic solutions)
	switch (this->m_solType)
	{
		case cosseratEquations::STATIC_SOLUTION:
			this->m_ODE->euler_ODE(y0, m_y, dyds, m_z, m_L, this);
			break;
		case cosseratEquations::DYNAMIC_SOLUTION:
			this->m_PDE->euler_PDE(y0, m_y, dyds, m_zh, m_z, m_L, this);
			break;
	}

    // variables for computing the internal loading at the distal end (prior to final plate)
    blaze::StaticVector<double, 3UL> vL, uL, vL_t, uL_t;
    vL = blaze::subvector<7UL, 3UL>(blaze::column<N - 1>(this->m_y));
    uL = blaze::subvector<10UL, 3UL>(blaze::column<N - 1>(this->m_y));

    switch (this->m_solType)
    {
		case cosseratEquations::STATIC_SOLUTION:
			uL_t = 0.0;
			vL_t = 0.0;
			break;
		case cosseratEquations::DYNAMIC_SOLUTION:
			double c0 = this->m_time->getC0();
			vL_t = c0 * vL + blaze::subvector<0UL, 3UL>(blaze::column<N - 1>(this->z_h));
			uL_t = c0 * uL + blaze::subvector<3UL, 3UL>(blaze::column<N - 1>(this->z_h));
			break;
    }

    blaze::StaticVector<double, 3UL> nb, mb;
    nb = this->K_se * (vL - e3) + this->Bse * vL_t; // internal force at distal end (prior to final plate)
    mb = this->m_Kbt * uL + this->m_Bbt * uL_t;     // internal moment at distal end (prior to final plate)

    // computing the equilibrium error
    blaze::StaticVector<double, 3UL> forceError(-nb), momentError(-mb), pbs_i, Fb_i;
    for (size_t idx = 0; idx < numTendons; ++idx)
    {
        pbs_i = blaze::cross(uL, this->m_r[idx]) + vL;
        Fb_i = -this->m_tau[idx]) * blaze::normalize(pbs_i);
        forceError += Fb_i;
        momentError += blaze::cross(this->m_r[idx], Fb_i);
    }

    // compute the violation in tendon Lengths
    blaze::StaticVector<double, numTendons> length, lengthStar, lengthError, stretch, q_t;
    length = blaze::subvector<20UL, numTendons>(blaze::column<N - 1>(this->m_y)); // integrated tendon lengths

    // WILL APPLY THE DYNAMIC ACTUATION TO THE FIRST TENDON ONLY
    // -->> (THINK ABOUT A MORE FLEXIBLE WAY TO IMPLEMENT THIS LATER)
    q_t[0UL] = this->m_time->getDynamicActuationProfile();
    // The reference tendon length is l∗(t) = L + basePlateToMotor + q(t)
    lengthStar = blaze::uniform(numTendons, this->m_L + this->m_basePlateToMotor) + q_t; // q(t) is being applied to tendon 1 only (as of now)
    stretch = this->m_tendonCompliance * lengthStar * this->m_tau;
    lengthError = length + this->m_tendonSlack - (lengthStar + stretch);

    // computing the residue at the distal end
    blaze::StaticVector<double, 6UL + numTendons> residue;
    blaze::subvector<0UL, 3UL>(residue) = std::move(forceError);
    blaze::subvector<3UL, 3UL>(residue) = std::move(momentError);
    blaze::subvector<6UL, numTendons>(residue) = std::move(lengthError);

    return residue;
}

// function that computes the finite-differences Jacobian for solving the BVP
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 6UL + numTendons, 6UL + numTendons> TendonDriven<N, numTendons>::jac_BVP(const blaze::StaticVector<double, 6UL + numTendons> &initGuess, const blaze::StaticVector<double, 6UL + numTendons> &residue)
{
	const size_t n = 6UL + numTendons;
	blaze::StaticMatrix<double, n, n, blaze::columnMajor> Jac_bvp;
	blaze::StaticVector<double, n> initGuessPerturbed(initGuess), residuePerturbed; 
	double eps = 1e-4;

	for (size_t iter = 0UL; iter < n; ++iter)
	{
		initGuessPerturbed[iter] += eps;
		// perturbed residue
		residuePerturbed = this->residueFunction(initGuessPerturbed);
		// building the finite-differences Residue jacobian
		blaze::column(jac_bvp, iter) = (residuePerturbed - residue) / eps;
		// restoring the original value of the array
		initGuessPerturbed[iter] = initGuess[iter];
	}

	return jac_bvp;
}

// function that computes the finite-differences Jacobian wrt actuation inputs
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 6UL> TendonDriven<N, numTendons>::jacobian(const blaze::StaticVector<double, 6UL + numTendons> &initGuess, const blaze::StaticVector<double, 3UL> &tipPos)
{
}

// manages which method to be used for solving the Boundary-value problem
template <std::size_t N, std::size_t numTendons>
bool TendonDriven<N, numTendons>::solveBVP(blaze::StaticVector<double, 6UL + numTendons> &initGuess)
{
	bool convergence; // flags the convergence (or lack thereof) of the selected nonlinear root-finding method

	switch (this->m_method)
	{
		case rootFindingMethod::POWELL_DOG_LEG:
			convergence = this->PowellDogLeg(initGuess);
			break;
		case rootFindingMethod::LEVENBERG_MARQUARDT:
			convergence = this->Levenberg_Marquardt(initGuess);
			break;
		case rootFindingMethod::BROYDEN:
			convergence = this->Broyden(initGuess);
			break;
		case rootFindingMethod::NEWTON_RAPHSON:
			convergence = this->Newton_Raphson(initGuess);
			break;
		case rootFindingMethod::MODIFIED_NEWTON_RAPHSON:
			convergence = this->Modified_Newton_Raphson(initGuess);
			break;
	}

	return convergrnce;
}

// function that implements Powell's Dog Leg Method (Nonlinear root-finding method for solving the BVP)
template <std::size_t N, std::size_t numTendons>
bool TendonDriven<N, numTendons>::PowellDogLeg(blaze::StaticVector<double, 6UL + numTendons> &initGuess)
{
    bool found;
    size_t k = 0UL;
    const size_t k_max = 200UL;
    double alpha, beta, delta, eps1, eps2, rho, c;
    blaze::StaticVector<double, 6UL + numTendons> g, f, f_new, x_new, h_sd, h_gn, h_dl;
    blaze::StaticMatrix<double, 6UL + numTendons, 6UL + numTendons> J;

    // initializing parameters
    delta = 1;
    eps1 = eps2 = 1e-22;

    f = this->residueFunction(initGuess);
    J = this->jac_BVP(initGuess, f);
    g = blaze::trans(J) * f;

    // checking if the initial guess satisfies the BVP without the need of any further refinement
    found = ((blaze::linfNorm(f) <= this->m_BVPAccuracy) || (blaze::linfNorm(g) <= eps1)) ? true : false;

    while (!found && (k < k_max))
    {
        k++;

        alpha = blaze::sqrNorm(g) / blaze::sqrNorm(J * g);
        h_sd = -alpha * g;         // steepest descend (this is a direction, not a step!)
        h_gn = -this->pInv(J) * f; // Gauss-Newton step (Least Square solution)

        // two candidates for the step to take from this point, a = alpha*h_sd & b = h_gn

        // computing the dog leg direction
        if (blaze::norm(h_gn) <= delta)
            h_dl = h_gn;
        else
        {
            if (blaze::norm(h_sd) >= delta)
                h_dl = delta * blaze::normalize(h_sd);
            else
            {
                c = blaze::trans(h_sd) * (h_gn - h_sd);

                if (c <= 0.0)
                {
                    beta = (-c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd)))) / blaze::sqrNorm(h_gn - h_sd);
                }
                else
                {
                    beta = (delta * delta - blaze::sqrNorm(h_sd)) / (c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd))));
                }

                h_dl = h_sd + beta * (h_gn - h_sd); // Dog Leg step
            }
        }

        if (blaze::norm(h_dl) <= eps2 * (blaze::norm(initGuess) + eps2))
            found = true;
        else
        {
            x_new = initGuess + h_dl;
            f_new = this->residueFunction(x_new);
            rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h_dl) * ((delta * h_dl) - g));

            if (rho > 0.0)
            {
                initGuess = std::move(x_new);
                f = std::move(f_new);
                J = this->jac_BVP(initGuess, f);
                g = blaze::trans(J) * f;

                if ((blaze::linfNorm(f) <= this->m_BVPAccuracy) || (blaze::linfNorm(g) <= eps1))
                    found = true;
            }

            if (rho > 0.75)
                delta = std::max(delta, 3 * blaze::norm(h_dl));
            else
            {
                if (rho < 0.25)
                    delta *= 0.5;
            }

            if (delta < eps2 * (blaze::norm(initGuess) + eps2))
                found = true;
        }
    }

    if (!found)
        std::cout << "Powell's Dog Leg failed! --> Residue = " << std::setprecision(12UL) << blaze::trans(f);

    return found;
}

// function that implements the Levenberg-Marquardt Method (Nonlinear root-finding method for solving the BVP)
template <std::size_t N, std::size_t numTendons>
bool TendonDriven<N, numTendons>::Levenberg_Marquardt(blaze::StaticVector<double, 6UL + numTendons> &initGuess)
{
    	size_t k = 0UL;
	const size_t k_max = 200UL;
	blaze::StaticVector<double, 6UL + numTendons> h, g, f, f_new;
	blaze::StaticMatrix<double, 6UL + numTendons, 6UL + numTendons> J, A;
	blaze::IdentityMatrix<double> I(6UL + numTendons);
	double rho, nu = 2, mu, tau = 1e-3, e1 = 1e-18, e2 = 1e-25;
	bool found;

	// computing the residue and residue Jacobian associated to initGuess
	f = this->residueFunction(initGuess);
	J = this->jac_BVP(initGuess, f);
	A = blaze::trans(J) * J;
	g = blaze::trans(J) * f;
	found = (blaze::linfNorm(g) <= e1) ? true : false;
	mu = tau * blaze::max(blaze::diagonal(A));

	// starting the iterative minimization loop
	while ((!found) && (k < k_max))
	{
		k++;
		blaze::solve(blaze::declsym(A + (mu * I)), h, -g);

		f_new = this->residueFunction(initGuess + h);
		rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h) * ((mu * h) - g));

		if (rho > 0.0)
		{
			// accept the decrease in the function
			initGuess += h;
			// computing the residue Jacobian at the new initial guess
			J = this->jac_BVP(initGuess, f_new);
			A = blaze::trans(J) * J;
			f = std::move(f_new);
			g = blaze::trans(J) * f;
			found = (blaze::linfNorm(g) <= e1) ? true : false;
			mu = mu * std::max(0.33333333, 1 - blaze::pow(2 * rho - 1, 3));
			nu = 2;
		}
		else
		{
			mu = mu * nu;
			nu = 2 * nu;
		}

		// checking if the tolerance has been satisfied
		if (blaze::linfNorm(f) <= this->m_BVPAccuracy)
			found = true;
	}

	if (!found)
		std::cout << "Levenberg-Marquardt failed! --> Residue = " << std::setprecision(12UL) << blaze::trans(f);

	return found;
}

// function that implements Broyden's Nonlinear root-finding method for solving the BVP (Jacobian inverse is estimated)
template <std::size_t N, std::size_t numTendons>
bool TendonDriven<N, numTendons>::Broyden(blaze::StaticVector<double, 6UL + numTendons> &initGuess)
{
    // found: returns true (false) when the root-finding method converges (does not converge) within k_max iterations
	bool found;

	// initial Hessian matrix --> computed via finite differences
	blaze::StaticMatrix<double, 6UL + numTendons, 6UL + numTendons> JacInv, JacInvNew;

	// setting up and starting my handmadeBFGS method
	blaze::StaticVector<double, 6UL + numTendons> F, Fold, X, Xold, deltaX, deltaF; // staticVectors are automatically initialized to 0

	// Residue yielded by the initial guess for the CTR BVP
	F = this->residueFunction(initGuess); // F(x_k)	: residue
	X = std::move(initGuess);		      // x_k		: initial guess
	JacInvNew = JacInv = this->pInv(this->jac_BVP(X, F));

	// checking if the initial guess already satisfies the BVP
	found = (blaze::linfNorm(F) <= this->m_BVPAccuracy) ? true : false;

	size_t k = 0UL;
	const size_t k_max = 200UL;
	while (!found && (k < k_max))
	{
		k++;

		deltaX = X - Xold; // dX := x_k - x_k-1
		deltaF = F - Fold; // dF := F(x_k) - F(x_k-1)

		JacInv = std::move(JacInvNew);
		if ((blaze::norm(deltaX) > 0.0) && (blaze::norm(deltaF) > 0.0))
			JacInvNew = JacInv + ((deltaX - JacInv * deltaF) / (blaze::trans(deltaX) * JacInv * deltaF)) * blaze::trans(deltaX) * JacInv;
		else
			JacInvNew = JacInv;

		Xold = std::move(X);
		Fold = std::move(F);

		// update the initial guess
		X = Xold - JacInv * F;
		F = this->residueFunction(X);

		while (blaze::isnan(F))
		{
			if (blaze::isnan(X))
				X = 0.0;
			else
				X /= blaze::max(blaze::abs(X));

			F = this->residueFunction(X);
			JacInv = JacInvNew = this->pInv(this->jac_BVP(X, F));
			Xold = std::move(X);
			X = Xold - JacInv * F;
		}

		if (k % 10 == 0.0)
		{
			JacInv = JacInvNew = this->pInv(this->jac_BVP(X, F));
			X = Xold - JacInv * F;
		}

		if (blaze::linfNorm(F) <= this->m_BVPAccuracy)
			found = true;
	}

	// std::cout << "Broyden finalized in " << k << " iterations. | Residue = " << std::setprecision(12UL) << blaze::trans(F);
	
	initGuess = std::move(X);
	return found;
}

// function that implements the Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
template <std::size_t N, std::size_t numTendons>
bool TendonDriven<N, numTendons>::Newton_Raphson(blaze::StaticVector<double, 6UL + numTendons> &initGuess)
{
    	bool found;
	// setting up and starting my handmade Newton-Raphson method
	blaze::StaticVector<double, 6UL + numTendons> Residue, Residue_new, d_Residue, dGuess; // staticVectors are automatically initialized to 0

	// Residue of the unperturbed initial guess for the CTR
	Residue = this->residueFunction(initGuess);

	found = (blaze::linfNorm(Residue) <= this->m_BVPAccuracy)? true : false;

	//  Weighing matrices for adjusting the initial guess iteratively (Implementing a PD regulator)
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 6UL + numTendons, 6UL + numTendons, blaze::rowMajor>> Kp, Kd;
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> jac_bvp;
	blaze::diagonal(Kp) = 0.5;	// 0.45 | 0.6  | 0.3
	blaze::diagonal(Kd) = 3e-3; // 3e-3 | 5e-3 | 2e-3

	size_t k = 0UL;
	const size_t k_max = 200UL;
	// starting iterations for adjusting the initial guess "u_guess ~ initGuess"
	while (!found && (k < k_max))
	{
		k++;
		jac_bvp = this->jac_BVP(initGuess, Residue);
		// error equation(globally asymptotically stable)
		dGuess = this->pInv(jac_bvp) * (Kp * Residue + Kd * d_Residue);
		// updating the initial guess(weighted negative gradient of the cost function)
		initGuess -= dGuess;
		// computing the new cost associated to the newly readjusted initial guess
		Residue_new = this->residueFunction(initGuess);

		// Checking if the Jacobian has large elements beyond machine precision
		while (blaze::isnan(Residue_new))
		{
			if (blaze::isnan(initGuess))
				initGuess = 0.0;
			else
				initGuess /= blaze::max(blaze::abs(initGuess));
			
			Residue_new = this->residueFunction(initGuess);
			d_Residue = Residue_new - Residue;
			Residue = std::move(Residue_new);
			continue;
		}

		// cost variation due to initial guess refinement
		d_Residue = Residue_new - Residue;
		// updating the cost
		Residue = std::move(Residue_new);		

		if (blaze::linfNorm(Residue) <= this->m_BVPAccuracy)
			found = true;
	}

	if (!found)
		std::cout << "Newton-Raphson failed! --> Residue = " << std::setprecision(12UL) << blaze::trans(Residue);

	return found;
}

// function that implements the Modified, globally convergent Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
template <std::size_t N, std::size_t numTendons>
bool TendonDriven<N, numTendons>::Modified_Newton_Raphson(blaze::StaticVector<double, 6UL + numTendons> &initGuess)
{
    	/*
		Algorithm extracted from page 309 of Introduction to Numerical Analysis 3rd edition by Josef Stoer & Roland Bulirsch
	*/

	bool found;
	// computes the residue associated to the initial guess
	blaze::StaticVector<double, 6UL + numTendons> f(this->residueFunction(initGuess)), d;
	blaze::StaticVector<double, 6UL + numTendons, blaze::rowVector> Dh;
	blaze::StaticMatrix<double, 6UL + numTendons, 6UL + numTendons> D, D_inv;
	double h, h_0, lambda, gamma, improvementFactor, d_norm, Dh_norm;
	size_t j = 0UL, k = 0UL;
	const size_t k_max = 250UL;
	std::vector<double> h_k; // vector to store all h_k's
	h_k.reserve(k_max);
	std::vector<double>::iterator result; // iterator for determining the min element in the vector h_k

	found = (blaze::linfNorm(f) <= this->m_BVPAccuracy) ? true : false;

	while (!found && (k < k_max))
	{
		k++;
		// computing the residue Jacobian
		D = this->jac_BVP(initGuess, f);
		D_inv = this->pInv(D);

		// search direction (directional derivative)
		d = D_inv * f;
		gamma = 1 / (blaze::norm(D_inv) * blaze::norm(D)); // gamma := 1/cond(Df)
		h_0 = blaze::sqrNorm(f);						   // h := f'f
		// Dh := D(f'f) = 2f'Df
		Dh = 2 * blaze::trans(f) * D;
		d_norm = blaze::norm(d);
		Dh_norm = blaze::norm(Dh);

		while (true)
		{
			f = this->residueFunction(initGuess - blaze::pow(0.5, j) * d);
			// std::cout << "Modified_Newton_Raphson -- j = : " << j << " | residue = " << blaze::trans(f);
			while (blaze::isnan(f))
			{
				j++;
				f = this->residueFunction(initGuess - blaze::pow(0.5, j) * d);
			}
			h = blaze::sqrNorm(f);
			improvementFactor = blaze::pow(0.5, j) * 0.25 * gamma * d_norm * Dh_norm;
			// storig the value of h_k to determine step size posteriorly
			h_k.push_back(h);

			if (h <= (h_0 - improvementFactor))
				break;
			else
				j++;
		}

		// determining the value of the step-size lambda
		result = std::min_element(std::execution::par, h_k.begin(), h_k.end());
		// retrieving the minimum h_k
		lambda = blaze::pow(0.5, std::distance(h_k.begin(), result));
		initGuess -= lambda * d;
		h_k.clear();

		// resets the exponent variable j
		j = 0UL;

		// compute the residue associated to the newly refined initGuess
		f = this->residueFunction(initGuess);

		// checking the terminating condition
		if (blaze::linfNorm(f) <= this->m_BVPAccuracy)
		{
			found = true;
		}
		else
		{
			if (blaze::isnan(f))
			{
				if (blaze::isnan(initGuess))
					initGuess = 0.0;
				else
					initGuess /= blaze::max(blaze::abs(initGuess));

				// recompute the residue of the readjusted initial guess
				f = this->residueFunction(initGuess);
			}
		}
	}
	
	return found;
}

// method for computing the Penrose Pseudoinverse via SVD decomposition
template <std::size_t N, std::size_t numTendons>
blaze::HybridMatrix<double, 6UL + numTendons, 6UL + numTendons> TendonDriven<N, numTendons>::pInv(const blaze::HybridMatrix<double, 6UL + numTendons, 6UL + numTendons> &M)
{
    // declaring the auxiliary matrices for pInv computation
    blaze::HybridMatrix<double, 6UL + numTendons, 6UL + numTendons> U; // The matrix for the left singular vectors
    blaze::HybridVector<double, 6UL + numTendons> s;                   // The vector for the singular values
    blaze::HybridMatrix<double, 6UL + numTendons, 6UL + numTendons> V; // The matrix for the right singular vectors

    // SVD decomposition of the matrix M
    try
    {
        blaze::svd(M, U, s, V);
    }
    catch (std::exception &e)
    {
        std::cerr << "Blaze SVD has failed: " << e.what() << std::endl;
        // std::cout << "CTR joints: " << blaze::trans(this->m_q) << std::endl;
        std::cerr << M << std::endl;
    }

    // computing the pseudoinverse of M via SVD decomposition
    blaze::DiagonalMatrix<blaze::HybridMatrix<double, 6UL + numTendons, 6UL + numTendons>> S_inv;
    // resizing S - as per the number of nonzero singular values
    const size_t len_S = s.size();
    S_inv.resize(len_S, len_S);

    // Creating a reference to the diagonal of matrix S
    auto diag = blaze::diagonal(S_inv);
    // applies a "damping factor" to the zero singular values of the matrix M
    diag = blaze::map(s, [](double d)
                      { return d / ((d * d) + 1e-25); }); // Damped least squares -- SVD pseudo inverse

    blaze::HybridMatrix<double, 6UL + numTendons, 6UL + numTendons> res = blaze::trans(U * S_inv * V);

    return res;
}