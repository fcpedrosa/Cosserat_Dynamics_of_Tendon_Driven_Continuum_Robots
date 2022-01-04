#pragma once

#include <blaze/Math.h>

template<std::size_t N, std::size_t numTendons>
class ODESystem {
    public:

    // default constructor
	ODESystem();

	// overloaded constructor
	ODESystem(blaze::StaticVector<double, 3UL> u, blaze::StaticVector<double, 3UL> v, blaze::StaticVector<double, 4UL> h);

	// copy constructor
	ODESystem(const ODESystem& rhs);

	// move constructor
	ODESystem(ODESystem&& rhs) noexcept;

	// CTR destructor
	~ODESystem();

	// copy assignment operator
	ODESystem& operator=(const ODESystem& rhs);

	// move assignment operator
	ODESystem& operator=(ODESystem&& rhs) noexcept;

    // functor that implements the system of ODEs governing an N-tendon catheter
    void operator()(const blaze::StaticVector<double, 19UL + numTendons>& y, blaze::StaticVector<double, 19UL + numTendons>& dyds, blaze::StaticVector<double, 12UL>& dydt, const tendonDriven &robot);

	// Integrates the ODE model equations using Euler's method
	void euler_ODE(const blaze::StaticVector<double, 19UL + numTendons>& y0, blaze::StaticMatrix<double, 19UL + numTendons, N>& Y, blaze::StaticVector<double, 19UL + numTendons>& dyds, blaze::StaticMatrix<double, 18UL + numTendons, N>& Z, const double L, tendonDriven &robot);

	// Integrates the ODE model equations using the classic 4th-order Runge-Kutta algorithm
	void rungeKutta_ODE(const blaze::StaticVector<double, 19UL + numTendons>& y0, blaze::StaticMatrix<double, 19UL + numTendons, N>& Y, blaze::StaticMatrix<double, 18UL + numTendons, N>& Z, tendonDriven &robot);

    // function that returns a rotation matrix in SO(3) from a set of non-unity quaternions
    blaze::StaticMatrix<double, 3UL, 3UL> getSO3(const blaze::StaticVector<double, 4UL>& h);

	// function that computes the squared hat operator
	blaze::StaticMatrix<double, 3UL, 3UL> hatSqr(const blaze::StaticVector<double, 3UL>& v);

	// function that computes the posmultiplication of a vector by the hat operator
	blaze::StaticMatrix<double, 3UL, 3UL> hatPostultiply(const blaze::StaticVector<double, 3UL>& v1, const blaze::StaticVector<double, 3UL>& v2);

    private:
        blaze::StaticMatrix<double, 3UL, 3UL> m_R;    // rotation matrix SO(3) -- orientation of local frame wrt global frame
        blaze::StaticVector<double, 4UL> m_h;         // non-unity quaternions -- parametrizing the rotation
        blaze::StaticVector<double, 3UL> m_u, m_v;    // u: curvature vector | v: linear rate of change of position wrt arc-length
		blaze::StaticMatrix<double, 6UL, 6UL> m_Phi;  // linear system matrix for determining vs, us
};