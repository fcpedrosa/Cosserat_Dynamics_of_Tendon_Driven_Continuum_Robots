#pragma once

#include <blaze/Math.h>

template <std::size_t N, std::size_t numTendons>
class PDESystem
{
public:
    // default constructor
    PDESystem();

    // overloaded constructor
    PDESystem(blaze::StaticVector<double, 3UL> u, blaze::StaticVector<double, 3UL> v, blaze::StaticVector<double, 4UL> h);

    // copy constructor
    PDESystem(const PDESystem &rhs);

    // move constructor
    ODESystem(ODESystem &&rhs) noexcept;

    // CTR destructor
    ~PDESystem();

    // copy assignment operator
    PDESystem &operator=(const PDESystem &rhs);

    // move assignment operator
    PDESystem &operator=(PDESystem &&rhs) noexcept;

    // functor that implements the system of ODEs governing an N-tendon catheter
    void operator()(const blaze::StaticVector<double, 19UL + numTendons> &y, blaze::StaticVector<double, 19UL + numTendons> &dyds, const blaze::StaticVector<double, 18UL> &z_h, blaze::StaticVector<double, 18UL> &z, const tendonDriven &robot);

    // Integrates the ODE model equations using Euler's method
    void euler_PDE(const blaze::StaticVector<double, 19UL + numTendons> &y0, double, 19UL + numTendons> &y, blaze::StaticVector<double, 19UL + numTendons> &dyds, const blaze::StaticVector<double, 18UL> &z_h, blaze::StaticVector<double, 18UL> &z, const tendonDriven &robot);

    // Integrates the ODE model equations using the classic 4th-order Runge-Kutta algorithm
    void rungeKutta_PDE(const blaze::StaticVector<double, 19UL + numTendons> &y0, double, 19UL + numTendons> &y, blaze::StaticVector<double, 19UL + numTendons> &dyds, const blaze::StaticVector<double, 18UL> &z_h, blaze::StaticVector<double, 18UL> &z, const tendonDriven &robot);

    // function that returns a rotation matrix in SO(3) from a set of non-unity quaternions
    blaze::StaticMatrix<double, 3UL, 3UL> getSO3(const blaze::StaticVector<double, 4UL> &h);

    // function that computes the squared hat operator
    blaze::StaticMatrix<double, 3UL, 3UL> hatSqr(const blaze::StaticVector<double, 3UL> &v);

    // function that computes the premultiplication of a vector by the hat operator
    blaze::StaticMatrix<double, 3UL, 3UL> PDESystem::hatPreMultiply(const blaze::StaticVector<double, 3UL> &v, const blaze::StaticMatrix<double, 3UL 3UL> &M);

    // function that computes the posmultiplication of a vector by the hat operator
    blaze::StaticMatrix<double, 3UL, 3UL> PDESystem::hatPostMultiply(const blaze::StaticMatrix<double, 3UL 3UL> &M, const blaze::StaticVector<double, 3UL> &v);

    // function that computes the differential quaternion evolution
    blaze::StaticVector<double, 4UL> PDESystem::quaternionDiff(const blaze::StaticVector<double, 3UL> &u, const blaze::StaticVector<double, 4UL> &h);

private:
    blaze::StaticMatrix<double, 3UL, 3UL> m_R;   // rotation matrix SO(3) -- orientation of local frame wrt global frame
    blaze::StaticVector<double, 4UL> m_h;        // non-unity quaternions -- parametrizing the rotation
    blaze::StaticVector<double, 3UL> m_u, m_v;   // u: curvature vector | v: linear rate of change of position wrt arc-length
    blaze::StaticMatrix<double, 6UL, 6UL> m_Phi; // linear system matrix for determining vs, us
    blaze::StaticVector<double, 3UL> m_us;       // first spacial derivative of the curvature vector
    blaze::StaticVector<double, 3UL> m_vs;       // first spacial derivative of the strain vector
};