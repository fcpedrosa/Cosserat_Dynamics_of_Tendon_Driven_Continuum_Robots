#include "PDESystem.hpp"
#include "tendonDriven.hpp"
#include <math.h>

// default constructor
template <std::size_t N, std::size_t numTendons>
PDESystem::PDESystem()
{
}

// overloaded constructor
template <std::size_t N, std::size_t numTendons>
PDESystem::PDESystem(blaze::StaticVector<double, 3UL> u, blaze::StaticVector<double, 3UL> v, blaze::StaticVector<double, 4UL> h)
{
}

// copy constructor
template <std::size_t N, std::size_t numTendons>
PDESystem::PDESystem(const PDESystem &rhs)
{
}

// move constructor
template <std::size_t N, std::size_t numTendons>
PDESystem::PDESystem(PDESystem &&rhs) noexcept
{
}

// CTR destructor
template <std::size_t N, std::size_t numTendons>
PDESystem::~PDESystem()
{
    // nothing to be done
}

// copy assignment operator
template <std::size_t N, std::size_t numTendons>
PDESystem &PDESystem::operator=(const PDESystem &rhs)
{
}

// move assignment operator
template <std::size_t N, std::size_t numTendons>
PDESystem &PDESystem::operator=(PDESystem &&rhs) noexcept
{
}

// functor that implements the system of PDEs governing an N-tendon catheter
template <std::size_t N, std::size_t numTendons>
void PDESystem::operator()(const blaze::StaticVector<double, 19UL + numTendons> &y, blaze::StaticVector<double, 19UL + numTendons> &dyds, const blaze::StaticVector<double, 18UL> &z_h, blaze::StaticVector<double, 18UL> &z, const tendonDriven &robot)
{
    // dyds := [ps, hs, vs, us, qs, ws, numTendons] ==>> First spatial derivative of the (19+numTendons)-dimensional state vector
    // zh := [vh uh qh wh]  ==>> history terms for the time implicit discretization
    blaze::StaticVector<double, 3UL> v, u, q, w;
    v = blaze::subvector<7UL, 3UL>(y);
    u = blaze::subvector<10UL, 3UL>(y);
    q = blaze::subvector<13UL, 3UL>(y);
    w = blaze::subvector<16UL, 3UL>(y);
    blaze::StaticVector<double, 4UL> h = blaze::subvector<3UL, 4UL>(y);
    blaze::StaticMatrix<double, 3UL, 3UL> R = this->getSO3(h);

    // obtaining the BDG-alpha time discretization coefficient C0
    double c0 = robot.m_time->getC0();

    // declaring auxiliar variables
    blaze::StaticVector<double, 3UL> a, a_i, b, pbs_i;
    blaze::StaticMatrix<double, 3UL, 3UL> A_i, G, G_i;
    blaze::StaticMatrix<double, 6UL, 6UL> Phi;
    blaze::StaticVector<double, numTendons> pbsNorm = blaze::subvector<20UL, numTendons>(dyds);

    auto Blk1 = blaze::submatrix<0UL, 0UL, 3UL, 3UL>(Phi); // block matrix corresponding to: Kse + A
    auto Blk2 = blaze::submatrix<3UL, 3UL, 3UL, 3UL>(Phi); // block matrix corresponding to: Kbt + H
    Blk1 = robot.m_Kse + (c0 * robot.m_Bse);
    Blk2 = robot.m_Kbt + (c0 * robot.m_Bbt);

    // iterating through the tendons in the robot design
    for (size_t idx = 0UL; idx < numTendons; ++idx)
    {
        pbs_i = blaze::cross(u, robot.m_r[idx]) + v;
        pbsNorm[idx] = blaze::norm(pbs_i);
        A_i = -this->hatSqr(pbs_i) * (robot.m_tau[idx] / pow(pbsNorm[idx], 3));
        G_i = -this->hatPostMultiply(A_i, robot.m_r[idx]);
        a_i = A_i * blaze::cross(u, pbs_i);

        a += a_i;
        b += blaze::cross(robot.m_r[idx], a_i);
        Blk1 += A_i;
        G += G_i;
        Blk2 += this->hatPreMultiply(robot.m_r[idx], G_i);
    }

    blaze::submatrix<0UL, 3UL, 3UL, 3UL>(Phi) = G;
    blaze::submatrix<3UL, 0UL, 3UL, 3UL>(Phi) = blaze::trans(G);

    blaze::StaticVector<double, 3UL> vs_h, us_h;
    vs_h = blaze::subvector<12UL, 3UL>(z_h);
    vs_h = blaze::subvector<15UL, 3UL>(z_h);

    // computing time derivatives of the state vector from time discretization
    blaze::StaticVector<double, 3UL> v_t, u_t, q_t, w_t;
    v_t = (c0 * v) + blaze::subvector<0UL, 3UL>(z_h);
    u_t = (c0 * u) + blaze::subvector<3UL, 3UL>(z_h);
    q_t = (c0 * q) + blaze::subvector<6UL, 3UL>(z_h);
    w_t = (c0 * w) + blaze::subvector<9UL, 3UL>(z_h);

    blaze::StaticVector<double, 3UL> nb, mb, e_3;
    e_3 = {0, 0, 1};
    nb = robot.m_Kse * (v - e_3) + robot.m_Bse * v_t;
    mb = robot.m_Kbt * u + robot.m_Bbt * u_t;

    // right hand side vector -- for determining us, vs
    blaze::StaticVector<double, 6UL> rhs;
    blaze::subvector<0UL, 3UL>(rhs) = -a + robot.m_rho * robot.m_A * (blaze::cross(w, q) + q_t) + robot.m_C * (q * blaze::abs(q)) - blaze::trans(R) * robot.m_rho * robot.m_A * robot.m_g - blaze::cross(u, nb) - robot.m_Bse * vs_h;
    blaze::subvector<3UL, 3UL>(rhs) = -b + blaze::cross(w, robot.m_rho * robot.m_J * w) + robot.m_rho * robot.m_J * w_t - blaze::cross(v, nb) - blaze::cross(u, mb) - robot.m_Bbt * us_h;

    // State differential equations
    blaze::subvector<0UL, 3UL>(dyds) = R * v;                                          // p_s
    blaze::subvector<3UL, 4UL>(dyds) = this->quaternionDiff(u, h);                     // h_s
    blaze::subvector<7UL, 6UL>(dyds) = blaze::solve(Phi, rhs);                         // v_s & u_s
    blaze::subvector<13UL, 3UL>(dyds) = v_t - blaze::cross(u, q) + blaze::cross(w, v); // q_s
    blaze::subvector<16UL, 3UL>(dyds) = u_t - blaze::cross(u, w);                      // w_s

    blaze::subvector<0UL, 3UL>(z) = v;                                 // v
    blaze::subvector<3UL, 3UL>(z) = u;                                 // u
    blaze::subvector<6UL, 3UL>(z) = q;                                 // q
    blaze::subvector<9UL, 3UL>(z) = w;                                 // w
    blaze::subvector<12UL, 6UL>(z) = blaze::subvector<7UL, 6UL>(dyds); // v_s & u_s
}

// Integrates the PDE model equations using Euler's method
template <std::size_t N, std::size_t numTendons>
void PDESystem::euler_PDE(const blaze::StaticVector<double, 19UL + numTendons> &y0, double, 19UL + numTendons > &y, blaze::StaticVector<double, 19UL + numTendons> &dyds, const blaze::StaticVector<double, 18UL> &z_h, blaze::StaticVector<double, 18UL> &z, const tendonDriven &robot)
{
    blaze::col<0UL>(Y) = y0;
    constexpr double ds = L / (N - 1);

    // Euler's method
    for (size_t i = 0UL; i < N - 1; ++i)
    {
        this->PDESystem(blaze::column(y, i), dyds, blaze::column(Z_h, i), blaze::column(Z, i), robot);
        blaze::column(Y, i + 1) = blaze::column(Y, i) + ds * dyds;
    }
}

// Integrates the ODE model equations using the classic 4th-order Runge-Kutta algorithm
template <std::size_t N, std::size_t numTendons>
void PDESystem::rungeKutta_PDE(const blaze::StaticVector<double, 19UL + numTendons> &y0, double, 19UL + numTendons > &y, blaze::StaticVector<double, 19UL + numTendons> &dyds, const blaze::StaticVector<double, 18UL> &z_h, blaze::StaticVector<double, 18UL> &z, const tendonDriven &robot)
{
}

// function that returns a rotation matrix in SO(3) from a set of non-unity quaternions
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 3UL> PDESystem::getSO3(const blaze::StaticVector<double, 4UL> &h)
{
    blaze::IdentityMatrix<double> I(3UL);
    blaze::StaticMatrix<double, 3UL, 3UL> R;

    R = {{-h[2UL] * h[2UL] - h[3UL] * h[3UL], h[1UL] * h[2UL] - h[3UL] * h[0UL], h[1UL] * h[3UL] + h[2UL] * h[0UL]},
         {h[1UL] * h[2UL] + h[3UL] * h[0UL], -h[1UL] * h[1UL] - h[3UL] * h[3UL], h[2UL] * h[3UL] - h[1UL] * h[0UL]},
         {h[1UL] * h[3UL] - h[2UL] * h[0UL], h[2UL] * h[3UL] + h[1UL] * h[0UL], -h[1UL] * h[1UL] - h[2UL] * h[2UL]}};

    R *= 2 / blaze::sqrNorm(h);
    R += I;

    return R;
}

// function that computes the squared hat operator
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 3UL> PDESystem::hatSqr(const blaze::StaticVector<double, 3UL> &v)
{
    blaze::StaticMatrix<double, 3UL, 3UL> hatSqr = {{-v[2UL] * v[2UL] - v[1UL] * v[1UL], v[1UL] * v[0UL]},
                                                    {v[0UL] * v[1UL], -v[2UL] * v[2UL] - v[0UL] * v[0UL], v[2UL] * v[1UL]},
                                                    {v[0UL] * v[2UL], -v[1UL] * v[2UL], -v[1UL] * v[1UL] - v[0UL] * v[0UL]}};

    return hatSqr;
}

// function that computes the premultiplication of a vector by the hat operator
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 3UL> PDESystem::hatPreMultiply(const blaze::StaticVector<double, 3UL> &v, const blaze::StaticMatrix<double, 3UL 3UL> &M)
{
    blaze::StaticMatrix<double, 3UL, 3UL> Res = {{-M(1UL, 0UL) * v[2UL] + M(2UL, 0UL) * v[1UL], -M(1UL, 1UL) * v[2UL] + M(2UL, 1UL) * v[1UL], -M(1UL, 2UL) * v[2UL] + M(2UL, 2UL) * v[1UL]},
                                                 {M(0UL, 0UL) * v[2UL] - M(2UL, 0UL) * v[0UL], M(0UL, 1UL) * v[2UL] - M(2UL, 1UL) * v[0UL], M(0UL, 2UL) * v[2UL] - M(2UL, 2UL) * v[0UL]},
                                                 {-M(0UL, 0UL) * v[1UL] + M(1UL, 0UL) * v[0UL], -M(0UL, 1UL) * v[1UL] + M(1UL, 1UL) * v[0UL], -M(0UL, 2UL) * v[1UL] + M(1UL, 2UL) * v[0UL]}};

    return Res;
}

// function that computes the posmultiplication of a vector by the hat operator
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 3UL> PDESystem::hatPostMultiply(const blaze::StaticMatrix<double, 3UL 3UL> &M, const blaze::StaticVector<double, 3UL> &v)
{
    blaze::StaticMatrix<double, 3UL, 3UL> Res = {{M(0UL, 1UL) * v[2UL] - M(0UL, 2UL) * v[1UL], -M(0UL, 0UL) * v[2UL] + M(0UL, 2UL) * v[0UL], M(0UL, 0UL) * v[1UL] - M(0UL, 1UL) * v[0UL]},
                                                 {M(1UL, 1UL) * v[2UL] - M(1UL, 2UL) * v[1UL], -M(1UL, 0UL) * v[2UL] + M(1UL, 2UL) * v[0UL], M(1UL, 0UL) * v[1UL] - M(1UL, 1UL) * v[0UL]},
                                                 {M(2UL, 1UL) * v[2UL] - M(2UL, 2UL) * v[1UL], -M(2UL, 0UL) * v[2UL] + M(2UL, 2UL) * v[0UL], M(2UL, 0UL) * v[1UL] - M(2UL, 1UL) * v[0UL]}};

    return Res;
}

// function that computes the differential quaternion evolution
template <std::size_t N, std::size_t numTendons>
blaze::StaticVector<double, 4UL> PDESystem::quaternionDiff(const blaze::StaticVector<double, 3UL> &u, const blaze::StaticVector<double, 4UL> &h)
{
    blaze::StaticMatrix<double, 4UL, 4UL> hs;
    hs = {{-u[0UL] * h[1UL] - u[1UL] * h[2UL] - u[2UL] * h[3UL]},
          {u[0UL] * h[0UL] + u[2UL] * h[2UL] - u[1UL] * h[3UL]},
          {u[1UL] * h[0UL] - u[2UL] * h[1UL] + u[0UL] * h[3UL]},
          {u[2UL] * h[0UL] + u[1UL] * h[1UL] - u[0UL] * h[2UL]}};

    return 0.5 * hs;
}