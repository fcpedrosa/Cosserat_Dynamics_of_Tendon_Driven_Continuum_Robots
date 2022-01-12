#include "ODESystem.hpp"
#include "tendonDriven.hpp"
#include <math.h>

// default constructor
template <std::size_t N, std::size_t numTendons>
ODESystem<N, numTendons>::ODESystem()
{
    m_u = 0.0;
    m_v = 0.0;
    m_h = {1.0, 0.0, 0.0, 0.0};
    m_R = this->getSO3(this->m_h);
    m_Phi = 0.0; 
}

// overloaded constructor
template <std::size_t N, std::size_t numTendons>
ODESystem<N, numTendons>::ODESystem(blaze::StaticVector<double, 3UL> u, blaze::StaticVector<double, 3UL> v, blaze::StaticVector<double, 4UL> h) : m_u(u), m_v(u), m_h(h)
{
    m_R = this->getSO3(this->m_h);    
}

// copy constructor
template <std::size_t N, std::size_t numTendons>
ODESystem<N, numTendons>::ODESystem(const ODESystem &rhs) : m_u(rhs.m_u), m_v(rhs.m_v), m_h(rhs.m_h), m_R(rhs.m_R), m_Phi(rhs.m_Phi) {}

// move constructor
template <std::size_t N, std::size_t numTendons>
ODESystem<N, numTendons>::ODESystem(ODESystem<N, numTendons> &&rhs) noexcept
{
    // handling self-assignment
    if (this != &rhs)
    {
        this->m_u = std::move(rhs.m_u);
        this->m_v = std::move(rhs.m_v);
        this->m_h = std::move(rhs.m_h);
        this->m_R = std::move(rhs.m_R);
        this->m_Phi = std::move(rhs.m_Phi);
    }
}

// CTR destructor
template <std::size_t N, std::size_t numTendons>
ODESystem<N, numTendons>::~ODESystem()
{
    // nothing to be done
}

// copy assignment operator
template <std::size_t N, std::size_t numTendons>
ODESystem<N, numTendons> &ODESystem<N, numTendons>::operator=(const ODESystem<N, numTendons> &rhs)
{
    // handling self-assignment
    if (this != &rhs)
    {
        this->m_u = rhs.m_u;
        this->m_v = rhs.m_v;
        this->m_h = rhs.m_h;
        this->m_R = rhs.m_R;
        this->m_Phi = rhs.m_Phi;
    }

    return *this;
}

// move assignment operator
template <std::size_t N, std::size_t numTendons>
ODESystem<N, numTendons> &ODESystem<N, numTendons>::operator=(ODESystem<N, numTendons> &&rhs) noexcept
{
    // handling self-assignment
    if (this != &rhs)
    {
        this->m_u = std::move(rhs.m_u);
        this->m_v = std::move(rhs.m_v);
        this->m_h = std::move(rhs.m_h);
        this->m_R = std::move(rhs.m_R);
        this->m_Phi = std::move(rhs.m_Phi);
    }

    return *this;
}

// functor that implements the system of ODEs governing an N-tendon catheter
template <std::size_t N, std::size_t numTendons>
void ODESystem<N, numTendons>::operator()(const blaze::StaticVector<double, 19UL + numTendons> &y, blaze::StaticVector<double, 19UL + numTendons> &dyds, blaze::StaticVector<double, 18UL> &z, const tendonDriven &robot)
{
    // dyds := [ps, hs, vs, us, qs, ws, numTendons] ==>> First spatial derivative of the (19+numTendons)-dimensional state vector 
    blaze::StaticVector<double, 3UL> v = blaze::subvector<7UL, 3UL>(y);
    blaze::StaticVector<double, 3UL> u = blaze::subvector<10UL, 3UL>(y);
    blaze::StaticVector<double, 4UL> h = blaze::subvector<3UL, 4UL>(y);
    blaze::StaticMatrix<double, 3UL, 3UL> R = this->getSO3(h);

    // declaring auxiliar variables
    blaze::StaticVector<double, 3UL> a, a_i, b, pbs_i;
    blaze::StaticMatrix<double, 3UL, 3UL> A_i, G, G_i;
    blaze::StaticMatrix<double, 6UL, 6UL> Phi;
    blaze::StaticVector<double, numTendons> pbsNorm = blaze::subvector<20UL, numTendons>(dyds);

    auto Blk1 = blaze::submatrix<0UL, 0UL, 3UL, 3UL>(Phi);  // block matrix corresponding to: Kse + A
    auto Blk2 = blaze::submatrix<3UL, 3UL, 3UL, 3UL>(Phi);  // block matrix corresponding to: Kbt + H
    Blk1 = robot.m_Kse;
    Blk2 = robot.m_Kbt;

    // iterating through the tendons in the robot design
    for (size_t idx = 0UL; idx < numTendons; ++idx)
    {
        pbs_i = blaze::cross(u, robot.m_r[idx]) + v;
        pbsNorm[idx] = blaze::norm(pbs_i);
        A_i = -this->hatSqr(pbs_i) * ( robot.m_tau[idx] / pow(pbsNorm[idx], 3) ); 
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

    blaze::StaticVector<double, 3UL> nb, mb, e3;
    e3 = { 0, 0, 1 };
    nb = robot.m_Kse * (v - e3);
    mb = robot.m_Kbt * u;

    blaze::StaticVector<double, 6UL> rhs;
    blaze::subvector<0UL, 3UL>(rhs) = -blaze::cross(u, nb) - blaze::trans(R) * robot.m_rho * robot.m_A * robot.m_g - a;
    blaze::subvector<3UL, 3UL>(rhs) = -blaze::trans(u, mb) - blaze::cross(v, nb) - b;

    // ODE state equations
    blaze::subvector<0UL, 3UL>(dyds) = R * v;                               // p_s
    blaze::subvector<0UL, 4UL>(dyds) = this->quaternionDiff(u, h);          // h_s
    blaze::subvector<7UL, 6UL>(dyds) = blaze::solve(Phi, rhs);              // v_s & u_s
    blaze::subvector<13UL, 6UL>(dyds) = 0.0;                                // q_s & w_s (zero linear & angular velocities)

    // outputting elements for variables with time derivatives
    blaze::subvector<0UL, 3UL>(z) = v;                                   // v
    blaze::subvector<3UL, 3UL>(z) = u;                                   // u
    blaze::subvector<6UL, 6UL>(z) = 0.0;                                 // q_s = w_s = 0.0 (zero linear & angular velocities)
    blaze::subvector<12UL, 6UL>(z) = blaze::subvector<7UL, 6UL>(dyds);   // v_s & u_s
}

// Integrates the ODE model equations using Euler's method
template <std::size_t N, std::size_t numTendons>
void ODESystem<N, numTendons>::euler_ODE(const blaze::StaticVector<double, 19UL + numTendons>& y0, blaze::StaticMatrix<double, 19UL + numTendons, N>& Y, blaze::StaticVector<double, 19UL + numTendons>& dyds, blaze::StaticMatrix<double, 18UL + numTendons, N>& Z, const double L, tendonDriven &robot)
{
    blaze::col<0UL>(Y) = y0;
    constexpr double ds = L / (N - 1);

    // Euler's method
    for(size_t i = 0UL; i < N - 1; ++i)
    {
        this->ODESystem(blaze::column(y, i), dyds, blaze::column(Z, i), robot);
        blaze::column(Y, i + 1) = blaze::column(Y, i) + ds * dyds;
    }
}

// Integrates the ODE model equations using the classic 4th-order Runge-Kutta algorithm
template <std::size_t N, std::size_t numTendons>
void ODESystem<N, numTendons>::rungeKutta_ODE(const blaze::StaticVector<double, 19UL + numTendons>& y0, blaze::StaticMatrix<double, 19UL + numTendons, N>& Y, blaze::StaticMatrix<double, 18UL + numTendons, N>& Z, const double L, tendonDriven &robot)
{
    blaze::column<0UL>(Y) = y0;

    constexpr double ds = L / (N - 1);
    constexpr double half_ds = ds / 2;
    constexpr double sixth_ds = ds / 6;
    
    // 4th-order Runge-Kutta
    blaze::StaticVector<double, 19UL + numTendons> k0, y1, k1, y2, k2, y3, k3;
    double s = s0;
    for(size_t i = 0UL; i < N-1; ++i){
        y0 = Y.col(i);
        k0 = ODE(s, y0);
        s += half_ds;
        y1 = Y.col(i) + k0*half_ds;
        k1 = ODE(s, y1);
        y2 = Y.col(i) + k1*half_ds;
        k2 = ODE(s, y2);
        s = s0 + (L*(i+1))/Nm1;
        y3 = Y.col(i) + k2*ds;
        k3 = ODE(s, y3);

        Y.col(i+1) = Y.col(i) + (k0 + 2*(k1 + k2) + k3) * sixth_ds;
    }

    return Y;
}

// function that returns a rotation matrix in SO(3) from a set of non-unity quaternions
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 3UL> ODESystem<N, numTendons>::getSO3(const blaze::StaticVector<double, 4UL> &h)
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
blaze::StaticMatrix<double, 3UL, 3UL> ODESystem<N, numTendons>::hatSqr(const blaze::StaticVector<double, 3UL> &v)
{
    blaze::StaticMatrix<double, 3UL, 3UL> hatSqr = {{-v[2UL] * v[2UL] - v[1UL] * v[1UL], v[1UL] * v[0UL]},
                                                    {v[0UL] * v[1UL], -v[2UL] * v[2UL] - v[0UL] * v[0UL], v[2UL] * v[1UL]},
                                                    {v[0UL] * v[2UL], -v[1UL] * v[2UL], -v[1UL] * v[1UL] - v[0UL] * v[0UL]}};

    return hatSqr;
}

// function that computes the premultiplication of a vector by the hat operator
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 3UL> ODESystem<N, numTendons>::hatPreMultiply(const blaze::StaticVector<double, 3UL> &v, const blaze::StaticMatrix<double, 3UL 3UL> &M)
{
    blaze::StaticMatrix<double, 3UL, 3UL> Res = {{-M(1UL, 0UL) * v[2UL] + M(2UL, 0UL) * v[1UL], -M(1UL, 1UL) * v[2UL] + M(2UL, 1UL) * v[1UL], -M(1UL, 2UL) * v[2UL] + M(2UL, 2UL) * v[1UL]},
                                                 {M(0UL, 0UL) * v[2UL] - M(2UL, 0UL) * v[0UL], M(0UL, 1UL) * v[2UL] - M(2UL, 1UL) * v[0UL], M(0UL, 2UL) * v[2UL] - M(2UL, 2UL) * v[0UL]},
                                                 {-M(0UL, 0UL) * v[1UL] + M(1UL, 0UL) * v[0UL], -M(0UL, 1UL) * v[1UL] + M(1UL, 1UL) * v[0UL], -M(0UL, 2UL) * v[1UL] + M(1UL, 2UL) * v[0UL]}};

    return Res;
}

// function that computes the posmultiplication of a vector by the hat operator
template <std::size_t N, std::size_t numTendons>
blaze::StaticMatrix<double, 3UL, 3UL> ODESystem<N, numTendons>::hatPostMultiply(const blaze::StaticMatrix<double, 3UL 3UL> &M, const blaze::StaticVector<double, 3UL> &v)
{
    blaze::StaticMatrix<double, 3UL, 3UL> Res = {{M(0UL, 1UL) * v[2UL] - M(0UL, 2UL) * v[1UL], -M(0UL, 0UL) * v[2UL] + M(0UL, 2UL) * v[0UL], M(0UL, 0UL) * v[1UL] - M(0UL, 1UL) * v[0UL]},
                                                 {M(1UL, 1UL) * v[2UL] - M(1UL, 2UL) * v[1UL], -M(1UL, 0UL) * v[2UL] + M(1UL, 2UL) * v[0UL], M(1UL, 0UL) * v[1UL] - M(1UL, 1UL) * v[0UL]},
                                                 {M(2UL, 1UL) * v[2UL] - M(2UL, 2UL) * v[1UL], -M(2UL, 0UL) * v[2UL] + M(2UL, 2UL) * v[0UL], M(2UL, 0UL) * v[1UL] - M(2UL, 1UL) * v[0UL]}};

    return Res;
}

// function that computes the differential quaternion evolution
template <std::size_t N, std::size_t numTendons>
blaze::StaticVector<double, 4UL> PDESystem<N, numTendons>::quaternionDiff(const blaze::StaticVector<double, 3UL> &u, const blaze::StaticVector<double, 4UL> &h)
{
    blaze::StaticMatrix<double, 4UL, 4UL> hs;
    hs = {{-u[0UL] * h[1UL] - u[1UL] * h[2UL] - u[2UL] * h[3UL]},
          {u[0UL] * h[0UL] + u[2UL] * h[2UL] - u[1UL] * h[3UL]},
          {u[1UL] * h[0UL] - u[2UL] * h[1UL] + u[0UL] * h[3UL]},
          {u[2UL] * h[0UL] + u[1UL] * h[1UL] - u[0UL] * h[2UL]}};

    return 0.5 * hs;
}