// *************************************************************************************************************** //
// *        This file is part of a C++ kinematics & Dynamics library for Tendon-Driven Continuum Robots          * //
// *																					                         * //
// *    ------------ # Copyright (C) 2021 Dr Mohammad Salehizadeh <msalehizadeh@bwh.harvard.edu> # -----------   * //
// *    ------------ #           Copyright (C) 2021 Filipe C. Pedrosa <fpedrosa@uwo.ca>          # -----------   * //
// *																					                         * //
// *       Project developed under the supervision of Dr Jayender Jagadeesan (❖) Dr Rajni Patel (◈)	            * //
// *  =======================================================================================================    * //
// *  (◈) CSTAR (Canadian Surgical Technologies & Advanced Robotics) @ Western University, London, ON  Canada   * //
// *            (❖) Surgical Planning Lab @ Brigham and Women's Hospital / Harvard Medical School				* //
// *************************************************************************************************************** //

#pragma once

#include <blaze/Math.h>
#include <array>
#include <math.h>
#include <memory>
#include <execution>
#include <iostream>

// enumeration for defining which set of equations to be used (static or dynamic solutions)
enum class cosseratEquations { STATIC_SOLUTION, DYNAMIC_SOLUTION };

// enumeration for defining which root finding method to be used
enum class rootFindingMethod { NEWTON_RAPHSON, LEVENBERG_MARQUARDT,	POWELL_DOG_LEG,	MODIFIED_NEWTON_RAPHSON, BROYDEN };

// class that implements the dynamics of tendon-driven robots (any arbitrary number of tendons)
template <std::size_t N, std::size_t numTendons>
class TendonDriven
{
public:
    // removes the default constructor
    TendonDriven() = delete;

    // overloaded constructor
    TendonDriven(double E, double rad, double mass, double rho, double L, double beta, double T, double dt, double alpha, double tendonOffset, blaze::StaticVector<double, 3UL> g = {0.0, 0.0, -9.81});

    // copy constructor
    TendonDriven(const TendonDriven &rhs);

    // move constructor
    TendonDriven(TendonDriven &&rhs) noexcept;

    // CTR destructor
    ~TendonDriven();

    // copy assignment operator
    TendonDriven &operator=(const TendonDriven &rhs);

    // move assignment operator
    TendonDriven &operator=(TendonDriven &&rhs) noexcept;

    // returns the residue (violation of boundary condition) for a particular initial guess
    blaze::StaticVector<double, 6UL + numTendons> residueFunction(const blaze::StaticVector<double, 6UL + numTendons> &initGuess);

    // function that computes the finite-differences Jacobian for solving the BVP
	blaze::StaticMatrix<double, 6UL + numTendons, 6UL + numTendons> jac_BVP(const blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 6UL + numTendons> &residue);

	// function that computes the finite-differences Jacobian wrt actuation inputs
	blaze::StaticMatrix<double, 3UL, 6UL> jacobian(const blaze::StaticVector<double, 6UL + numTendons> &initGuess, const blaze::StaticVector<double, 3UL> &tipPos);

    // manages which method to be used for solving the Boundary-value problem
    bool solveBVP(blaze::StaticVector<double, 6UL + numTendons> &initGuess);

    // function that implements Powell's Dog Leg Method (Nonlinear root-finding method for solving the BVP)
    bool PowellDogLeg(blaze::StaticVector<double, 6UL + numTendons> &initGuess);

    // function that implements the Levenberg-Marquardt Method (Nonlinear root-finding method for solving the BVP)
    bool Levenberg_Marquardt(blaze::StaticVector<double, 6UL + numTendons> &initGuess);

    // function that implements Broyden's Nonlinear root-finding method for solving the BVP (Jacobian inverse is estimated)
    bool Broyden(blaze::StaticVector<double, 6UL + numTendons> &initGuess);

    // function that implements the Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
    bool Newton_Raphson(blaze::StaticVector<double, 6UL + numTendons> &initGuess);

    // function that implements the Modified, globally convergent Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
    bool Modified_Newton_Raphson(blaze::StaticVector<double, 6UL + numTendons> &initGuess);

    // method for computing the Penrose Pseudoinverse via SVD decomposition
    blaze::HybridMatrix<double, 6UL, 6UL> TendonDriven<N, numTendons>::pInv(const blaze::HybridMatrix<double, 6UL + numTendons, 6UL + numTendons> &M)

    // defining friend classes
    friend class ODESystem;
    friend class PDESystem;

private:
    double m_E;                                                         // Young's modulus [GPa]
    double m_rad;                                                       // backbone radius [m]
    double m_A;                                                         // cross-sectional area [m^2]
    double m_mass;                                                      // total mass (backbone + disks) [kg]
    double m_rho;                                                       // material density [kgm^-3]
    blaze::StaticVector<double, 3UL> m_g;                               // acceleration of gravity [m/s^2]
    double m_L;                                                         // backbone total length [m]
    double m_beta;                                                      // length of tendons retracted within the actuation unit due to revolute actuator [m]
    double m_T, m_dt;                                                   // total time horizon & time step for dynamics simulation [s]
    double m_alpha;                                                     // backward differentiation parameter (BDF_alpha)
    double m_tendonOffset;                                              // offset distance between the tendons and backbone [m]
    blaze::StaticVector<double, numTendons> m_tendonSlack;              // eventual slack buildup on tendons
    double m_basePlateToMotor;                                          // distance between base plate and motor (tendon length within the actuation unit)
    blaze::StaticVector<double, numTendons> m_tau;                      // tensile actuation force on tendons
    blaze::StaticVector<double, numTendons> m_tendonCompliance;         // tendon compliance to tensile forces
    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL>> m_Bse; // shear & extension viscous-like damping [Ns]
    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL>> m_Bbt; // bending & torsion viscous-like damping [Nm^2s]
    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL>> m_C;   // air drag viscous-like damping [kg/m^2]
    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL>> m_Kse; // stiffness matrix for sear & extension [N]
    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL>> m_Kbt; // stiffness matrix for bending & torsion [Nm^2]
    blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL>> m_J;   // second mass moment of inertia tensor [m^4]
    std::array<blaze::StaticVector<double, 3UL>, numTendons> m_r;       // array of vectors describing the radial position of each tendon
    std::array<blaze::StaticVector<double, 3UL>, numTendons> m_p;       // array of vectors describing the shape of each tendon
    std::unique_ptr<timeManager> m_time;                                // implements the time manager for dynamics simulation
    std::unique_ptr<ODESystem<N, numTendons>> m_ODE;                    // implements the ODE state equations for a static tendon-driven robot (steady-state)
    std::unique_ptr<PDESystem<N, numTendons>> m_PDE;                    // implements the PDE state equations for a dynamic tendon-driven robot
    blaze::StaticMatrix<double, 19UL + numTendons, N> m_y;              // set of N state vectors
    blaze::StaticMatrix<double, 18UL, N> m_z, m_zh;                     // set of N
    cosseratEquations m_solType;                                        // enum defining if ODE or PDE equations for solving BVP
    rootFindingMethod m_method;                                         // enum defining which nonlinear root-finder to be used
    double m_BVPAccuracy;                                               // accuracy of BVP solution
};