// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "TendonDriven.hpp"
#include <chrono>
#include <future> // needed in order to invoke async "execute functions asynchronously"


int main()
{
	//  # # # # # # # # ---- Properties of Nitinol Tubes ---- # # # # # # # #
	double E = 207.0e9; 				// Young's modulus [GPa]
	double G = 79.62e9; 				// Shear modulus [GPa]

	double radius = 0.00135 * 0.5;		// radius of the sheath [m]
	double mass = 0.034;				// total mass of the sheath [kg]
	double Length = 0.7144;				// total length of the sheath [m]
	// blaze::staticVector<double, 3UL> g = {0.0, 0.0, -9.81};	// acceleration of gravity [m/s/s]
	double T = 25;						// time horizon of the dynamic simulation [s]
	double dt = 5e-3;					// time infinitesimal [s]
	double basePlateToMotor = 0.0518;   // distance between base plate and motor [m]
	double tendonOffset = 0.01506;		// offset distance between the tendons and backbone [m]
	double beta = 0.0;					// length of tendons retracted within the actuation unit due to revolute actuator [m]
	double alpha = 0.0;					// backward differentiation parameter BDF-alpha

	double tendonCompliance = 1.6e-3;	// tendon compliance to tensile forces -- MADE EQUAL FOR ALL TENDONS
	

	// instantiating a tendon-driven robot (4 tendons, 200 discretized arc-length points)
	tendonDriven rbt<200UL, 4UL>();



	return 0;
}