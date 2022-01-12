#include <blaze/Math.h>

template<std::size_t N>
class timeManager {
    public:

    // default constructor
	timeManager();

	// overloaded constructor
	timeManager(const double alpha, const double dt);

	// copy constructor
	timeManager(const timeManager& rhs);

	// move constructor
	timeManager(timeManager&& rhs) noexcept;

	// CTR destructor
	~timeManager();

	// copy assignment operator
	timeManager& operator=(const timeManager& rhs);

	// move assignment operator
	timeManager& operator=(timeManager&& rhs) noexcept;

    // getter method for returning the implicit time differentiation C0
    double getC0();

    // method for stepping time (updates the backward time differentiation)
    void stepTime(const blaze::StaticVector<double, 12UL> &currentState, const double& dt);

    // getter method for retrieving the dynamic actuation input profile for the tendons
    double getDynamicActuationProfile();

    private:
    // simulation run time
    double m_t;
    // parameters for implicit time differentiation (BDF-alpha)
    double m_c0;
    double m_c1;
    double m_c2;
    double m_d1;

    blaze::StaticMatrix<double, 12UL, N> m_yLag1;        // state vector one time step behind y(t-dt)   -- lag 1
    blaze::StaticMatrix<double, 12UL, N> m_yLag2;        // state vector two time steos behind y(t-2dt) -- lag 2
    blaze::StaticMatrix<double, 12UL, N> m_yTimeDot;     // time derivative of the state vector 
};