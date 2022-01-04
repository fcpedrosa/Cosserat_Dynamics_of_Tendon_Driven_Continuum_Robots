#include "timeManager.hpp"


 // CHECK THIS ON HOW TO IMPLEMENT TEMPLATED CONSTRUCTORS: https://stackoverflow.com/questions/8305952/overloading-assignment-operator-in-a-class-template-that-can-cast-to-another-tem

// default constructor
template<std::size_t N>
timeManager<N>::timeManager()
{
    m_c0 = m_c1 = m_c2 = m_d1 = 0.0;
}

// overloaded constructor
template<std::size_t N>
timeManager<N>::timeManager(const double alpha, const double dt)
{
    this->m_c0 = (1.5 + alpha) / (dt * (1 + alpha));
    this->m_c1 = -2 / dt;
    this->m_c2 = (0.5 + alpha) / (dt * (1 + alpha));
    this->m_d1 = alpha / (1 + alpha);
}

// copy constructor
template<std::size_t N>
timeManager<N>::timeManager(const timeManager<N> &rhs) : m_c0(rhs.m_c0), m_c1(rhs.m_c1), m_c2(rhs.m_c2), m_d1(rhs.m_d1),
                                                   m_yLag1(rhs.m_yLag1), m_yLag2(rhs.m_yLag2), m_yTimeDot(rhs.m_yTimeDot) {}

// move constructor
template<std::size_t N>
timeManager<N>::timeManager(timeManager<N> &&rhs) noexcept
{
    // handling self-assignment
    if (this != &rhs)
    {
        this->m_c0 = rhs.m_c0;
        this->m_c1 = rhs.m_c1;
        this->m_c2 = rhs.m_c2;
        this->m_d1 = rhs.m_d1;
        this->m_yLag1 = std::move(rhs.m_yLag1);
        this->m_yLag2 = std::move(rhs.m_yLag2);
        this->m_yTimeDot = std::move(rhs.m_yTimeDot);
    }
}

// CTR destructor
template<std::size_t N>
timeManager<N>::~timeManager()
{
    // nothing to be done
}

// copy assignment operator
template<std::size_t N>
timeManager &timeManager<N>::operator=(const timeManager<N> &rhs)
{
    // handling self-assignment
    if (this != &rhs)
    {
        this->m_c0 = rhs.m_c0;
        this->m_c1 = rhs.m_c1;
        this->m_c2 = rhs.m_c2;
        this->m_d1 = rhs.m_d1;
        this->m_yLag1 = rhs.m_yLag1;
        this->m_yLag2 = rhs.m_yLag2;
        this->m_yTimeDot = rhs.m_yTimeDot;
    }

    return *this;
}

// move assignment operator
template<std::size_t N>
timeManager &timeManager<N>::operator=(timeManager<N> &&rhs) noexcept
{
    // handling self-assignment
    if (this != &rhs)
    {
        this->m_c0 = rhs.m_c0;
        this->m_c1 = rhs.m_c1;
        this->m_c2 = rhs.m_c2;
        this->m_d1 = rhs.m_d1;
        this->m_yLag1 = std::move(rhs.m_yLag1);
        this->m_yLag2 = std::move(rhs.m_yLag2);
        this->m_yTimeDot = std::move(rhs.m_yTimeDot);
    }

    return *this;
}

// getter method for returning the implicit time differentiation C0
template<std::size_t N>
double timeManager<N>::getC0()
{
    return this->m_c0;
}

// method for stepping time (updates the backward time differentiation)
template<std::size_t N>
void timeManager<N>::stepTime(const blaze::StaticVector<double, 12UL> &currentState)
{
    // updates the time derivative of the state vector
    this->m_yTimeDot = (this->m_c0 * currentState) + (this->m_c1 * this->m_yLag1) + (this->m_c2 * this->m_yLag2) + (this->m_d1 * this->m_yTimeDot);
    // updates the lagged state vectors
    this->m_yLag2 = std::move(this->m_yLag1);
    this->m_yLag1 = currentState;
}

// getter method for retrieving the dynamic actuation input profile for the tendons
template<std::size_t N>
double timeManager<N>::getDynamicActuationProfile()
{
    /*
                    *** THE SIMULATION STARTS WITH THE LOWEST TENDON DISPLACEMENT (qL) ***
        -- After t1 seconds: the displacement increases linearly toward its highest (qH) within (t2 - t1) seconds
        -- The tendon displacement remain at its highest value for (t3 - t2) seconds
        -- At t3 seconds: the displacement decreases linearly toward its lowest (qL) within (t4 - t3) seconds
        -- After t4 seconds: the displacement is kept unaltered at its lowest (qL) value
    */

    const double qL = 0.0;         // lowest tendon displacement value [meters]
    const double qH = -0.01619;    // lowest tendon displacement value [meters]
    const double t1 = 3.652;        
    const double t2 = 3.967;
    const double t3 = 17.63;
    const double t4 = 3.94;
    double q_t;                     // tendon displacement at time t [meters] 

    // computing the tendon displacment at each time
    if ( (this->m_t <= t1) || (this->m_t > t4) )
        q_t = qL;                                                   // region of lowest displacement
    else
    {
        if ( (this->m_t > t1) && (this->m_t <= t2) )
            q_t = qL + (qH - qL)*(this->m_t - t1)/(t2 - t1);        // ramp linearly toward the highest displacement
        else{
            if ( (this->m_t > t3) && (this->m_t <= t4) )
                q_t = qH + (qL - qH)*(this->m_t - t3)/(t4 - t3);    // ramp linearly toward lowest displacement
            else
                q_t = qH;                                           // region of highest displacement
        }
    }

    return q_t;
}