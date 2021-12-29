#pragma once

#include <blaze/Math.h>

typedef blaze::StaticVector<double, 24UL> state_type;

class PDESystem {
    public:


    private:
    blaze::StaticVector<double, 3UL> m_us;                  // first spacial derivative of the curvature vector 
    blaze::StaticVector<double, 3UL> m_vs;                  // first spacial derivative of the strain vector
};