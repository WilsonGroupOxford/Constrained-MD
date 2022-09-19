// Filthy awful hack to test.
#include "gtest/gtest.h"
#include <Eigen/Dense>
#define private public
#include "bond.h"


TEST(EnergyTest, ZeroEnergy)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 1.0, 0.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 0.0);
}

TEST(EnergyTest, QuadraticEnergy)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 2.0, 0.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 1.0 / 2.0);
    
    positions << 0.0, 0.0,
                 3.0, 0.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 4.0 / 2.0);
    
    positions << 0.0, 0.0,
                 4.0, 0.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 9.0 / 2.0);
}

TEST(EnergyTest, QuadraticEnergyOpposite)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 0.0, 0.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 1.0 / 2.0);
    
    positions << 0.0, 0.0,
                 -1.0, 0.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 4.0 / 2.0);
    
    positions << 0.0, 0.0,
                 -2.0, 0.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 9.0 / 2.0);
}

TEST(EnergyTest, ZeroEnergyOffRail)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 1.0, 1.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 0.0);
}

TEST(EnergyTest, QuadraticEnergyOffRail)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 2.0, 1.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 1.0 / 2.0);
    
    positions << 0.0, 0.0,
                 3.0, 1.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 4.0 / 2.0);
    
    positions << 0.0, 0.0,
                 4.0, 1.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 9.0 / 2.0);
}

TEST(EnergyTest, QuadraticEnergyOppositeOffRail)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 0.0, 1.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 1.0 / 2.0);
    
    positions << 0.0, 0.0,
                 -1.0, 1.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 4.0 / 2.0);
    
    positions << 0.0, 0.0,
                 -2.0, 1.0;
                 
    ASSERT_FLOAT_EQ(bond.energy(positions), 9.0 / 2.0);
}

TEST(PeriodTest, SimpleFreq)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    Eigen::VectorXd masses(2);
    masses << 1.0, 1.0;
    
    ASSERT_FLOAT_EQ(bond.frequency(masses), std::sqrt(2.0));
}

TEST(PeriodTest, DoubleFreq)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {2.0, 1.0, 1, 0, rail};
    Eigen::VectorXd masses(2);
    masses << 1.0, 1.0;
    
    ASSERT_FLOAT_EQ(bond.frequency(masses), 2.0);
}

TEST(PeriodTest, SimplePeriod)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    Eigen::VectorXd masses(2);
    masses << 1.0, 1.0;
    
    ASSERT_FLOAT_EQ(bond.period(masses), std::sqrt(2.0) * M_PI);
}

TEST(PeriodTest, DoublePeriod)
{
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {2.0, 1.0, 1, 0, rail};
    Eigen::VectorXd masses(2);
    masses << 1.0, 1.0;
    
    ASSERT_FLOAT_EQ(bond.period(masses), M_PI);
}
