// Filthy awful hack to test.
#include "gtest/gtest.h"
#include <Eigen/Dense>
#define private public
#include "bond.h"


TEST(ExcitementTest, ZeroExcitement)
{
    // Test the excitement factor of a single bond that isn't stretched
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 1.0, 0.0;

    ASSERT_FLOAT_EQ(bond.get_excitement_factor(positions), 1.0);
}

TEST(ExcitementTest, NegativeExcitement)
{
    // Test the excitement factor of a single bond that is stretched the other way
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 1.0, 0.0,
                 0.0, 0.0;

    ASSERT_FLOAT_EQ(bond.get_excitement_factor(positions), -1.0);
}

TEST(ExcitementTest, DoubleNegativeExcitement)
{
    // Test the excitement factor of a single bond that is stretched the other way
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 2.0, 0.0,
                 0.0, 0.0;

    ASSERT_FLOAT_EQ(bond.get_excitement_factor(positions), -2.0);
}

TEST(ExcitementTest, DoubleExcitement)
{
    // Test the excitement factor of a single bond that is stretched the other way
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 2.0, 0.0;

    ASSERT_FLOAT_EQ(bond.get_excitement_factor(positions), 2.0);
}

TEST(ExcitementTest, OffRailDoubleExcitement)
{
    // Test the excitement factor of a single bond that is off the rail
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 1, 0, rail};
    
    Eigen::MatrixXd positions(2, 2);
    positions << 0.0, 0.0,
                 2.0, 1.0;

    ASSERT_FLOAT_EQ(bond.get_excitement_factor(positions), 2.0);
}
