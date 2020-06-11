// Filthy awful hack to test.
#include "gtest/gtest.h"
#include <Eigen/Dense>
#define private public
#include "bond.h"


TEST(ProjectionTest, OnRailUnchanged)
{
    // Test if projecting a vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << 1.0, 0.0;

    auto projected_vec = bond.project_onto_rail(new_vec);
    
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), new_vec(i));
    }
}

TEST(ProjectionTest, UnnormalisedRail)
{
    // Test if projecting a vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << 10.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << 1.0, 0.0;

    auto projected_vec = bond.project_onto_rail(new_vec);
    
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), new_vec(i));
    }
}

TEST(ProjectionTest, OnRailOpposite)
{
    // Test if projecting negative vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << -1.0, 0.0;

    auto projected_vec = bond.project_onto_rail(new_vec);
    
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), new_vec(i));
    }
}

TEST(ProjectionTest, OnRailLarger)
{
    // Test if projecting a larger vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << 2.0, 0.0;

    auto projected_vec = bond.project_onto_rail(new_vec);
    
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), new_vec(i));
    }
}

TEST(ProjectionTest, OnRailOppositeLarger)
{
    // Test if projecting a larger negative vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << -2.0, 0.0;

    auto projected_vec = bond.project_onto_rail(new_vec);
    
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), new_vec(i));
    }
}

TEST(ProjectionTest, OffRail)
{
    // Test if projecting a vector that's off a rail back onto the rail
    // changes it.
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << 1.0, 1.0;

    Eigen::VectorXd expected_vec(2);
    expected_vec << 1.0, 0.0;
    auto projected_vec = bond.project_onto_rail(new_vec);

    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), expected_vec(i));
    }
}

TEST(ProjectionTest, OffRailNegative)
{
    // Test if projecting a vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << 1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << -1.0, -1.0;

    Eigen::VectorXd expected_vec(2);
    expected_vec << -1.0, 0.0;
    
    auto projected_vec = bond.project_onto_rail(new_vec);
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), expected_vec(i));
    }
}

TEST(ProjectionTest, NegativeRail)
{
    // Test if projecting negative vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << -1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << -1.0, 0.0;

    auto projected_vec = bond.project_onto_rail(new_vec);
    
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), new_vec(i));
    }
}

TEST(ProjectionTest, NegativeRailOppositeVec)
{
    // Test if projecting negative vector onto itself leaves the vector unchanged.
    Eigen::VectorXd rail(2);
    rail << -1.0, 0.0;
    const HarmonicBond bond {1.0, 1.0, 0, 1, rail};
    
    Eigen::VectorXd new_vec(2);
    new_vec << 1.0, 0.0;

    auto projected_vec = bond.project_onto_rail(new_vec);
    
    for (int i = 0; i < new_vec.rows(); ++i) {
        ASSERT_FLOAT_EQ(projected_vec(i), new_vec(i));
    }
}
