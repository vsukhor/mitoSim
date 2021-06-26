#include "gtest/gtest.h"

#include "../edge.h"

namespace edge_test {

TEST(EdgeTest, Constructor)
{
  mitosim::Edge<3> e {3, 4, 5};

  EXPECT_EQ(3, e.get_ind());
  EXPECT_EQ(4, e.get_indcl());
  EXPECT_EQ(5, e.get_cl());
  EXPECT_EQ(e.get_fin(0), 0.);
  EXPECT_EQ(e.get_fin(1), 0.);
}

TEST(EdgeTest, Reflect)
{
  mitosim::Edge<3> e {3, 4, 5};
  e.set_fin(0, 10.);
  e.reflect();

  EXPECT_EQ(e.get_fin(0), 0.);
  EXPECT_EQ(e.get_fin(1), 10.);
}

}  // namespace edge_test
