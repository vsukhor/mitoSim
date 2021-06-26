#include "gtest/gtest.h"

#include "../definitions.h"
#include "../segment.h"
#include "../structure.h"

namespace structure_test {

class StructureTest
    : public testing::Test {

protected:

    using Mt = mitosim::Segment<3>;
    using Structure = mitosim::Structure<Mt>;
    using Msgr = mitosim::Msgr;
    using szt = mitosim::szt;

    StructureTest()
        : msgr {}
        , stc {msgr}
    {}

    void SetUp() override {
        stc.add_disconnected_segment(4);
        stc.add_disconnected_segment(3);
    }

    Msgr msgr;
    Structure stc;
};

TEST_F(StructureTest, Constructor)
{
    Structure s {msgr};

    ASSERT_TRUE(s.clagl.empty());
    ASSERT_TRUE(s.glm.empty());
    ASSERT_TRUE(s.gla.empty());
    ASSERT_TRUE(s.mt.empty());
    for (const auto& n : s.nn)
        EXPECT_EQ(n, 0);
    EXPECT_EQ(s.mtnum, 0);
    EXPECT_EQ(s.clnum, 0);
    EXPECT_EQ(s.mtmass, 0);
    ASSERT_TRUE(s.clmt.empty());
    ASSERT_TRUE(s.mt11.empty());
    ASSERT_TRUE(s.mtc11.empty());
    ASSERT_TRUE(s.mt22.empty());
    ASSERT_TRUE(s.mtc22.empty());
    ASSERT_TRUE(s.mt33.empty());
    ASSERT_TRUE(s.mtc33.empty());
    ASSERT_TRUE(s.mt13.empty());
    ASSERT_TRUE(s.mtc13.empty());
}

TEST_F(StructureTest, AddDisconnectedSegment)
{
    szt len {4};
    Structure s {msgr};
    s.add_disconnected_segment(len);

    ASSERT_TRUE(s.clagl.empty());
    ASSERT_TRUE(s.glm.empty());
    ASSERT_TRUE(s.gla.empty());
    for (const auto& n : s.nn)
        EXPECT_EQ(n, 0);
    EXPECT_EQ(s.mtnum, 1);
    EXPECT_EQ(s.clnum, 1);
    EXPECT_EQ(s.mtmass, len);
    EXPECT_EQ(s.mt.size(), s.mtnum + 1);
    ASSERT_TRUE(s.clmt.empty());
    ASSERT_TRUE(s.mt11.empty());
    ASSERT_TRUE(s.mtc11.empty());
    ASSERT_TRUE(s.mt22.empty());
    ASSERT_TRUE(s.mtc22.empty());
    ASSERT_TRUE(s.mt33.empty());
    ASSERT_TRUE(s.mtc33.empty());
    ASSERT_TRUE(s.mt13.empty());
    ASSERT_TRUE(s.mtc13.empty());
}

TEST_F(StructureTest, UpdateNn1)
{
    stc.update_nn<1>();
    EXPECT_EQ(stc.nn[0], 2 * stc.mtnum);
    EXPECT_EQ(stc.nn[1], 0);
    EXPECT_EQ(stc.nn[2], 0);
}

TEST_F(StructureTest, UpdateNn2)
{
    stc.update_nn<2>();
    EXPECT_EQ(stc.nn[0], 0);
    EXPECT_EQ(stc.nn[1], stc.mt[1].g.size() + stc.mt[2].g.size() - 2);
    EXPECT_EQ(stc.nn[2], 0);
}

TEST_F(StructureTest, UpdateNn3)
{
    stc.update_nn<3>();
    EXPECT_EQ(stc.nn[0], 0);
    EXPECT_EQ(stc.nn[1], 0);
    EXPECT_EQ(stc.nn[2], 0);
}

TEST_F(StructureTest, UpdateNodeNumbers)
{
    stc.update_node_numbers();
    EXPECT_EQ(stc.nn[0], 2 * stc.mtnum);
    EXPECT_EQ(stc.nn[1], stc.mt[1].g.size() + stc.mt[2].g.size() - 2);
    EXPECT_EQ(stc.nn[2], 0);
}

}  // namespace structure_test
