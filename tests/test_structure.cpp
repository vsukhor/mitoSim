#include "gtest/gtest.h"

#include "../definitions.h"
#include "../segment.h"
#include "../structure.h"

namespace structure_test {

class StructureTest
    : public testing::Test {

protected:

    using Msgr = mitosim::Msgr;
    using Mt = mitosim::Segment<3>;
    using Structure = mitosim::Structure<Mt>;
    using szt = mitosim::szt;

    StructureTest()
        : msgr {}
    {}


    Msgr msgr;
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
    constexpr szt len {4};
    Structure s {msgr};
    s.add_disconnected_segment(len);

    ASSERT_EQ(s.clagl.size(), 1);
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
    Structure s {msgr};
    s.add_disconnected_segment(4);
    s.add_disconnected_segment(3);

    s.update_nn<1>();
    EXPECT_EQ(s.nn[0], 2 * s.mtnum);
    EXPECT_EQ(s.nn[1], 0);
    EXPECT_EQ(s.nn[2], 0);
}

TEST_F(StructureTest, UpdateNn2)
{
    Structure s {msgr};
    s.add_disconnected_segment(4);
    s.add_disconnected_segment(3);

    s.update_nn<2>();
    EXPECT_EQ(s.nn[0], 0);
    EXPECT_EQ(s.nn[1], s.mt[1].g.size() + s.mt[2].g.size() - 2);
    EXPECT_EQ(s.nn[2], 0);
}

TEST_F(StructureTest, UpdateNn3)
{
    Structure s {msgr};
    s.add_disconnected_segment(4);
    s.add_disconnected_segment(3);

    s.update_nn<3>();
    EXPECT_EQ(s.nn[0], 0);
    EXPECT_EQ(s.nn[1], 0);
    EXPECT_EQ(s.nn[2], 0);
}

TEST_F(StructureTest, UpdateNodeNumbers)
{
    Structure s {msgr};
    s.add_disconnected_segment(4);
    s.add_disconnected_segment(3);

    s.update_node_numbers();
    EXPECT_EQ(s.nn[0], 2 * s.mtnum);
    EXPECT_EQ(s.nn[1], s.mt[1].g.size() + s.mt[2].g.size() - 2);
    EXPECT_EQ(s.nn[2], 0);
}

TEST_F(StructureTest, MakeIndma)
{
    constexpr std::array<szt,2> len {4, 3};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    Structure s {msgr};

    for (const auto u : len)
        s.add_disconnected_segment(u);

    s.make_indma();

    EXPECT_EQ(s.mtnum, 2);
    EXPECT_EQ(s.clnum, 2);
    EXPECT_EQ(s.mtmass, lensum);

    EXPECT_EQ(s.cls.size(), 2);
    for (szt j=0; j<s.clnum; j++)
        EXPECT_EQ(s.cls[j], len[j]);

    EXPECT_EQ(s.glm.size(), s.mtmass);
    EXPECT_EQ(s.gla.size(), s.mtmass);
    for (szt j=0; j<len[0]; j++) {
        EXPECT_EQ(s.glm[j], 1);
        EXPECT_EQ(s.gla[j], j);
    }
    for (szt j=len[0]; j<len[1]; j++) {
        EXPECT_EQ(s.glm[j], 2);
        EXPECT_EQ(s.gla[j], j - len[0]);
    }
}


TEST_F(StructureTest, PopulateClusterVectors)
{
    constexpr std::array<szt,2> len {4, 3};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    Structure s {msgr};
    for (const auto u : len)
        s.add_disconnected_segment(u);

    s.populate_cluster_vectors();

    EXPECT_EQ(s.nn[0], 2 * len.size());
    EXPECT_EQ(s.nn[1], lensum - len.size());
    EXPECT_EQ(s.nn[2], 0);

    ASSERT_EQ(s.mtc11.size(), len.size());
    ASSERT_EQ(s.mtc22.size(), len.size());
    ASSERT_EQ(s.mtc33.size(), len.size());
    ASSERT_EQ(s.mtc13.size(), len.size());

    ASSERT_EQ(s.mt11.size(), 2);
    EXPECT_EQ(s.mt11[0], 1);
    EXPECT_EQ(s.mt11[1], 2);
    ASSERT_TRUE(s.mt22.empty());
    ASSERT_TRUE(s.mt33.empty());
    ASSERT_TRUE(s.mt13.empty());

    ASSERT_EQ(s.clmt.size(), 2);
    ASSERT_EQ(s.clmt[0].size(), 1);
    EXPECT_EQ(s.clmt[0][0], 1);
    ASSERT_EQ(s.clmt[1].size(), 1);
    EXPECT_EQ(s.clmt[1][0], 2);
}

TEST_F(StructureTest, MakeAJL)
{
    constexpr std::array<szt,2> len {4, 5};

    Structure s {msgr};
    for (const auto u : len)
        s.add_disconnected_segment(u);

    s.make_indma();                 // to initialize 'cls'
    s.populate_cluster_vectors();   // to initialize 'clmt'
    
    utils::vec2<szt> ajl;
    s.make_adjacency_list_edges(1, ajl);

    EXPECT_EQ(ajl.size(), len[1]);
    EXPECT_EQ(ajl[0].size(), 1);
    EXPECT_EQ(ajl[0][0], 1);
    for (szt i=1; i<ajl.size()-2; i++) {
        EXPECT_EQ(ajl[i].size(), 2);
        EXPECT_EQ(ajl[i][0], i-1);
        EXPECT_EQ(ajl[i][1], i+1);
    }
    EXPECT_EQ(ajl.back().size(), 1);
    EXPECT_EQ(ajl.back()[0], len[1]-2);

}

}  // namespace structure_test
