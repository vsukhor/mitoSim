#include "gtest/gtest.h"

#include "../definitions.h"
#include "../segment.h"

namespace segment_test {

class SegmentTest
    : public testing::Test {

protected:

    using Msgr = mitosim::Msgr;
    using Segment = mitosim::Segment<3>;
    using szt = mitosim::szt;

    struct Config {

        static constexpr szt segmass = 4;
        static constexpr szt cl = 34;
        static constexpr szt mtmass0 = 3;
        static constexpr szt ei0 = 8;

    };

    SegmentTest()
        : msgr {}
        , conf {}
    {}

//    void SetUp() override {
//    }

    static constexpr auto maxDegree = Segment::maxDegree;

    Msgr msgr;
    Config conf;
};

TEST_F(SegmentTest, Constructor1)
{
    const Segment sg {msgr};

    EXPECT_EQ(sg.neig[1].size(), maxDegree);
    EXPECT_EQ(sg.neig[2].size(), maxDegree);
    EXPECT_EQ(sg.neen[1].size(), maxDegree);
    EXPECT_EQ(sg.neen[2].size(), maxDegree);
    ASSERT_TRUE(sg.g.empty());
    EXPECT_EQ(sg.nn[0], 0);
    EXPECT_EQ(sg.nn[1], 0);
    EXPECT_EQ(sg.nn[1], 0);
    EXPECT_EQ(sg.get_cl(), 0);
}

TEST_F(SegmentTest, Constructor2)
{
    const Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};

    EXPECT_EQ(sg.neig[1].size(), maxDegree);
    EXPECT_EQ(sg.neig[2].size(), maxDegree);
    EXPECT_EQ(sg.neen[1].size(), maxDegree);
    EXPECT_EQ(sg.neen[2].size(), maxDegree);
    EXPECT_EQ(sg.g.size(), conf.segmass);
    EXPECT_EQ(sg.nn[0], 0);
    EXPECT_EQ(sg.nn[1], 0);
    for (szt i = 1; i<sg.g.size(); i++) {
        EXPECT_EQ(sg.g[i-1].get_ind() + 1, sg.g[i].get_ind());
        EXPECT_EQ(sg.g[i-1].get_indcl() + 1, sg.g[i].get_indcl());
        EXPECT_EQ(sg.g[i-1].get_cl(), sg.g[i].get_cl());
    }
    for (szt i = 1; i<sg.g.size(); i++) {
        EXPECT_EQ(sg.g[i].get_fin(0), 0.);
        EXPECT_EQ(sg.g[i].get_fin(1), 0.);
    }
}

TEST_F(SegmentTest, reflectG)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};
    sg.reflect_g();

    for (szt i = 0; i < sg.g.size(); i++) {
        EXPECT_EQ(sg.g[i].get_ind(), conf.ei0 + Config::segmass - i - 1);
        EXPECT_EQ(sg.g[i].get_indcl(), sg.g.size() - i - 1);
        EXPECT_EQ(sg.g[i].get_cl(), Config::cl);
    }
}

TEST_F(SegmentTest, setGCl)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};
    const auto newCl = Config::cl + 100;
    const auto newIndcl = conf.ei0 + 100;

    const auto res = sg.set_gCl(newCl, newIndcl);

    for (szt i = 0; i < sg.g.size(); i++) {
        EXPECT_EQ(sg.g[i].get_indcl(), newIndcl + i);
        EXPECT_EQ(sg.g[i].get_ind(), conf.ei0 + i);
        EXPECT_EQ(sg.g[i].get_cl(), newCl);
        EXPECT_EQ(sg.g[i].get_fin(0), 0.);
        EXPECT_EQ(sg.g[i].get_fin(1), 0.);
        EXPECT_EQ(res, newIndcl + sg.g.size());
    }
}

TEST_F(SegmentTest, setCl)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};
    const auto newCl = Config::cl + 100;
    const auto newIndcl = conf.ei0 + 100;

    const auto res = sg.setCl(newCl, newIndcl);

    EXPECT_EQ(sg.get_cl(), newCl);
    EXPECT_EQ(res, sg.set_gCl(newCl, newIndcl));
}

TEST_F(SegmentTest, End2A)
{
    const Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};

    EXPECT_EQ(sg.end2a(1), 0);
    EXPECT_EQ(sg.end2a(2), sg.g.size() - 1);
//    ASSERT_DEATH(sg.end2a(0), ::testing::Eq("Incorrect end index."));
}

TEST_F(SegmentTest, HasOneFreeEnd)
{
    const Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};

    EXPECT_FALSE(sg.has_one_free_end());
}

TEST_F(SegmentTest, SingleNeigIndex)
{
    const Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};

    EXPECT_EQ(sg.single_neig_index(1), mitosim::huge<szt>);
    EXPECT_EQ(sg.single_neig_index(2), mitosim::huge<szt>);
}

TEST_F(SegmentTest, DoubleNeigIndexes)
{
    const Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};

    ASSERT_EXIT(sg.double_neig_indexes(1),
                testing::ExitedWithCode(EXIT_FAILURE), "");
}

TEST_F(SegmentTest, IsCycle)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};

    EXPECT_FALSE(sg.is_cycle());
}

TEST_F(SegmentTest, NumNodes)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};

    EXPECT_EQ(sg.num_nodes(1), 2);
    EXPECT_EQ(sg.num_nodes(2), Config::segmass - 1);
    EXPECT_EQ(sg.num_nodes(3), 0);
//  ASSERT_DEATH(sg.num_nodes(4), ::testing::Eq("Error in Segment::num_nodes(). Not implemented for degree 4."));
}

TEST_F(SegmentTest, SetEndFin1)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};
    sg.template set_end_fin<1>();

    for (const auto& g : sg.g) {
        EXPECT_EQ(g.get_fin(0), 0.);
        EXPECT_EQ(g.get_fin(1), 0.);
    }
}

TEST_F(SegmentTest, SetEndFin2)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};
    sg.template set_end_fin<2>();

    for (const auto& g : sg.g) {
        EXPECT_EQ(g.get_fin(0), 0.);
        EXPECT_EQ(g.get_fin(1), 0.);
    }
}

TEST_F(SegmentTest, SetBulkFin)
{
    Segment sg {Config::segmass, Config::cl, conf.ei0, msgr};
    const szt a = 100;
    sg.set_bulk_fin(a);

    EXPECT_EQ(sg.g[0].get_fin(0), 0.);
    for (szt i = 0; i < sg.g.size()-1; i++) {
        if (i == a) {
            EXPECT_EQ(sg.g[i].get_fin(1), 1.);
            EXPECT_EQ(sg.g[i+1].get_fin(0), 1.);
        }
        else {
            EXPECT_EQ(sg.g[i].get_fin(1), 0.);
            EXPECT_EQ(sg.g[i+1].get_fin(0), 0.);
        }
    }
    EXPECT_EQ(sg.g.back().get_fin(1), 0.);
}

}  // namespace segment_test
