#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"
#include "../segment.h"
#include "../ability_for_fusion.h"

namespace ability_fission_test {

// Subclass to make protected members accessible for testing:
// Use AbilityForFusion rather than AbilityForFission to be able to produce
// branched fusion constracts needed for testing 'fiss3'.
class AF
    : public mitosim::AbilityForFusion<mitosim::Segment<3>> {

public:

    using Msgr = mitosim::Msgr;
    using Mt = mitosim::Segment<3>;
    using AbilityForFusion = mitosim::AbilityForFusion<Mt>;

    using AbilityForFusion::clnum;
    using AbilityForFusion::clagl;
    using AbilityForFusion::copy_neigs;
    using AbilityForFusion::fuse_antiparallel;
    using AbilityForFusion::fuse_parallel;
    using AbilityForFusion::fuse_to_loop;
    using AbilityForFusion::gla;
    using AbilityForFusion::glm;
    using AbilityForFusion::mt;
    using AbilityForFusion::mtnum;
    using AbilityForFusion::rename_mito;
    using AbilityForFusion::update_neigs;
    using AbilityForFission::fiss2;
    using AbilityForFission::fiss3;

    AF(Msgr* msgr) : AbilityForFusion {*msgr} {}
};


class AbilityFissionTest
    : public testing::Test {

protected:

    using Msgr = mitosim::Msgr;
    using Mt = mitosim::Segment<3>;
    using szt = mitosim::szt;

    AbilityFissionTest()
        : msgr {&std::cout, nullptr, 6}
    {}

    mitosim::Msgr msgr;

};

TEST_F(AbilityFissionTest, Constructor)
{
    AF ct {&msgr};

    ASSERT_TRUE(ct.clagl.empty());
    ASSERT_TRUE(ct.glm.empty());
    ASSERT_TRUE(ct.gla.empty());
    ASSERT_TRUE(ct.mt.empty());
    for (const auto& n : ct.nn)
        EXPECT_EQ(n, 0);
    EXPECT_EQ(ct.mtnum, 0);
    EXPECT_EQ(ct.clnum, 0);
    EXPECT_EQ(ct.mtmass, 0);
    ASSERT_TRUE(ct.clmt.empty());
    ASSERT_TRUE(ct.mt11.empty());
    ASSERT_TRUE(ct.mtc11.empty());
    ASSERT_TRUE(ct.mt22.empty());
    ASSERT_TRUE(ct.mtc22.empty());
    ASSERT_TRUE(ct.mt33.empty());
    ASSERT_TRUE(ct.mtc33.empty());
    ASSERT_TRUE(ct.mt13.empty());
    ASSERT_TRUE(ct.mtc13.empty());
}

TEST_F(AbilityFissionTest, Fiss2)
{
    // Tests fission of a separate linear segment.

    constexpr szt len {4};

    constexpr szt w1 = 1;
    constexpr szt w2 = 2;

    for (szt a=1; a<len; a++) {
        AF ct {&msgr};
        ct.add_disconnected_segment(len);

        ct.fiss2(w1, a);

        ASSERT_EQ(ct.mtmass, len);
        ASSERT_EQ(ct.mtnum, 2);
        ASSERT_EQ(ct.clnum, 2);

        ct.update_node_numbers();
        ASSERT_EQ(ct.nn[0], 4);
        ASSERT_EQ(ct.nn[1], len - 2);
        ASSERT_EQ(ct.nn[2], 0);

        ASSERT_EQ(ct.mt[w1].g.size(), a);
        ASSERT_EQ(ct.mt[w2].g.size(), len - a);

        for (szt c=0, j=1; j<=ct.mtnum; j++) {
            const auto& m = ct.mt[j];
            ASSERT_EQ(m.get_cl(), j - 1);
            for (szt i=0; i<m.g.size(); i++) {
                ASSERT_EQ(m.g[i].get_cl(), m.get_cl());
                ASSERT_EQ(m.g[i].get_indcl(), i);
                ASSERT_EQ(m.g[i].get_ind(), c++);
            }
            ASSERT_EQ(m.nn[1], 0);
            ASSERT_EQ(m.nn[2], 0);
        }
    }
}

TEST_F(AbilityFissionTest, Fiss3a)
{
    // fiss3 of a 13 segment at end 2.

    constexpr std::array<szt,2> len {4, 4};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    constexpr szt a = 2;

    AF ct {&msgr};
    for (const auto u : len)
        ct.add_disconnected_segment(u);

    ct.fuse12(1, 2, 2, a);

    ct.fiss3(1, 2);

    ASSERT_EQ(ct.mtmass, lensum);
    ASSERT_EQ(ct.mtnum, 2);
    ASSERT_EQ(ct.clnum, 2);

    ct.update_node_numbers();
    ASSERT_EQ(ct.nn[0], 4);
    ASSERT_EQ(ct.nn[1], lensum - 2);
    ASSERT_EQ(ct.nn[2], 0);

    ASSERT_EQ(ct.mt[1].g.size(), len[0]);
    ASSERT_EQ(ct.mt[2].g.size(), len[1]);

    for (szt c=0, j=1; j<=ct.mtnum; j++) {
        const auto& m = ct.mt[j];
        ASSERT_EQ(m.get_cl(), j - 1);
        for (szt i=0; i<m.g.size(); i++) {
            ASSERT_EQ(m.g[i].get_cl(), m.get_cl());
            ASSERT_EQ(m.g[i].get_indcl(), i);
            ASSERT_EQ(m.g[i].get_ind(), c++);
        }
        ASSERT_EQ(m.nn[1], 0);
        ASSERT_EQ(m.nn[2], 0);
    }
}

TEST_F(AbilityFissionTest, Fiss3b)
{
    // fiss3 of a 13 segment at end 1.

    constexpr std::array<szt,2> len {4, 4};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    constexpr szt a = 2;

        AF ct {&msgr};
    for (const auto u : len)
        ct.add_disconnected_segment(u);


    ct.fuse12(1, 2, 2, a);

    ct.fiss3(3, 1);

    ASSERT_EQ(ct.mtmass, lensum);
    ASSERT_EQ(ct.mtnum, 2);
    ASSERT_EQ(ct.clnum, 2);

    ct.update_node_numbers();
    ASSERT_EQ(ct.nn[0], 4);
    ASSERT_EQ(ct.nn[1], lensum - 2);
    ASSERT_EQ(ct.nn[2], 0);

    ASSERT_EQ(ct.mt[1].g.size(), len[0] + a);
    ASSERT_EQ(ct.mt[2].g.size(), len[1] - a);

    ASSERT_EQ(ct.mt[1].get_cl(), 1);
    ASSERT_EQ(ct.mt[2].get_cl(), 0);

    for (szt c=0, j=ct.mtnum; j<0; j--) {
        const auto& m = ct.mt[j];
        for (szt i=0; i<m.g.size(); i++) {
            ASSERT_EQ(m.g[i].get_cl(), m.get_cl());
            ASSERT_EQ(m.g[i].get_indcl(), i);
            ASSERT_EQ(m.g[i].get_ind(), c++);
        }
        ASSERT_EQ(m.nn[1], 0);
        ASSERT_EQ(m.nn[2], 0);
    }
}

TEST_F(AbilityFissionTest, Fiss3c)
{
    // Tests fission of circular segment at end 1.
    constexpr szt len = 10;

    constexpr szt w = 1;
    const szt e = 1;
    constexpr szt a = 6;

    AF ct {&msgr};
    ct.add_disconnected_segment(len);

    ct.fuse12(w, e, w, a);

    ct.fiss3(w, 1);

    const auto& m = ct.mt[1];

    ASSERT_EQ(m.g.size(), len);
    ASSERT_EQ(m.get_cl(), 0);
    for (szt i=0; i<m.g.size(); i++) {
        ASSERT_EQ(m.g[i].get_cl(), m.get_cl());
        ASSERT_EQ(m.g[i].get_indcl(), i);
        ASSERT_EQ(m.g[i].get_ind(), i);
    }
    ASSERT_EQ(m.nn[1], 0);
    ASSERT_EQ(m.nn[2], 0);

}

TEST_F(AbilityFissionTest, Fiss3d)
{
    // Tests fission of circular segment at end 2.
    constexpr szt len = 10;

    constexpr szt w = 1;
    const szt e = 1;
    constexpr szt a = 6;

    AF ct {&msgr};
    ct.add_disconnected_segment(len);

    ct.fuse12(w, e, w, a);

    ct.fiss3(w, 2);

    const auto& m = ct.mt[1];

    ASSERT_EQ(m.g.size(), len);
    ASSERT_EQ(m.get_cl(), 0);
    for (szt i=0; i<m.g.size(); i++) {
        ASSERT_EQ(m.g[i].get_cl(), m.get_cl());
        ASSERT_EQ(m.g[i].get_indcl(), i);
        if (i < a)
            ASSERT_EQ(m.g[i].get_ind(), a - 1 - i);
        else
            ASSERT_EQ(m.g[i].get_ind(), i);
    }
    ASSERT_EQ(m.nn[1], 0);
    ASSERT_EQ(m.nn[2], 0);

}

TEST_F(AbilityFissionTest, Fiss3e)
{
    // Tests fission of linear segment at end 1 from a circular one.
    constexpr szt len = 10;

    constexpr szt w = 1;
    constexpr szt v = 2;    // segment resulting from the fusion.
    const szt e = 1;
    constexpr szt a = 6;

    AF ct {&msgr};
    ct.add_disconnected_segment(len);

    ct.fuse12(w, e, w, a);

    ct.fiss3(v, 1);

    const auto& m = ct.mt[1];
    const auto& n = ct.mt[2];

    ASSERT_EQ(m.g.size(), a);
    ASSERT_EQ(n.g.size(), len - a);

    ASSERT_EQ(m.nn[1], 1);
    ASSERT_EQ(m.neig[1][1], 1);
    ASSERT_EQ(m.neen[1][1], 2);
    ASSERT_EQ(m.nn[2], 1);
    ASSERT_EQ(m.neig[2][1], 1);
    ASSERT_EQ(m.neen[2][1], 1);

    ASSERT_EQ(n.nn[1], 0);
    ASSERT_EQ(n.nn[2], 0);

    for (szt c=0, j=1; j<=ct.mtnum; j++) {
        const auto& x = ct.mt[j];
        ASSERT_EQ(x.get_cl(), ct.clnum - j);
        for (szt i=0; i<x.g.size(); i++) {
            ASSERT_EQ(x.g[i].get_cl(), x.get_cl());
            ASSERT_EQ(x.g[i].get_indcl(), i);
            ASSERT_EQ(x.g[i].get_ind(), c++);
        }
    }
}

}  // namespace ability_fission_test
