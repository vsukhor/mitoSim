#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"
#include "../segment.h"
#include "../ability_for_fusion.h"

namespace ability_fusion_test {

// Subclass to make protected members accessible for testing:
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

    AF(Msgr* msgr) : AbilityForFusion {*msgr} {}
};

class AbilityFusionTest
    : public ::testing::Test {

protected:

    using Msgr = mitosim::Msgr;
    using Mt = mitosim::Segment<3>;
    using szt = mitosim::szt;

    void SetUp() override {
        msgr = new Msgr{&std::cout, nullptr, 6};
    }

    void TearDown() override {
        delete msgr;
    }

    mitosim::Msgr* msgr;

};

TEST_F(AbilityFusionTest, Constructor)
{
    AF ct {msgr};

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

TEST_F(AbilityFusionTest, Fuse12a)
{
    // Tests fusion to a separate linear segment.
    constexpr std::array<szt,2> len {{4, 4}};

    constexpr szt w1 = 2;
    constexpr szt w2 = 1;
    constexpr szt e = 2;
    constexpr szt oe = e == 1 ? 2 : 1;

    for (szt a=1; a<len[w2-1]; a++) {
        AF ct {msgr};
        for (const auto u : len)
            ct.add_disconnected_segment(u);

        ct.fuse12(w1, e, w2, a);

        const auto clmass = len[w1-1] + len[w2-1];
        ASSERT_EQ(ct.mtmass, clmass);
        ASSERT_EQ(ct.mtnum, 3);
        ASSERT_EQ(ct.clnum, 1);

        ct.update_node_numbers();
        ASSERT_EQ(ct.nn[0], 3);
        ASSERT_EQ(ct.nn[1], clmass - 3);
        ASSERT_EQ(ct.nn[2], 1);

        for (szt i=0; i<ct.mt.size(); i++) {
            ASSERT_EQ(ct.mt[i].get_cl(), 0);
            for (const auto& g : ct.mt[i].g)
                ASSERT_EQ(g.get_cl(), ct.mt[i].get_cl());
        }

        szt indcl {};
        ASSERT_EQ(ct.mt[w2].g.size(), a);
        for (szt i=0; i<ct.mt[w2].g.size(); i++) {
            ASSERT_EQ(ct.mt[w2].g[i].get_indcl(), indcl++);
            ASSERT_EQ(ct.mt[w2].g[i].get_ind(), i);
        }
        ASSERT_EQ(ct.mt[w1].g.size(), len[w1-1]);
        for (szt i=0; i<ct.mt[w1].g.size(); i++) {
            ASSERT_EQ(ct.mt[w1].g[i].get_indcl(), indcl++);
            ASSERT_EQ(ct.mt[w1].g[i].get_ind(),
                      i + ct.mt[w2].g.size() + ct.mt.back().g.size());
        }
        ASSERT_EQ(ct.mt.back().g.size(), len[w2-1] - a);
        for (szt i=0; i<ct.mt.back().g.size(); i++) {
            ASSERT_EQ(ct.mt.back().g[i].get_indcl(), indcl++);
            ASSERT_EQ(ct.mt.back().g[i].get_ind(),
                      i + ct.mt[w2].g.size());
        }

        ASSERT_EQ(ct.mt[w1].nn[oe], 0);
        ASSERT_EQ(ct.mt[w1].nn[e], 2);
        ASSERT_EQ(ct.mt[w1].neig[e][1], w2);
        ASSERT_EQ(ct.mt[w1].neen[e][1], 2);
        ASSERT_EQ(ct.mt[w1].neig[e][2], ct.mt.size()-1);
        ASSERT_EQ(ct.mt[w1].neen[e][2], 1);

        ASSERT_EQ(ct.mt[w2].nn[1], 0);
        ASSERT_EQ(ct.mt[w2].nn[2], 2);
        ASSERT_EQ(ct.mt[w2].neig[2][1], w1);
        ASSERT_EQ(ct.mt[w2].neen[2][1], e);
        ASSERT_EQ(ct.mt[w2].neig[2][2], ct.mt.size()-1);
        ASSERT_EQ(ct.mt[w2].neen[2][2], 1);

        ASSERT_EQ(ct.mt.back().nn[1], 2);
        ASSERT_EQ(ct.mt.back().neig[1][1], w1);
        ASSERT_EQ(ct.mt.back().neen[1][1], 2);
        ASSERT_EQ(ct.mt.back().neig[1][2], w2);
        ASSERT_EQ(ct.mt.back().neen[1][2], e);
        ASSERT_EQ(ct.mt.back().nn[2], 0);
    }
}

TEST_F(AbilityFusionTest, Fuse12b)
{
    // Tests fusion to itself.
    constexpr szt len = 10;

    constexpr szt w = 1;
    constexpr szt v = 2;    // segment resulting from the fusion.
    for (const szt e : {1, 2}) {
        const szt oe = e == 1 ? 2 : 1;
        for (szt a=1; a<len; a++) {
            AF ct {msgr};
            ct.add_disconnected_segment(len);

            ct.fuse12(w, e, w, a);

            const auto clmass = len;
            ASSERT_EQ(ct.mtmass, clmass);
            ASSERT_EQ(ct.mtnum, 2);
            ASSERT_EQ(ct.clnum, 1);

            ct.update_node_numbers();
            ASSERT_EQ(ct.nn[0], 1);
            ASSERT_EQ(ct.nn[1], clmass - 2);
            ASSERT_EQ(ct.nn[2], 1);

            for (szt i=0; i<ct.mt.size(); i++) {
                ASSERT_EQ(ct.mt[i].get_cl(), 0);
                for (const auto& g : ct.mt[i].g)
                    ASSERT_EQ(g.get_cl(), ct.mt[i].get_cl());
            }

            szt indcl {};
            ASSERT_EQ(ct.mt[w].g.size(), a);
            for (szt i=0; i<ct.mt[w].g.size(); i++) {
                ASSERT_EQ(ct.mt[w].g[i].get_indcl(), indcl++);
                ASSERT_EQ(ct.mt[w].g[i].get_ind(), i);
            }
            ASSERT_EQ(ct.mt[v].g.size(), len - a);
            for (szt i=0; i<ct.mt[v].g.size(); i++) {
                ASSERT_EQ(ct.mt[v].g[i].get_indcl(), indcl++);
                ASSERT_EQ(ct.mt[v].g[i].get_ind(),
                          i + ct.mt[w].g.size());
            }
            if (e == 1) {
                ASSERT_EQ(ct.mt[w].nn[oe], 2);
                ASSERT_EQ(ct.mt[w].neig[oe][1], w);
                ASSERT_EQ(ct.mt[w].neen[oe][1], e);
                ASSERT_EQ(ct.mt[w].neig[oe][2], v);
                ASSERT_EQ(ct.mt[w].neen[oe][2], e);
                ASSERT_EQ(ct.mt[w].nn[e], 2);
                ASSERT_EQ(ct.mt[w].neig[e][1], w);
                ASSERT_EQ(ct.mt[w].neen[e][1], oe);
                ASSERT_EQ(ct.mt[w].neig[e][2], v);
                ASSERT_EQ(ct.mt[w].neen[e][2], e);

                ASSERT_EQ(ct.mt[v].nn[oe], 0);
                ASSERT_EQ(ct.mt[v].nn[e], 2);
                ASSERT_EQ(ct.mt[v].neig[e][1], w);
                ASSERT_EQ(ct.mt[v].neen[e][1], e);
                ASSERT_EQ(ct.mt[v].neig[e][2], w);
                ASSERT_EQ(ct.mt[v].neen[e][2], oe);
            }
            else {
                ASSERT_EQ(ct.mt[v].nn[oe], 2);
                ASSERT_EQ(ct.mt[v].neig[oe][1], w);
                ASSERT_EQ(ct.mt[v].neen[oe][1], e);
                ASSERT_EQ(ct.mt[v].neig[oe][2], v);
                ASSERT_EQ(ct.mt[v].neen[oe][2], e);
                ASSERT_EQ(ct.mt[v].nn[e], 2);
                ASSERT_EQ(ct.mt[v].neig[e][1], w);
                ASSERT_EQ(ct.mt[v].neen[e][1], e);
                ASSERT_EQ(ct.mt[v].neig[e][2], v);
                ASSERT_EQ(ct.mt[v].neen[e][2], oe);

                ASSERT_EQ(ct.mt[w].nn[oe], 0);
                ASSERT_EQ(ct.mt[w].nn[e], 2);
                ASSERT_EQ(ct.mt[w].neig[e][1], v);
                ASSERT_EQ(ct.mt[w].neen[e][1], oe);
                ASSERT_EQ(ct.mt[w].neig[e][2], v);
                ASSERT_EQ(ct.mt[w].neen[e][2], e);
            }
        }
    }
}

TEST_F(AbilityFusionTest, Fuse12c)
{
    // Fusion to a cycle segment.
    constexpr std::array<szt,2> len {{4, 4}};

    constexpr szt w1 = 2;
    constexpr szt w2 = 1;
    constexpr szt e = 2;
    constexpr szt oe = e == 1 ? 2 : 1;

    for (szt a=1; a<len[w2-1]; a++) {
        AF ct {msgr};
        for (const auto u : len)
            ct.add_disconnected_segment(u);

        ct.fuse_to_loop(w2);
        ct.fuse12(w1, e, w2, a);

        const auto clmass = len[w1-1] + len[w2-1];
        ASSERT_EQ(ct.mtmass, clmass);
        ASSERT_EQ(ct.mtnum, 2);
        ASSERT_EQ(ct.clnum, 1);

        ct.update_node_numbers();
        ASSERT_EQ(ct.nn[0], 1);
        ASSERT_EQ(ct.nn[1], clmass - 2);
        ASSERT_EQ(ct.nn[2], 1);

        ASSERT_EQ(ct.mt[w1].g.size(), len[w1-1]);
        ASSERT_EQ(ct.mt[w2].g.size(), len[w2-1]);
        szt c {};

        std::vector<szt> t(ct.mt[w2].g.size());
        std::iota(t.begin(), t.end(), 0);
        std::rotate(t.begin(), t.begin() + a, t.end());
        for (szt i=0; i<ct.mt[w2].g.size(); i++) {
            ASSERT_EQ(ct.mt[w2].g[i].get_indcl(), c);
            ASSERT_EQ(ct.mt[w2].g[i].get_ind(), t[i]);
            c++;
        }
        for (const auto& g : ct.mt[w1].g) {
            ASSERT_EQ(g.get_indcl(), c);
            ASSERT_EQ(g.get_ind(), c);
            c++;
        }

        ASSERT_EQ(ct.mt[w1].nn[oe], 0);
        ASSERT_EQ(ct.mt[w1].nn[e], 2);
        ASSERT_EQ(ct.mt[w1].neig[e][1], w2);
        ASSERT_EQ(ct.mt[w1].neen[e][1], 2);
        ASSERT_EQ(ct.mt[w1].neig[e][2], w2);
        ASSERT_EQ(ct.mt[w1].neen[e][2], 1);

        ASSERT_EQ(ct.mt[w2].nn[1], 2);
        ASSERT_EQ(ct.mt[w2].nn[2], 2);
        ASSERT_EQ(ct.mt[w2].neig[1][1], w1);
        ASSERT_EQ(ct.mt[w2].neen[1][1], e);
        ASSERT_EQ(ct.mt[w2].neig[1][2], w2);
        ASSERT_EQ(ct.mt[w2].neen[1][2], 2);
        ASSERT_EQ(ct.mt[w2].neig[2][1], w1);
        ASSERT_EQ(ct.mt[w2].neen[2][1], e);
        ASSERT_EQ(ct.mt[w2].neig[2][2], w2);
        ASSERT_EQ(ct.mt[w2].neen[2][2], 1);
    }
}

TEST_F(AbilityFusionTest, Fuse1L)
{
    // Fusion to a cycle segment.
    constexpr std::array<szt,2> len {{4, 4}};

    constexpr szt w1 = 2;
    constexpr szt w2 = 1;
    constexpr szt e = 2;
    constexpr szt oe = e == 1 ? 2 : 1;

    AF ct {msgr};
    for (const auto u : len)
        ct.add_disconnected_segment(u);

    ct.fuse_to_loop(w2);

    ct.fuse1L(w1, e, w2);

    const auto clmass = len[w1-1] + len[w2-1];
    ASSERT_EQ(ct.mtmass, clmass);
    ASSERT_EQ(ct.mtnum, 2);
    ASSERT_EQ(ct.clnum, 1);

    ct.update_node_numbers();
    ASSERT_EQ(ct.nn[0], 1);
    ASSERT_EQ(ct.nn[1], clmass - 2);
    ASSERT_EQ(ct.nn[2], 1);

    ASSERT_EQ(ct.mt[w1].g.size(), len[w1-1]);
    ASSERT_EQ(ct.mt[w2].g.size(), len[w2-1]);
    szt c {};

    std::vector<szt> v(ct.mt[w2].g.size());
    std::iota(v.begin(), v.end(), 0);
    std::rotate(v.begin(), v.begin()+0, v.end());
    for (szt i=0; i<ct.mt[w2].g.size(); i++) {
        ASSERT_EQ(ct.mt[w2].g[i].get_indcl(), c);
        ASSERT_EQ(ct.mt[w2].g[i].get_ind(), v[i]);
        c++;
    }
    for (const auto& g : ct.mt[w1].g) {
        ASSERT_EQ(g.get_indcl(), c);
        ASSERT_EQ(g.get_ind(), c);
        c++;
    }

    ASSERT_EQ(ct.mt[w1].nn[oe], 0);
    ASSERT_EQ(ct.mt[w1].nn[e], 2);
    ASSERT_EQ(ct.mt[w1].neig[e][1], w2);
    ASSERT_EQ(ct.mt[w1].neen[e][1], 1);
    ASSERT_EQ(ct.mt[w1].neig[e][2], w2);
    ASSERT_EQ(ct.mt[w1].neen[e][2], 2);

    ASSERT_EQ(ct.mt[w2].nn[1], 2);
    ASSERT_EQ(ct.mt[w2].nn[2], 2);
    ASSERT_EQ(ct.mt[w2].neig[1][1], w2);
    ASSERT_EQ(ct.mt[w2].neen[1][1], 2);
    ASSERT_EQ(ct.mt[w2].neig[1][2], w1);
    ASSERT_EQ(ct.mt[w2].neen[1][2], e);
    ASSERT_EQ(ct.mt[w2].neig[2][1], w2);
    ASSERT_EQ(ct.mt[w2].neen[2][1], 1);
    ASSERT_EQ(ct.mt[w2].neig[2][2], w1);
    ASSERT_EQ(ct.mt[w2].neen[2][2], e);
}

}  // namespace ability_fusion_test
