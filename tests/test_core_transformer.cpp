#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"
#include "../segment.h"
#include "../core_transformer.h"

namespace core_transformer_test {


// Subclass to make protected members accessible for testing:
class CT
    : public mitosim::CoreTransformer<mitosim::Segment<3>> {

public:

    using CoreTransformer::clnum;
    using CoreTransformer::copy_neigs;
    using CoreTransformer::fuse_antiparallel;
    using CoreTransformer::fuse_parallel;
    using CoreTransformer::fuse_to_loop;
    using CoreTransformer::mt;
    using CoreTransformer::mtnum;
    using CoreTransformer::rename_mito;
    using CoreTransformer::update_neigs;
    using Msgr = mitosim::Msgr;

    CT(Msgr& msgr)
        : CoreTransformer {msgr}
    {}
};


class CoreTransformerTest
    : public testing::Test {

protected:

    using Msgr = mitosim::Msgr;
    using Mt = mitosim::Segment<3>;
    using szt = mitosim::szt;

    CoreTransformerTest()
        : msgr {}
    {}

    mitosim::Msgr msgr;

};

TEST_F(CoreTransformerTest, Constructor)
{
    CT ct {msgr};

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

TEST_F(CoreTransformerTest, FuseAntiparE1)
{
    constexpr std::array<szt,6> len {{4, 3, 6, 5, 2, 7}};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    szt w1=2;
    szt w2=3;
//    for (szt w1=1; w1<=len.size(); w1++)
//        for (szt w2=1; w2<=len.size(); w2++)
            if (w1 != w2) {

                CT ct {msgr};
                for (const auto u : len)
                    ct.add_disconnected_segment(u);
            }
}

TEST_F(CoreTransformerTest, FuseAntiparE2)
{
    constexpr std::array<szt,6> len {{4, 3, 6, 5, 2, 7}};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    szt w1=2;
    szt w2=3;
//    for (szt w1=1; w1<=len.size(); w1++)
//        for (szt w2=1; w2<=len.size(); w2++)
            if (w1 != w2) {

                CT ct {msgr};
                for (const auto u : len)
                    ct.add_disconnected_segment(u);

                for (szt c0=0, i=1; i<=ct.mtnum; i++)
                    for (const auto& g : ct.mt[i].g)
                        ASSERT_EQ(g.get_ind(), c0++);

                decltype(ct.mt[w1].g) g1;
                std::copy(ct.mt[w1].g.begin(),
                          ct.mt[w1].g.end(), std::back_inserter(g1));
                decltype(ct.mt[w2].g) g2;
                std::copy(ct.mt[w2].g.begin(),
                          ct.mt[w2].g.end(), std::back_inserter(g2));

                ct.fuse_antiparallel(2, w1, w2);

                ASSERT_EQ(ct.mtmass, lensum);
                ASSERT_EQ(ct.mtnum, len.size()-1);
                ASSERT_EQ(ct.clnum, len.size()-1);

                // 'res' is the fused segment:
                const auto res = w1 != len.size() ? w1 : w2;

                for (szt i=1; i<=ct.mtnum; i++) {
                    ASSERT_EQ(ct.mt[i].get_cl(), i-1);

                    if (i == res)
                        ASSERT_EQ(ct.mt[i].g.size(), len[w1-1]+len[w2-1]);
                    else if (w1 != len.size() && i == w2)
                        // w2 is occupied by the previously last segment
                        ASSERT_EQ(ct.mt[i].g.size(), len[ct.mtnum]);
                    else    // not affected
                        ASSERT_EQ(ct.mt[i].g.size(), len[i-1]);

                    for (const auto& gg : ct.mt[i].g)
                        ASSERT_EQ(gg.get_cl(), ct.mt[i].get_cl());
                }

                int c {};
                for (szt i=0; i<len[w1-1]; c++, i++)
                    ASSERT_EQ(ct.mt[res].g[i].get_ind(), g1[c].get_ind());
                c = len[w2-1] - 1;
                for (szt i=len[w1-1]; c>=0; c--, i++)
                    ASSERT_EQ(ct.mt[res].g[i].get_ind(), g2[c].get_ind());

                for (const auto i : {1, 2}) {
                    ASSERT_EQ(ct.mt[res].nn[1], 0);
                    for (const auto& n : {ct.mt[res].neig, ct.mt[res].neen}) {
                        ASSERT_EQ(n[i].size(), Mt::maxDegree);
                        for (const auto& nn : n[i])
                            ASSERT_EQ(nn, 0);
                    }
                }
            }
}

TEST_F(CoreTransformerTest, FuseParallel)
{
    constexpr std::array<szt,6> len {{4, 3, 6, 5, 2, 7}};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    szt w1=2;
    szt w2=3;
//    for (szt w1=1; w1<=len.size(); w1++)
//        for (szt w2=1; w2<=len.size(); w2++)
            if (w1 != w2) {

                CT ct {msgr};
                for (const auto u : len)
                    ct.add_disconnected_segment(u);

                for (szt c0=0, i=1; i<=ct.mtnum; i++)
                    for (const auto& g : ct.mt[i].g)
                        ASSERT_EQ(g.get_ind(), c0++);

                decltype(ct.mt[w1].g) g1;
                std::copy(ct.mt[w1].g.begin(),
                          ct.mt[w1].g.end(), std::back_inserter(g1));
                decltype(ct.mt[w2].g) g2;
                std::copy(ct.mt[w2].g.begin(),
                          ct.mt[w2].g.end(), std::back_inserter(g2));

                ct.fuse_parallel(w1, w2);

                ASSERT_EQ(ct.mtmass, lensum);
                ASSERT_EQ(ct.mtnum, len.size()-1);
                ASSERT_EQ(ct.clnum, len.size()-1);

                // 'res' is the fused segment:
                const auto res = w1 != len.size() ? w1 : w2;

                for (szt i=1; i<=ct.mtnum; i++) {
                    ASSERT_EQ(ct.mt[i].get_cl(), i-1);

                    if (i == res)
                        ASSERT_EQ(ct.mt[i].g.size(), len[w1-1]+len[w2-1]);
                    else if (w1 != len.size() && i == w2)
                        // w2 is occupied by the previously last segment
                        ASSERT_EQ(ct.mt[i].g.size(), len[ct.mtnum]);
                    else    // not affected
                        ASSERT_EQ(ct.mt[i].g.size(), len[i-1]);

                    for (const auto& gg : ct.mt[i].g)
                        ASSERT_EQ(gg.get_cl(), ct.mt[i].get_cl());
                }

                for (szt c=0, i=0; i<len[w2-1]; c++, i++)
                    ASSERT_EQ(ct.mt[res].g[i].get_ind(), g2[c].get_ind());
                for (szt c=0, i=len[w2-1]; i<ct.mt[res].g.size(); c++, i++)
                    ASSERT_EQ(ct.mt[res].g[i].get_ind(), g1[c].get_ind());

                for (const auto i : {1, 2}) {
                    ASSERT_EQ(ct.mt[res].nn[1], 0);
                    for (const auto& n : {ct.mt[res].neig, ct.mt[res].neen}) {
                        ASSERT_EQ(n[i].size(), Mt::maxDegree);
                        for (const auto& nn : n[i])
                            ASSERT_EQ(nn, 0);
                    }
                }
            }
}

TEST_F(CoreTransformerTest, FuseToLoop)
{
    constexpr std::array<szt,6> len {{4, 3, 6, 5, 2, 7}};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    for (szt w=1; w<=len.size(); w++) {
        CT ct {msgr};
        for (const auto u : len)
            ct.add_disconnected_segment(u);
        const auto& m = ct.mt;

        for (szt c0=0, i=1; i<=ct.mtnum; i++)
            for (const auto& g : m[i].g)
                ASSERT_EQ(g.get_ind(), c0++);

        decltype(m[w].g) g;
        std::copy(m[w].g.begin(),
                  m[w].g.end(), std::back_inserter(g));

        ct.fuse_to_loop(w);

        ASSERT_EQ(ct.mtmass, lensum);
        ASSERT_EQ(ct.mtnum, len.size());
        ASSERT_EQ(ct.clnum, len.size());

        for (szt i=1; i<=ct.mtnum; i++) {
            ASSERT_EQ(m[i].g.size(), len[i-1]);
            ASSERT_EQ(m[i].get_cl(), i-1);
        }
        for (szt i=0; i<len[w-1]; i++) {
            ASSERT_EQ(m[w].g[i].get_ind(), g[i].get_ind());
            ASSERT_EQ(m[w].g[i].get_cl(), m[w].get_cl());
        }

        ASSERT_EQ(m[w].nn[1], 1);
        ASSERT_EQ(m[w].nn[2], 1);
        ASSERT_EQ(m[w].neig[1].size(), Mt::maxDegree);
        ASSERT_EQ(m[w].neig[1][0], 0);
        ASSERT_EQ(m[w].neig[1][1], w);
        ASSERT_EQ(m[w].neig[1][2], 0);
        ASSERT_EQ(m[w].neig[2].size(), Mt::maxDegree);
        ASSERT_EQ(m[w].neig[2][0], 0);
        ASSERT_EQ(m[w].neig[2][1], w);
        ASSERT_EQ(m[w].neig[2][2], 0);
        ASSERT_EQ(m[w].neen[1].size(), Mt::maxDegree);
        ASSERT_EQ(m[w].neen[1][0], 0);
        ASSERT_EQ(m[w].neen[1][1], 2);
        ASSERT_EQ(m[w].neen[1][2], 0);
        ASSERT_EQ(m[w].neen[2].size(), Mt::maxDegree);
        ASSERT_EQ(m[w].neen[2][0], 0);
        ASSERT_EQ(m[w].neen[2][1], 1);
        ASSERT_EQ(m[w].neen[2][2], 0);
    }
}

}  // namespace core_transformer_test
