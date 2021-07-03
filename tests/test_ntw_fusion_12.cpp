#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"
#include "../segment.h"
#include "../network.h"
#include "ntw_fusion12.h"

namespace ntw_fusion12_test {

class NtwFusion12Test
    : public testing::Test {

protected:

    using Mt = mitosim::Segment<3>;
    using Network = mitosim::Network<Mt>;
    using RandFactory = mitosim::RandFactory;
    using real = mitosim::real;
    using szt = mitosim::szt;
    using NtwFusion12 = mitosim::NtwFusion12<Network>;

    static constexpr auto minLL =
        mitosim::Structure<typename Network::ST>::minLoopLength;

    const std::string workingDir {std::filesystem::current_path()
                                  / "tests" / "data/"};
    const std::string fnameSuffix {"sample"};
    const std::string runName {"42"};

    NtwFusion12Test()
        : msgr {}
        , conf {workingDir, fnameSuffix, runName, msgr}
        , rnd {std::make_unique<mitosim::RandFactory>(10, msgr)}
        , ntw {conf, *rnd, msgr}
    {}

    mitosim::Msgr msgr;
    mitosim::Config<real> conf;
    std::unique_ptr<RandFactory> rnd;
    Network ntw;
};

TEST_F(NtwFusion12Test, Constructor)
{
    NtwFusion12 nf {ntw};

    EXPECT_EQ(nf.get_cnd().size(), 0);
}

TEST_F(NtwFusion12Test, SetProp1)
{
    // Tests propensity contribution from 11-segments.
    constexpr std::array<szt,2> len {4, 3};
    const auto lensum = std::accumulate(len.begin(), len.end(), 0);

    for (const auto u : len)
        ntw.add_disconnected_segment(u);

    NtwFusion12 nf {ntw};
    ntw.populate_cluster_vectors();
    nf.set_prop();

    EXPECT_EQ(nf.get_cnd().size(), // 16
              Mt::numEnds * ((lensum - 2) + (lensum - 2*minLL)));
}

TEST_F(NtwFusion12Test, SetProp2)
{
    // Tests propensity in a system of 1x11 + 3*13 segments.
    constexpr std::array<szt,3> len {4, 4, 4};

    for (const auto u : len)
        ntw.add_disconnected_segment(u);
    ntw.fuse12(1, 1, 2, 2);

    NtwFusion12 nf {ntw};
    ntw.populate_cluster_vectors();
    nf.set_prop();

    EXPECT_EQ(nf.get_cnd().size(), 35);
}


}  // namespace ntw_fusion12_test
