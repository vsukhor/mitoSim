#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"
#include "../segment.h"
#include "../network.h"
#include "ntw_fusion11.h"

namespace ntw_fusion11_test {

class NtwFusion11Test
    : public testing::Test {

protected:

    using Mt = mitosim::Segment<3>;
    using Network = mitosim::Network<Mt>;
    using RandFactory = mitosim::RandFactory;
    using real = mitosim::real;
    using szt = mitosim::szt;
    using NtwFusion11 = mitosim::NtwFusion11<Network>;

    static constexpr auto minLL =
        mitosim::Structure<typename Network::ST>::minLoopLength;

    const std::string workingDir {std::filesystem::current_path()
                                  / "tests" / "data/"};
    const std::string fnameSuffix {"sample"};
    const std::string runName {"42"};

    NtwFusion11Test()
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

TEST_F(NtwFusion11Test, Constructor)
{
    NtwFusion11 nf {ntw};

    EXPECT_EQ(nf.get_cnd().size(), 0);
}

TEST_F(NtwFusion11Test, SetProp)
{
    constexpr std::array<szt,4> len {4, 9, 5, 6};

    for (const auto u : len)
        ntw.add_disconnected_segment(u);
    ntw.fuse12(1, 1, 2, 3);
    // Now there are 3 13-segments and 2 11-segments in the network.

    NtwFusion11 nf {ntw};
    ntw.populate_cluster_vectors();
    nf.set_prop();

    EXPECT_EQ(nf.get_cnd().size(), 7*6/2);
}

}  // namespace ntw_fusion11_test
