#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"
#include "../segment.h"
#include "../network.h"
#include "ntw_fusion1u.h"

namespace ntw_fusion1u_test {

class NtwFusion1UTest
    : public testing::Test {

protected:

    using Mt = mitosim::Segment<3>;
    using Network = mitosim::Network<Mt>;
    using RandFactory = mitosim::RandFactory;
    using Msgr = mitosim::Msgr;
    using Config = mitosim::Config<mitosim::real>;
    using real = mitosim::real;
    using szt = mitosim::szt;
    using NtwFusion1U = mitosim::NtwFusion1U<Network>;

    class Ntw : public Network {

    public:

        using Network::fuse_to_loop;

        explicit Ntw(
            const Config& cfg,
            RandFactory& rnd,
            Msgr& msgr
        ) : Network {cfg, rnd, msgr}
        {}
    };

    static constexpr auto minLL =
        mitosim::Structure<typename Network::ST>::minLoopLength;

    const std::string workingDir {std::filesystem::current_path()
                                  / "tests" / "data/"};
    const std::string fnameSuffix {"sample"};
    const std::string runName {"42"};

    NtwFusion1UTest()
        : msgr {}
        , conf {workingDir, fnameSuffix, runName, msgr}
        , rnd {std::make_unique<mitosim::RandFactory>(10, msgr)}
        , ntw {conf, *rnd, msgr}
    {}

    Msgr msgr;
    Config conf;
    std::unique_ptr<RandFactory> rnd;
    Ntw ntw;
};

TEST_F(NtwFusion1UTest, Constructor)
{
    NtwFusion1U nf {ntw};

    EXPECT_EQ(nf.get_cnd().size(), 0);
}

TEST_F(NtwFusion1UTest, SetProp)
{
    // Tests propensity contribution from 11-segments to separate cycles.
    constexpr std::array<szt,4> len {4, 8, 3, 5};

    for (const auto u : len)
        ntw.add_disconnected_segment(u);
    ntw.fuse_to_loop(1);
    ntw.fuse_to_loop(2);

    NtwFusion1U nf {ntw};
    ntw.populate_cluster_vectors();
    nf.set_prop();

    EXPECT_EQ(nf.get_cnd().size(), 8);
}


}  // namespace ntw_fusion1u_test
