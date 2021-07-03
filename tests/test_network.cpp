#include <filesystem>
#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"
#include "../segment.h"
#include "../network.h"

namespace network_test {

class NetworkTest
    : public testing::Test {

protected:

    using Mt = mitosim::Segment<3>;
    using Network = mitosim::Network<Mt>;
    using real = mitosim::real;
    using szt = mitosim::szt;

    const std::string workingDir {std::filesystem::current_path()
                                  / "tests" / "data/"};
    const std::string fnameSuffix {"sample"};
    const std::string runName {"42"};

    NetworkTest()
        : msgr {}
        , conf {workingDir, fnameSuffix, runName, msgr}
        , rnd {std::make_unique<mitosim::RandFactory>(10, msgr)}
    {}

    mitosim::Msgr msgr;
    mitosim::Config<real> conf;
    std::unique_ptr<mitosim::RandFactory> rnd;
};

TEST_F(NetworkTest, Constructor)
{
    Network ntw {conf, *rnd, msgr};

    ASSERT_TRUE(ntw.clagl.empty());
    ASSERT_TRUE(ntw.glm.empty());
    ASSERT_TRUE(ntw.gla.empty());
    ASSERT_TRUE(ntw.mt.empty());
    for (const auto& n : ntw.nn)
        EXPECT_EQ(n, 0);
    EXPECT_EQ(ntw.mtnum, 0);
    EXPECT_EQ(ntw.clnum, 0);
    EXPECT_EQ(ntw.mtmass, 0);
    ASSERT_TRUE(ntw.clmt.empty());
    ASSERT_TRUE(ntw.mt11.empty());
    ASSERT_TRUE(ntw.mtc11.empty());
    ASSERT_TRUE(ntw.mt22.empty());
    ASSERT_TRUE(ntw.mtc22.empty());
    ASSERT_TRUE(ntw.mt33.empty());
    ASSERT_TRUE(ntw.mtc33.empty());
    ASSERT_TRUE(ntw.mt13.empty());
    ASSERT_TRUE(ntw.mtc13.empty());
}

}  // namespace network_test
