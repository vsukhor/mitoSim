#include <filesystem>

#include "gtest/gtest.h"

#include "../config.h"
#include "../definitions.h"

namespace config_test {

class ConfigTest
    : public testing::Test {

protected:
    using real = mitosim::real;
    using Config = mitosim::Config<real>;

    const std::string workingDir {std::filesystem::current_path() / "tests" / "data/"};
    const std::string fnameSuffix {"sample"};
    const std::string runName {"42"};

    mitosim::Msgr msgr;
    const Config cfg;

    ConfigTest()
        : msgr {}
        , cfg {workingDir, fnameSuffix, runName, msgr}
    {}

};

TEST_F(ConfigTest, Constructor)
{
    ASSERT_STREQ(cfg.workingDirOut.c_str(), workingDir.c_str());
    ASSERT_STREQ(cfg.fnameSuffix.c_str(), fnameSuffix.c_str());
    ASSERT_STREQ(cfg.runName.c_str(), runName.c_str());
    ASSERT_FLOAT_EQ(cfg.timeTotal, 20.f);
    ASSERT_EQ(cfg.logFrequency, 100);
    ASSERT_EQ(cfg.saveFrequency, 20000);
    ASSERT_FLOAT_EQ(cfg.edgeLength, 0.2f);
    ASSERT_EQ(cfg.mtmassini, 1000);
    ASSERT_EQ(cfg.segmassini, 20);
    ASSERT_TRUE(cfg.use_fission);
    ASSERT_FLOAT_EQ(cfg.rate_fission, 1.f);
    ASSERT_TRUE(cfg.use_11_fusion);
    ASSERT_FLOAT_EQ(cfg.fusion_rate_11, .1f);
    ASSERT_TRUE(cfg.use_12_fusion);
    ASSERT_FLOAT_EQ(cfg.fusion_rate_12, .002f);
    ASSERT_TRUE(cfg.use_1L_fusion);
    ASSERT_FLOAT_EQ(cfg.fusion_rate_1L, .005f);
}

}
