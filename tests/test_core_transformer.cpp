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

    using Msgr = mitosim::Msgr;
    using Mt = mitosim::Segment<3>;

    using CoreTransformer::clnum;
    using CoreTransformer::copy_neigs;
    using CoreTransformer::fuse_antiparallel;
    using CoreTransformer::fuse_parallel;
    using CoreTransformer::fuse_to_loop;
    using CoreTransformer::mt;
    using CoreTransformer::mtnum;
    using CoreTransformer::rename_mito;
    using CoreTransformer::update_neigs;
    using mitosim::Structure<Mt>::add_disconnected_segment;

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
    mitosim::CoreTransformer<Mt> ct {msgr};
    ct.add_disconnected_segment(3);
//    ct.add_disconnected_segment(6);
}

}  // namespace core_transformer_test
