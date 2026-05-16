#include <gtest/gtest.h>
#include <mphys/mphys.hpp>

TEST(Version, Major)
{
    EXPECT_EQ(mphys::version_major, 1);
}
