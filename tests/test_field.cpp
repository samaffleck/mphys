#include <gtest/gtest.h>

#include "mphys/field.hpp"

namespace {

using mphys::Field;

TEST(Field, DefaultConstructorIsEmpty) {
  Field f;
  EXPECT_TRUE(f.name.empty());
  EXPECT_EQ(f.NCells(), 0);
}

TEST(Field, SizedConstructorFillsInitValue) {
  Field f("c", 3, 1.5);
  EXPECT_EQ(f.name, "c");
  EXPECT_EQ(f.NCells(), 3);
  for (int i = 0; i < 3; ++i) EXPECT_DOUBLE_EQ(f[i], 1.5);
}

TEST(Field, DefaultInitValueIsZero) {
  Field f("c", 2);
  EXPECT_DOUBLE_EQ(f[0], 0.0);
  EXPECT_DOUBLE_EQ(f[1], 0.0);
}

TEST(Field, SubscriptReadWrite) {
  Field f("c", 2);
  f[0] = 7.0;
  f[1] = -3.0;
  EXPECT_DOUBLE_EQ(f[0], 7.0);
  EXPECT_DOUBLE_EQ(f[1], -3.0);

  const Field& cf = f;
  EXPECT_DOUBLE_EQ(cf[0], 7.0);
}

TEST(Field, Addition) {
  Field a("a", 3);
  Field b("b", 3);
  a.values = {1, 2, 3};
  b.values = {10, 20, 30};
  Field c = a + b;
  EXPECT_DOUBLE_EQ(c[0], 11);
  EXPECT_DOUBLE_EQ(c[1], 22);
  EXPECT_DOUBLE_EQ(c[2], 33);
}

TEST(Field, Subtraction) {
  Field a("a", 2);
  Field b("b", 2);
  a.values = {5, 8};
  b.values = {1, 10};
  Field c = a - b;
  EXPECT_DOUBLE_EQ(c[0], 4);
  EXPECT_DOUBLE_EQ(c[1], -2);
}

TEST(Field, ScalarMultiplyBothOrders) {
  Field a("a", 2);
  a.values = {2, -3};
  Field left = a * 4.0;
  Field right = 4.0 * a;
  EXPECT_DOUBLE_EQ(left[0], 8);
  EXPECT_DOUBLE_EQ(left[1], -12);
  EXPECT_DOUBLE_EQ(right[0], 8);
  EXPECT_DOUBLE_EQ(right[1], -12);
}

TEST(Field, UnaryNegation) {
  Field a("a", 2);
  a.values = {2, -3};
  Field n = -a;
  EXPECT_DOUBLE_EQ(n[0], -2);
  EXPECT_DOUBLE_EQ(n[1], 3);
}

TEST(Field, CompoundAssignment) {
  Field a("a", 2);
  Field b("b", 2);
  a.values = {1, 1};
  b.values = {2, 3};
  a += b;
  EXPECT_DOUBLE_EQ(a[0], 3);
  EXPECT_DOUBLE_EQ(a[1], 4);
  a -= b;
  EXPECT_DOUBLE_EQ(a[0], 1);
  EXPECT_DOUBLE_EQ(a[1], 1);
}

}  // namespace
