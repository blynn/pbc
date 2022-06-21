#include "pbcxx.h"

#include <memory>
#include <string>
#include <fstream>
#include <gtest/gtest.h>
#include <utility>
#include <pbc.h>

using std::string;
using std::unique_ptr;
using std::make_unique;

// Just need some pairing params for testing
#define DEFAULT_PAIRING_PARAM "type a\n" \
"q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n" \
"h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776\n" \
"r 730750818665451621361119245571504901405976559617\n" \
"exp2 159\n" \
"exp1 107\n" \
"sign1 1\n" \
"sign0 1\n" \


namespace pbc {
namespace testing {

class ZrTest : public ::testing::Test {
 protected:
  ZrTest() : pairing_(DEFAULT_PAIRING_PARAM), test_(pairing_) {}

  const Pairing pairing_;
  ZrElement test_;
};

TEST_F(ZrTest, Construct) {
  EXPECT_TRUE(test_.is0());
  EXPECT_FALSE(test_.is1());
}

TEST_F(ZrTest, Set0) {
  test_.set_si(42);
  EXPECT_FALSE(test_.is0());
  test_.set0();
  EXPECT_TRUE(test_.is0());
}

TEST_F(ZrTest, Set1) {
  test_.set1();
  EXPECT_FALSE(test_.is0());
  EXPECT_TRUE(test_.is1());
}

TEST_F(ZrTest, SetRandom) {
  pbc_random_set_deterministic(12345678);
  test_.set_random();
  EXPECT_FALSE(test_.is0());
  EXPECT_FALSE(test_.is1());
}

TEST_F(ZrTest, FromHash) {
  test_.from_hash("AA5555AA");
  EXPECT_FALSE(test_.is0());
  EXPECT_FALSE(test_.is1());
}

TEST_F(ZrTest, Compare) {
  ZrElement test1(pairing_), test2(pairing_);
  test1.set_si(123);
  test2.set_si(123);
  EXPECT_TRUE(test1 == test2);
  EXPECT_FALSE(test1 != test2);
  EXPECT_TRUE(test1.compare(test2) == 0);
  test2.set_si(122);
  EXPECT_FALSE(test1 == test2);
  EXPECT_TRUE(test1 != test2);
  EXPECT_TRUE(test1.compare(test2) != 0);
  test1.set_si(-1);
  EXPECT_FALSE(test1 == test2);
  EXPECT_TRUE(test1 != test2);
  EXPECT_TRUE(test1.compare(test2) != 0);
}

TEST_F(ZrTest, CopyConstruct) {
  test_.set_si(555);
  ZrElement test2(test_);
  EXPECT_TRUE(test_ == test2);
}

TEST_F(ZrTest, Assignment) {
  ZrElement test1(pairing_), test2(pairing_);
  test1.set_si(76484);
  test2 = test1;
  EXPECT_TRUE(test1 == test2);
}

TEST_F(ZrTest, PlusEquals) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  test1.set_si(123);
  test2.set_si(456);
  test1 += test2;
  EXPECT_TRUE(test2 == exp.set_si(456));
  EXPECT_TRUE(test1 == exp.set_si(579));
}

TEST_F(ZrTest, Addition) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  EXPECT_TRUE(test1.set_si(123) + test2.set_si(456) == exp.set_si(579));
  ZrElement r = test1.set_si(456) + test2.set_si(-123);
  EXPECT_TRUE(r == exp.set_si(333));
}

TEST_F(ZrTest, MinusEquals) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  test1.set_si(555);
  test2.set_si(222);
  test1 -= test2;
  EXPECT_TRUE(test2 == exp.set_si(222));
  EXPECT_TRUE(test1 == exp.set_si(333));
}

TEST_F(ZrTest, Subtraction) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  EXPECT_TRUE(test1.set_si(1111) - test2.set_si(111) == exp.set_si(1000));
  ZrElement r = test1 - test2;
  EXPECT_TRUE(r == exp);
}

TEST_F(ZrTest, StarEquals) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  test1.set_si(364565);
  test2.set_si(646456);
  test1 *= test2;
  EXPECT_TRUE(test2 == exp.set_si(646456));
  EXPECT_TRUE(test1 == exp.set_si(364565L * 646456L));
}

TEST_F(ZrTest, Multiplication) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  EXPECT_TRUE(test1.set_si(476878368LL) * test2.set_si(12345787LL)
    == exp.set_si(476878368LL * 12345787LL));
  ZrElement r = test1 * test2;
  EXPECT_TRUE(r == exp);
}

TEST_F(ZrTest, DivEquals) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  test1.set_si(6873645351LL * 23);
  test2.set_si(23);
  test1 /= test2;
  EXPECT_TRUE(test2 == exp.set_si(23));
  EXPECT_TRUE(test1 == exp.set_si(6873645351LL));
}

TEST_F(ZrTest, Division) {
  ZrElement test1(pairing_), test2(pairing_), exp(pairing_);
  EXPECT_TRUE(test1.set_si(476878368LL * 17) / test2.set_si(17)
    == exp.set_si(476878368LL));
  ZrElement r = test1 / test2;
  EXPECT_TRUE(r == exp);
}

TEST_F(ZrTest, Twice) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  EXPECT_TRUE(test.twice() == exp.set_si(15686222LL));
  EXPECT_TRUE(test == exp.set_si(7843111LL));
}

TEST_F(ZrTest, TwiceInplace) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  EXPECT_TRUE(test.twice_inplace() == exp.set_si(15686222LL));
  EXPECT_TRUE(test == exp.set_si(15686222LL));
}

TEST_F(ZrTest, Square) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  EXPECT_TRUE(test.square() == exp.set_si(7843111LL * 7843111LL));
  EXPECT_TRUE(test == exp.set_si(7843111LL));
}

TEST_F(ZrTest, SquareInPlace) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  EXPECT_TRUE(test.square_inplace() == exp.set_si(7843111LL * 7843111LL));
  EXPECT_TRUE(test == exp.set_si(7843111LL * 7843111LL));
}

TEST_F(ZrTest, Neg) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  EXPECT_TRUE(test.neg() == exp.set_si(-7843111LL));
  EXPECT_TRUE(test == exp.set_si(7843111LL));
}

TEST_F(ZrTest, NegInplace) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  EXPECT_TRUE(test.neg_inplace() == exp.set_si(-7843111LL));
  EXPECT_TRUE(test == exp.set_si(-7843111LL));
}

TEST_F(ZrTest, Invert) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  exp.set1();
  exp /= test;
  EXPECT_TRUE(test.invert() == exp);
  EXPECT_TRUE(test == exp.set_si(7843111LL));
}

TEST_F(ZrTest, InvertInplace) {
  ZrElement test(pairing_), exp(pairing_);
  test.set_si(7843111LL);
  exp.set1();
  exp /= test;
  EXPECT_TRUE(test.invert_inplace() == exp);
  EXPECT_TRUE(test == exp);
}

TEST_F(ZrTest, ToFromBytes) {
  test_.from_hash("ABCD1234");
  string buf = test_.to_bytes();
  ZrElement test2(pairing_);
  EXPECT_TRUE(test_ == test2.from_bytes(buf));
  EXPECT_TRUE(test_ == test2);
}

TEST_F(ZrTest, ToFromBytesTooShort) {
  test_.from_hash("ABCD1234");
  string buf = test_.to_bytes();
  buf.resize(buf.size() - 1);
  ZrElement test2(pairing_);
  EXPECT_THROW(test2.from_bytes(buf), std::invalid_argument);
}

TEST_F(ZrTest, ToFromBytesTooLong) {
  test_.from_hash("ABCD1234");
  string buf = test_.to_bytes();
  buf += "blah";
  ZrElement test2(pairing_);
  EXPECT_THROW(test2.from_bytes(buf), std::invalid_argument);
}



template <class T>
class GElementTest : public ::testing::Test {
 protected:
  GElementTest() : pairing_(DEFAULT_PAIRING_PARAM), test_(pairing_) {}

  const Pairing pairing_;
  T test_;
};

typedef ::testing::Types<G1Element, G2Element, GTElement> Implementations;
TYPED_TEST_SUITE(GElementTest, Implementations);

TYPED_TEST(GElementTest, Construct) {
  EXPECT_TRUE(this->test_.is0());
  EXPECT_TRUE(this->test_.is1());
}

TYPED_TEST(GElementTest, FromHash) {
  this->test_.from_hash("AA5555AA");
  EXPECT_FALSE(this->test_.is0());
  EXPECT_FALSE(this->test_.is1());
}

TYPED_TEST(GElementTest, Set0) {
  this->test_.from_hash("AA5555AA");
  EXPECT_FALSE(this->test_.is0());
  EXPECT_FALSE(this->test_.is1());
  this->test_.set0();
  EXPECT_TRUE(this->test_.is0());
}

TYPED_TEST(GElementTest, Set1) {
  this->test_.from_hash("AA5555AA");
  EXPECT_FALSE(this->test_.is0());
  EXPECT_FALSE(this->test_.is1());
  this->test_.set1();
  EXPECT_TRUE(this->test_.is1());
}

TYPED_TEST(GElementTest, SetRandom) {
  pbc_random_set_deterministic(12345678);
  this->test_.set_random();
  EXPECT_FALSE(this->test_.is0());
  EXPECT_FALSE(this->test_.is1());
}

TYPED_TEST(GElementTest, Compare) {
  TypeParam test1(this->pairing_), test2(this->pairing_);
  test1.from_hash("AA55");
  test2.from_hash("AA55");
  EXPECT_TRUE(test1 == test2);
  EXPECT_FALSE(test1 != test2);
  EXPECT_TRUE(test1.compare(test2) == 0);
  test2.from_hash("4321");
  EXPECT_FALSE(test1 == test2);
  EXPECT_TRUE(test1 != test2);
  EXPECT_TRUE(test1.compare(test2) != 0);
}

TYPED_TEST(GElementTest, CopyConstruct) {
  TypeParam test1(this->pairing_);
  test1.from_hash("AA5555AA");
  TypeParam test2(test1);
  EXPECT_TRUE(test1 == test2);
}

TYPED_TEST(GElementTest, Assignment) {
  TypeParam test1(this->pairing_), test2(this->pairing_);
  test1.from_hash("ABCDEF01");
  test2 = test1;
  EXPECT_TRUE(test1 == test2);
}

TYPED_TEST(GElementTest, PowZn) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  a.set_si(2);
  EXPECT_TRUE(g.pow_zn(a) == g.square());
  a.set_si(4);
  TypeParam g1 = g.pow_zn(a);
  EXPECT_TRUE(g1 == g.square().square());
}

TYPED_TEST(GElementTest, PowZnInPlace) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam exp = g.square().square();
  a.set_si(4);
  g.pow_zn_inplace(a);
  EXPECT_TRUE(g == exp);
}

TYPED_TEST(GElementTest, MultiplyZn) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  exp = g.twice();
  a.set_si(2);
  EXPECT_TRUE(g * a == exp);
  exp.twice_inplace();
  a.set_si(4);
  TypeParam g1 = g * a;
  EXPECT_TRUE(g1 == exp);
  EXPECT_TRUE(g*a == a*g);
}

TYPED_TEST(GElementTest, StarEqualsZn) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  exp = g;
  exp.twice_inplace().twice_inplace();
  a.set_si(4);
  g *= a;
  EXPECT_TRUE(g == exp);
}

TYPED_TEST(GElementTest, PowSameAsMulZn) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  a.set_si(5);
  EXPECT_TRUE(g.pow_zn(a) == g * a);
}

TYPED_TEST(GElementTest, PlusEquals) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g * a.set_si(2);
  TypeParam g2 = g * a.set_si(3);
  exp = g2;
  g1 += g2;
  EXPECT_TRUE(g2 == exp);
  EXPECT_TRUE(g1 == g * a.set_si(5));
}

TYPED_TEST(GElementTest, Addition) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g * a.set_si(2);
  TypeParam g2 = g * a.set_si(3);
  TypeParam gr = g1 + g2;
  EXPECT_TRUE(gr == g * a.set_si(5));
  EXPECT_TRUE(g1 == g * a.set_si(2));
  EXPECT_TRUE(g2 == g * a.set_si(3));
  EXPECT_TRUE(g1 + g2 == g * a.set_si(5));
}

TYPED_TEST(GElementTest, MinusEquals) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g * a.set_si(5);
  TypeParam g2 = g * a.set_si(3);
  exp = g2;
  g1 -= g2;
  EXPECT_TRUE(g2 == exp);
  EXPECT_TRUE(g1 == g * a.set_si(2));
}

TYPED_TEST(GElementTest, Subtraction) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g * a.set_si(5);
  TypeParam g2 = g * a.set_si(3);
  TypeParam gr = g1 - g2;
  EXPECT_TRUE(gr == g * a.set_si(2));
  EXPECT_TRUE(g1 == g * a.set_si(5));
  EXPECT_TRUE(g2 == g * a.set_si(3));
  EXPECT_TRUE(g1 - g2 == g * a.set_si(2));
}

TYPED_TEST(GElementTest, StarEquals) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g.pow_zn(a.set_si(2));
  TypeParam g2 = g.pow_zn(a.set_si(3));
  exp = g2;
  g1 *= g2;
  EXPECT_TRUE(g2 == exp);
  EXPECT_TRUE(g1 == g.pow_zn(a.set_si(5)));
}

TYPED_TEST(GElementTest, Multiplication) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g.pow_zn(a.set_si(2));
  TypeParam g2 = g.pow_zn(a.set_si(3));
  TypeParam gr = g1 * g2;
  EXPECT_TRUE(gr == g.pow_zn(a.set_si(5)));
  EXPECT_TRUE(g1 == g.pow_zn(a.set_si(2)));
  EXPECT_TRUE(g2 == g.pow_zn(a.set_si(3)));
  EXPECT_TRUE(g1 * g2 == g.pow_zn(a.set_si(5)));
}

TYPED_TEST(GElementTest, MultiplicationSameAsAddition) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g.pow_zn(a.set_si(7));
  TypeParam g2 = g.pow_zn(a.set_si(5));
  EXPECT_TRUE(g1 + g2 == g1 * g2);
}

TYPED_TEST(GElementTest, DivEquals) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g.pow_zn(a.set_si(5));
  TypeParam g2 = g.pow_zn(a.set_si(3));
  exp = g2;
  g1 /= g2;
  EXPECT_TRUE(g2 == exp);
  EXPECT_TRUE(g1 == g.pow_zn(a.set_si(2)));
}

TYPED_TEST(GElementTest, Division) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g.pow_zn(a.set_si(5));
  TypeParam g2 = g.pow_zn(a.set_si(3));
  TypeParam gr = g1 / g2;
  EXPECT_TRUE(gr == g.pow_zn(a.set_si(2)));
  EXPECT_TRUE(g1 == g.pow_zn(a.set_si(5)));
  EXPECT_TRUE(g2 == g.pow_zn(a.set_si(3)));
  EXPECT_TRUE(g1 / g2 == g.pow_zn(a.set_si(2)));
}

TYPED_TEST(GElementTest, DivisionSameAsSubtraction) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g.pow_zn(a.set_si(7));
  TypeParam g2 = g.pow_zn(a.set_si(5));
  EXPECT_TRUE(g1 - g2 == g1 / g2);
}

TYPED_TEST(GElementTest, Twice) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g * a.set_si(5);
  EXPECT_TRUE(g1.twice() == g * a.set_si(10));
  EXPECT_TRUE(g1 == g * a.set_si(5));
}

TYPED_TEST(GElementTest, TwiceInplace) {
  TypeParam g(this->pairing_);
  ZrElement a(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g * a.set_si(5);
  EXPECT_TRUE(g1.twice_inplace() == g * a.set_si(10));
  EXPECT_TRUE(g1 == g * a.set_si(10));
}

TYPED_TEST(GElementTest, Square) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g;
  EXPECT_TRUE(g1.square() == g * g);
  EXPECT_TRUE(g1 == g);
}

TYPED_TEST(GElementTest, SquareInplace) {
  TypeParam g(this->pairing_), exp(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g;
  EXPECT_TRUE(g1.square_inplace() == g * g);
  EXPECT_TRUE(g1 == g * g);
}

TYPED_TEST(GElementTest, Neg) {
  TypeParam g(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g;
  EXPECT_TRUE((g1.neg() + g).is0());
  EXPECT_TRUE(g1 == g);
}

TYPED_TEST(GElementTest, NegInplace) {
  TypeParam g(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g;
  EXPECT_TRUE((g1.neg_inplace() + g).is0());
  EXPECT_TRUE((g1 + g).is0());
}

TYPED_TEST(GElementTest, Invert) {
  TypeParam g(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g;
  EXPECT_TRUE((g1.invert() * g).is1());
  EXPECT_TRUE(g1 == g);
}

TYPED_TEST(GElementTest, InvertInplace) {
  TypeParam g(this->pairing_);
  g.from_hash("ABCDEF01");
  TypeParam g1 = g;
  EXPECT_TRUE((g1.invert_inplace() * g).is1());
  EXPECT_TRUE((g1 * g).is1());
}

TYPED_TEST(GElementTest, ToFromBytes) {
  this->test_.from_hash("ABCD1234");
  string buf = this->test_.to_bytes();
  TypeParam test2(this->pairing_);
  EXPECT_TRUE(this->test_ == test2.from_bytes(buf));
  EXPECT_TRUE(this->test_ == test2);
}

TYPED_TEST(GElementTest, ToFromBytesTooShort) {
  this->test_.from_hash("ABCD1234");
  string buf = this->test_.to_bytes();
  buf.resize(buf.size() - 1);
  ZrElement test2(this->pairing_);
  EXPECT_THROW(test2.from_bytes(buf), std::invalid_argument);
}

TYPED_TEST(GElementTest, ToFromBytesTooLong) {
  this->test_.from_hash("ABCD1234");
  string buf = this->test_.to_bytes();
  buf += "blah";
  ZrElement test2(this->pairing_);
  EXPECT_THROW(test2.from_bytes(buf), std::invalid_argument);
}



TEST(Pairing, IsSymmetric) {
  Pairing pairing(DEFAULT_PAIRING_PARAM);
  EXPECT_TRUE(pairing.is_symmetric());
}

TEST(Pairing, PairingConsistency) {
  Pairing pairing(DEFAULT_PAIRING_PARAM);
  G1Element g(pairing);
  g.from_hash("FEDCBA9876");
  G2Element h(pairing);  
  h.from_hash("0123456789");
  ZrElement a(pairing);
  GTElement gt(pairing);
  a.set_si(42);
  gt = pairing(g * a, h);
  EXPECT_TRUE(gt == pairing(g, h * a));
  EXPECT_TRUE(gt == pairing(g, h).pow_zn(a));
  gt.set1();
  ZrElement a2(pairing);
  a.set_si(17);
  a2.set_si(101);
  pairing.apply(gt, g * a, h * a2);
  EXPECT_TRUE(gt == pairing(g * a2, h * a));
  EXPECT_TRUE(gt == pairing(g, h).pow_zn(a * a2));
}

TEST(Pairing, PairingConsistencySymmetricG1) {
  Pairing pairing(DEFAULT_PAIRING_PARAM);
  G1Element g(pairing), h(pairing);
  g.from_hash("FEDCBA9876");
  h.from_hash("0123456789");
  ZrElement a(pairing);
  GTElement gt(pairing);
  a.set_si(42);
  gt = pairing(g * a, h);
  EXPECT_TRUE(gt == pairing(g, h * a));
  EXPECT_TRUE(gt == pairing(g, h).pow_zn(a));
  gt.set1();
  ZrElement a2(pairing);
  a.set_si(17);
  a2.set_si(101);
  pairing.apply(gt, g * a, h * a2);
  EXPECT_TRUE(gt == pairing(g * a2, h * a));
  EXPECT_TRUE(gt == pairing(g, h).pow_zn(a * a2));
}

TEST(Pairing, PairingConsistencySymmetricG2) {
  Pairing pairing(DEFAULT_PAIRING_PARAM);
  G2Element g(pairing), h(pairing);
  g.from_hash("FEDCBA9876");
  h.from_hash("0123456789");
  ZrElement a(pairing);
  GTElement gt(pairing);
  a.set_si(42);
  gt = pairing(g * a, h);
  EXPECT_TRUE(gt == pairing(g, h * a));
  EXPECT_TRUE(gt == pairing(g, h).pow_zn(a));
  gt.set1();
  ZrElement a2(pairing);
  a.set_si(17);
  a2.set_si(101);
  pairing.apply(gt, g * a, h * a2);
  EXPECT_TRUE(gt == pairing(g * a2, h * a));
  EXPECT_TRUE(gt == pairing(g, h).pow_zn(a * a2));
}

} // namespace testing
} // namespace pbc


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
