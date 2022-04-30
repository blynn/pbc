#ifndef PBC_PBCXX_H
#define PBC_PBCXX_H

#include <iostream>
#include <string>
#include <pbc/pbc.h>

// WARNING WARNING WARNING!!!
// Elements must be destroyed before pairings, otherwise this will SEGFAULT!
// TODO: find a way to make sure elements are destroyed before pairings.

namespace pbc {

using std::string;

// Forward declaration of Pairing class so it can be used to initialize elements
class Pairing;

// Use an anonymous namespace to prevent internals from leaking
namespace {

// Function pointer used for initializing elements. Used as a template parameter
// Should not be used anywhere outside of this file.
typedef void (*__element_initializer_t)(element_t, pairing_t);

/**
 * Generic wrapper for PBC elements. First template parameter is a function
 * pointer for how to initialize an element, serves as differentiation element
 * types. Second template parameter is the type of the derived class. This is so
 * that the derived class can be returned from operators, instead of the base
 * class. Doing so provides a cleaner interface for end users, we just have to
 * be very careful in this class with how we use it.
 *
 * A note on const. PBC basically has no notion of const. If you define a const
 * element_t, you can't use it anywhere. You'll notice that the element being
 * wrapped is marked mutable. This is so that even when the instance of the
 * ElementWrapper is const, we can still call PBC functions. We try here to
 * never modify the element in a const method, but that's best effort and relies
 * on PBC being well behaved.
 *
 * NB: You probably don't want to instantiate this class. See below for the
 * classes end users would actually use.
 */
template <__element_initializer_t init, class Derived>
class ElementWrapper {
  // GElementBase needs to be a friend so that pow_zn and mul_zn can work
  // access elements across templates.
  template<__element_initializer_t, class>
  friend class GElementBase;

  // Allow Pairing to access elements for computing pairings.
  friend class ::pbc::Pairing;

  // Shortcut for the type of the current template instatiation.
  typedef ElementWrapper<init, Derived> MyType;

 protected:
  // Only used in very special cases, inits element same as other but doesn't
  // copy the value. Used in some operators.
  ElementWrapper(element_t e) {
    element_init_same_as(v_, e);
  }

  // PBC doesn't use const. So we need this to be mutable for cases when we
  // don't actually modify it, but the compiler doesn't know that.
  mutable element_t v_;

 public:
  // Don't allow default construction, need to ensure element is initialized.
  ElementWrapper() = delete;

  // Initialize an element given a pairing.
  ElementWrapper(const Pairing& pairing);

  // Initializes and copies an element.
  ElementWrapper(const MyType& e) {
    element_init_same_as(v_, e.v_);
    element_set(v_, e.v_);
  }

  // Make sure to clean up!
  ~ElementWrapper() {
    element_clear(v_);
  }

  Derived& operator=(const MyType& other) {
    if (this != &other) {
      element_set(v_, other.v_);
    }
    return *static_cast<Derived*>(this);
  }

  int compare(const MyType& e) const { return element_cmp(v_, e.v_); }

  bool is0() const { return element_is0(v_); }

  bool is1() const { return element_is1(v_); }

  // Here's an example of returning the derived class and not the base class.
  // If we didn't have Derived, users would get back a reference to
  // ElementWrapper and not the class they actually created.
  Derived& set0() {
    element_set0(v_);
    return *static_cast<Derived*>(this);
  }

  Derived& set1() {
    element_set1(v_);
    return *static_cast<Derived*>(this);
  }

  Derived& set_random() {
    element_random(v_);
    return *static_cast<Derived*>(this);
  }

  Derived& from_hash(const string& data) {
    element_from_hash(v_, (void*)data.data(), data.size());
    return *static_cast<Derived*>(this);
  }

  // Utility function to print a prefix string followed by the element.
  void debug_print(const char* str) const {
    element_printf("%s %B\n", str, v_);
  }

  inline bool operator==(const MyType& other) const {
    return compare(other) == 0;
  }

  inline bool operator!=(const MyType& other) const {
    return compare(other) != 0;
  }

  Derived& operator+=(const MyType& rhs) {
    element_add(v_, v_, rhs.v_);
    return *static_cast<Derived*>(this);
  }

  friend Derived operator+(const MyType& lhs, const MyType& rhs) {
    Derived result(lhs.v_);
    element_add(result.v_, lhs.v_, rhs.v_);
    return result;
  }

  Derived& operator-=(const MyType& rhs) {
    element_sub(v_, v_, rhs.v_);
    return *static_cast<Derived*>(this);
  }

  friend Derived operator-(const MyType& lhs, const MyType& rhs) {
    Derived result(lhs.v_);
    element_sub(result.v_, lhs.v_, rhs.v_);
    return result;
  }

  Derived& operator*=(const MyType& rhs) {
    element_mul(v_, v_, rhs.v_);
    return *static_cast<Derived*>(this);
  }

  friend Derived operator*(const MyType& lhs, const MyType& rhs) {
    Derived result(lhs.v_);
    element_mul(result.v_, lhs.v_, rhs.v_);
    return result;
  }

  Derived& operator/=(const MyType& rhs) {
    element_div(v_, v_, rhs.v_);
    return *static_cast<Derived*>(this);
  }

  friend Derived operator/(const MyType& lhs, const MyType& rhs) {
    Derived result(lhs.v_);
    element_div(result.v_, lhs.v_, rhs.v_);
    return result;
  }

  // Double the element and return the new value.
  Derived twice() {
    Derived result(v_);
    element_double(result.v_, v_);
    return result;
  }

  // Double the element in place.
  Derived& twice_inplace() {
    element_double(v_, v_);
    return *static_cast<Derived*>(this);
  }

  // Halve the element in place.
  // Apparently, halve doesn't actually work in anything other than Zr.
  // Derived& halve() {
  //   element_halve(v_, v_);
  //   return *static_cast<Derived*>(this);
  // }

  // Square the element and return the new value.
  Derived square() {
    Derived result(v_);
    element_square(result.v_, v_);
    return result;
  }

  // Square the element in place.
  Derived& square_inplace() {
    element_square(v_, v_);
    return *static_cast<Derived*>(this);
  }

  // Negative of the element and return the new value.
  Derived neg() {
    Derived result(v_);
    element_neg(result.v_, v_);
    return result;
  }

  // Negative of the element in place.
  Derived& neg_inplace() {
    element_neg(v_, v_);
    return *static_cast<Derived*>(this);
  }

  // Invert the element and return the new value.
  Derived invert() {
    Derived result(v_);
    element_invert(result.v_, v_);
    return result;
  }

  // Invert the element in place.
  Derived& invert_inplace() {
    element_invert(v_, v_);
    return *static_cast<Derived*>(this);
  }

  // Get a byte representation of the current element.
  string to_bytes() const {
    string buf(element_length_in_bytes(v_), 0);
    element_to_bytes((unsigned char*)&buf[0], v_);
    return buf;
  }

  // Extract element data from byte string. Note that PBC method is very unsafe
  // as it could read past the end of the buffer!
  Derived& from_bytes(const string& buf) {
    if (element_length_in_bytes(v_) > buf.size()) {
      throw std::invalid_argument("Buffer too small for element");
    }
    auto n = element_from_bytes(v_, (unsigned char*)&buf[0]);
    if (n != buf.size()) {
      // TODO: not sure if this should actually throw. Might be ok?
      throw std::invalid_argument("wrong number of bytes read");
    }
    return *static_cast<Derived*>(this);
  }
};

} // namespace

/**
 * Wrapper for elements in the Zr integer field. Allows for setting from signed
 * integers and mpz data types. Can also be used a pow_zn exponents.
 */
class ZrElement : public ElementWrapper<element_init_Zr, ZrElement> {
 protected:
  // Needed so that operators in ElementWrapper can do partial constructions of
  // this derived class, without defining operators here.
  friend class ElementWrapper<element_init_Zr, ZrElement>;
  ZrElement(element_t e) : ElementWrapper(e) {}

 public:
  // Construct Zr element from definition in pairing.
  ZrElement(const Pairing& pairing) : ElementWrapper(pairing) {}

  // Copy another Zr element, both parameters and data.
  ZrElement(const ZrElement& e) : ElementWrapper(e) {}

  // Set Zr element to a signed integer value.
  ZrElement& set_si(signed long int i) {
    element_set_si(v_, i);
    return *this;
  }

  // Set Zr element to a gmp mpz value.
  ZrElement& set_mpz(mpz_t z) {
    element_set_mpz(v_, z);
    return *this;
  }
};

// Keep GElementBase anonymous to prevent leakage
namespace {

/**
 * Generic wrapper for G1, G2, GT group elements. Extensions for doing pow_zn
 * with a Zr element.
 */
template <__element_initializer_t init, class Derived>
class GElementBase : public ElementWrapper<init, Derived> {
 protected:
  typedef ElementWrapper<init, Derived> MyWrapper;
  // Needed so that operators in ElementWrapper can do partial constructions of
  // this derived class, without defining operators here.
  friend class ElementWrapper<init, Derived>;
  GElementBase(element_t e) : MyWrapper(e) {}

 public:
  // Construct group element from definition in pairing.
  GElementBase(const Pairing& pairing) : MyWrapper(pairing) {}

  // Copy another group element of the same type, including data.
  GElementBase(const GElementBase& e) : MyWrapper(e) {}

  // Raise this group element to the power of zn and return result.
  // Does not modify this instance.
  Derived pow_zn(const ZrElement& zn) const {
    Derived result(this->v_);
    element_pow_zn(result.v_, this->v_, zn.v_);
    return result;
  }

  // Raise this group element to the power of zn in place.
  Derived& pow_zn_inplace(const ZrElement& zn) {
    element_pow_zn(this->v_, this->v_, zn.v_);
    return *static_cast<Derived*>(this);
  }

  // Similar to pow_zn_inplace above, multiply by Zn in place.
  Derived& operator*=(const ZrElement& rhs) {
    element_mul_zn(this->v_, this->v_, rhs.v_);
    return *static_cast<Derived*>(this);
  }

  // Because we define a new operator*= have to pass through to base class.
  Derived& operator*=(const GElementBase& rhs) {
    return this->MyWrapper::operator*=(rhs);
  }

  // Similar to pow_zn above, multiply by a Zn element without modifying this.
  friend Derived operator*(const GElementBase& g, const ZrElement& z) {
    Derived result(g.v_);
    element_mul_zn(result.v_, g.v_, z.v_);
    return result;
  }

  // Same thing just with order swapped.
  friend Derived operator*(const ZrElement& z, const GElementBase& g) {
    return g * z;
  }

  // Again, pass through to base class.
  friend Derived operator*(const GElementBase& lhs,
      const GElementBase& rhs) {
    return static_cast<MyWrapper>(lhs) * static_cast<MyWrapper>(rhs);
  }
};

} // namespace

#define __DEFINE_ELEMENT_CLASS(name, init)                   \
class name : public GElementBase<init, name> {               \
 protected:                                                  \
  friend class ElementWrapper<init, name>;                   \
  friend class GElementBase<init, name>;                     \
  name(element_t e) : GElementBase(e) {}                     \
 public:                                                     \
  name(const Pairing& pairing) : GElementBase(pairing) {}    \
  name(const name& e) : GElementBase(e) {}                   \
};

__DEFINE_ELEMENT_CLASS(G1Element, element_init_G1);
__DEFINE_ELEMENT_CLASS(G2Element, element_init_G2);
__DEFINE_ELEMENT_CLASS(GTElement, element_init_GT);

#undef __DEFINE_ELEMENT_CLASS

/**
 * Wrapper for PBC Pairings. This makes it much easier to construct and manage
 * pairings and elements. Pairing is cleared automatically on destruction.
 */
class Pairing {
  template <__element_initializer_t init, class Derived>
  friend class ElementWrapper;

 private:
  // Just like elements, this needs to be mutable as PBC doesn't use const.
  mutable pairing_t pairing_;

 public:

  // Initialize pairing based on given params.
  Pairing(const string& params) {
    pairing_init_set_buf(pairing_, params.c_str(), params.size());
  }

  ~Pairing() {
    pairing_clear(pairing_);
  }

  // Returns true if the pairing is symmetric, false otherwise.
  bool is_symmetric() const {
    return pairing_is_symmetric(pairing_);
  }

  // Apply pairing to g and h and save result in gt. Useful if gt is already
  // allocated.
  void apply(GTElement& gt, const G1Element& g, const G2Element& h) const {
    pairing_apply(gt.v_, g.v_, h.v_, pairing_);
  }

  // Apply pairing to two elements of the first group. Throws invalid_argument
  // if the pairing is not symmetric.
  void apply(GTElement& gt, const G1Element& g, const G1Element& h) const {
    if (!is_symmetric()) {
      throw std::invalid_argument("pairing is not symmetric");
    }
    pairing_apply(gt.v_, g.v_, h.v_, pairing_);
  }

  // Apply pairing to two elements of the second group. Throws invalid_argument
  // if the pairing is not symmetric.
  void apply(GTElement& gt, const G2Element& g, const G2Element& h) const {
    if (!is_symmetric()) {
      throw std::invalid_argument("pairing is not symmetric");
    }
    pairing_apply(gt.v_, g.v_, h.v_, pairing_);
  }

  // Apply pairing to g and h and return the result.
  GTElement operator()(const G1Element& g, const G2Element& h) const {
    GTElement gt(*this);
    apply(gt, g, h);
    return gt;
  }

  // Apply pairing to two elements of the first group. Throws invalid_argument
  // if the pairing is not symmetric.
  GTElement operator()(const G1Element& g, const G1Element& h) const {
    GTElement gt(*this);
    apply(gt, g, h);
    return gt;
  }

  // Apply pairing to two elements of the second group. Throws invalid_argument
  // if the pairing is not symmetric.
  GTElement operator()(const G2Element& g, const G2Element& h) const {
    GTElement gt(*this);
    apply(gt, g, h);
    return gt;
  }
};

template <__element_initializer_t init, class Derived>
ElementWrapper<init, Derived>::ElementWrapper(const Pairing& pairing) {
  (*init)(v_, pairing.pairing_);
}

} // namespace pbc

#endif  // PBC_PBCXX_H
