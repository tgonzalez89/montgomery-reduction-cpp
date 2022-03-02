// https://www.nayuki.io/page/montgomery-reduction-algorithm

#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>

uint32_t bit_length(uint32_t n)
{
  uint32_t result = 0;
  while (n > 0) {
    n >>= 1;
    ++result;
  }
  return result;
}

uint32_t mod(const int32_t x, const int32_t n)
{
  int32_t result = x % n;
  if (result < 0) {
    result += n;
  }
  return result;
}

// Based on a simplification of the extended Euclidean algorithm
uint32_t reciprocal_mod(const uint32_t n, const uint32_t r)
{
  uint32_t x = n;
  uint32_t y = r % n;
  int32_t a = 0;
  int32_t b = 1;
  // std::cout << "a=" << a << ", b=" << b << ", x=" << x << ", y=" << y << "\n";

  while (y != 0) {
    const auto tmp_b = b;
    b = a - static_cast<int32_t>(x / y) * b;
    a = tmp_b;

    const auto tmp_y = y;
    y = x % y;
    x = tmp_y;
    // std::cout << "a=" << a << ", b=" << b << ", x=" << x << ", y=" << y << "\n";
  }

  if (x != 1) {
    std::cout << "n=" << n << ", r=" << r << "\n";
    throw std::runtime_error("Reciprocal does not exist.");
  }

  return mod(a, n);
}

class Montgomery {
public:
  Montgomery(const uint32_t _n) : n(_n)
  {
    if (n < 3) {
      std::cout << "n=" << n << "\n";
      throw std::invalid_argument("Modulus must be >= 3.");
    }
    if (n % 2 == 0) {
      std::cout << "n=" << n << "\n";
      throw std::invalid_argument("Modulus must be odd.");
    }
    if (n > INT32_MAX) {
      std::cout << "n=" << n << "\n";
      throw std::invalid_argument("Modulus must be less than 2^31.");
    }

    r_bit_len = bit_length(n);
    assert(r_bit_len <= 31);

    uint32_t r = 1U << r_bit_len;
    r_mask = r - 1;
    r_reciprocal = reciprocal_mod(n, r); // r^-1 mod n

    k = (static_cast<uint64_t>(r) * static_cast<uint64_t>(r_reciprocal) - 1) / n;

    // std::cout << "r_bit_len=" << r_bit_len << "\n";
    // std::cout << "r=" << r << "\n";
    // std::cout << "r_reciprocal=" << r_reciprocal << "\n";
    // std::cout << "r_mask=" << r_mask << "\n";
    // std::cout << "k=" << k << "\n\n";
  }

  uint32_t convert_in(const uint32_t x)
  {
    return (static_cast<uint64_t>(x) << r_bit_len) % n;
  }

  uint32_t convert_out(const uint32_t x)
  {
    return (static_cast<uint64_t>(x) * static_cast<uint64_t>(r_reciprocal)) % n;
  }

  uint32_t multiply(const uint32_t a, const uint32_t b)
  {
    const uint64_t x = static_cast<uint64_t>(a) * static_cast<uint64_t>(b);
    const uint32_t s = ((x & r_mask) * static_cast<uint64_t>(k)) & r_mask;
    const uint64_t t = x + static_cast<uint64_t>(s) * static_cast<uint64_t>(n);
    uint32_t u = t >> r_bit_len;
    if (u >= n) {
      u -= n;
    }
    return u;
  }

private:
  uint32_t n;
  uint32_t r_bit_len;
  uint32_t r_reciprocal;
  uint32_t r_mask;
  uint32_t k;
};

int main()
{
  // uint32_t n = 1280541179;
  // uint32_t a = 1115177062;
  // uint32_t b = 95490452;
  // Montgomery mont(n);
  // uint32_t a_ = mont.convert_in(a);
  // uint32_t b_ = mont.convert_in(b);
  // uint32_t c_ = mont.multiply(a_, b_);
  // uint32_t c = mont.convert_out(c_);
  // uint32_t expected = static_cast<uint64_t>(a) * static_cast<uint64_t>(b) % n;
  // if (c != expected) {
  //   std::cout << "ERROR: results differ, expected=" << expected << ", result=" << c << "\n";
  // }
  // return 0;


  std::random_device rd;
  std::mt19937 gen(rd());
  for (uint32_t bitlen = 1; bitlen <= 30; ++bitlen) {
    std::cout << "bitlen=" << bitlen+1 << "\n";
    const uint32_t min_n = (1U << bitlen) + 1;
    const uint32_t max_n = UINT32_MAX >> (31 - bitlen);
    std::uniform_int_distribution<uint32_t> distr_n(min_n, max_n);
    for (size_t i = 0; i < 1000; ++i) {
      uint32_t n = 0;
      while (n % 2 == 0) {
        n = distr_n(gen);
      }
      Montgomery mont(n);
      std::uniform_int_distribution<uint32_t> distr_ops(0, n - 1);
      const uint32_t a = distr_ops(gen);
      const uint32_t b = distr_ops(gen);
      const uint32_t a_ = mont.convert_in(a);
      const uint32_t b_ = mont.convert_in(b);
      const uint32_t c_ = mont.multiply(a_, b_);
      const uint32_t c = mont.convert_out(c_);
      const uint32_t expected = static_cast<uint64_t>(a) * static_cast<uint64_t>(b) % n;
      if (c != expected) {
        std::cout << "res=" << c << ", ref=" << expected << "\n";
        std::cout << "a=" << a << ", b=" << b << ", n=" << n << "\n";
        throw std::runtime_error("Montgomery multiplication test failed.");
      }
    }
  }

  return 0;
}
