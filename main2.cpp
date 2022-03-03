// https://www.nayuki.io/page/montgomery-reduction-algorithm

#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <cstdio>

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

// Calculate the modular multiplicative inverse
// Based on a simplification of the extended Euclidean algorithm
uint32_t mod_mult_inv(const uint32_t n, const uint32_t r)
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

/// @brief Hensel's Lemma for 2-adic numbers
/// Find solution for qX + 1 = 0 mod 2^r
/// @param[in] r
/// @param[in] q such that gcd(2, q) = 1
/// @return Unsigned long int in [0, 2^r − 1] such that q*x ≡ −1 mod 2^r
uint64_t HenselLemma2adicRoot(uint32_t r, uint64_t q) {
  uint64_t a_prev = 1;
  uint64_t c = 2;
  uint64_t mod_mask = 3;

  // Root:
  //    f(x) = qX + 1 and a_(0) = 1 then f(1) ≡ 0 mod 2
  // General Case:
  //    - a_(n) ≡ a_(n-1) mod 2^(n)
  //      => a_(n) = a_(n-1) + 2^(n)*t
  //    - Find 't' such that f(a_(n)) = 0 mod  2^(n+1)
  // First case in for:
  //    - a_(1) ≡ 1 mod 2 or a_(1) = 1 + 2t
  //    - Find 't' so f(a_(1)) ≡ 0 mod 4  => q(1 + 2t) + 1 ≡ 0 mod 4
  for (uint64_t k = 2; k <= r; k++) {
    uint64_t f = 0;
    uint64_t t = 0;
    uint64_t a = 0;

    do {
      a = a_prev + c * t++;
      f = q * a + 1ULL;
    } while (f & mod_mask);  // f(a) ≡ 0 mod 2^(k)

    // Update vars
    mod_mask = mod_mask * 2 + 1ULL;
    c *= 2;
    a_prev = a;
  }

  return a_prev;
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
    r_inv_mod = mod_mult_inv(n, r); // r^-1 mod n

    n_inv_mod = (static_cast<uint64_t>(r) * static_cast<uint64_t>(r_inv_mod) - 1) / n;

    r2_mod_n = (static_cast<uint64_t>(r) * static_cast<uint64_t>(r)) % n;

    // std::cout << "r_bit_len=" << r_bit_len << "\n";
    // std::cout << "r=" << r << "\n";
    // std::cout << "r_inv_mod=" << r_inv_mod << "\n";
    // std::cout << "r_mask=" << r_mask << "\n";
    // std::cout << "k=" << k << "\n\n";
  }

  uint32_t convert_in(uint32_t x)
  {
    if (x >= n) {
      x %= n; // Maybe this should be a warning or an error...
    }
    return REDC(static_cast<uint64_t>(x) * static_cast<uint64_t>(r2_mod_n));
  }

  uint32_t convert_out(const uint32_t x)
  {
    return REDC(x);
  }

  uint32_t multiply(const uint32_t a, const uint32_t b)
  {
    return REDC(static_cast<uint64_t>(a) * static_cast<uint64_t>(b));
  }

  uint32_t REDC(const uint64_t x)
  {
    const uint32_t s = ((x & r_mask) * static_cast<uint64_t>(n_inv_mod)) & r_mask;
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
  uint32_t r_inv_mod;
  uint32_t r_mask;
  uint32_t n_inv_mod;
  uint32_t r2_mod_n;
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

  // int32_t n1 = 2345;
  // int32_t bl = bit_length(n1);
  // int32_t n2 = 1U << bl;
  // int32_t n1_ = 0;
  // int32_t n2_ = 0;
  // int32_t gdc = gcdExtended(n1, n2, &n1_, &n2_);
  // int32_t n2_r = mod_mult_inv(n1, n2);
  // uint64_t n2_r2 = HenselLemma2adicRoot(bl, n1);
  // int32_t k = (static_cast<int64_t>(n2) * static_cast<int64_t>(n2_r) - 1) / n1;

  // printf("n1=%d, n2=%d, *n1_=%d, *n2_=%d, gdc=%d, n2_r=%d, n2_r2=%lu, k=%d\n", n1, n2, n1_, n2_, gdc, n2_r, n2_r2, k);

  return 0;
}
