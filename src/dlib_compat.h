// dlib_compat.h -- dlib's optimizer as this package includes it.
//
// Marked a system header, as -isystem would mark dlib's include directory, so
// that warnings from dlib's own headers, which the package cannot change, are
// not reported. The dlib package's headers (dlib 19.2) keep their unicode
// strings as std::basic_string<unsigned int>, whose character traits libc++
// no longer supplies from version 19 on, so a specialization is given here
// for libc++; libstdc++ still has the general template. The std::auto_ptr
// and std::binary_function they also use, which libc++ drops under C++17,
// are restored by src/Makevars.
//
// dlib's always-on checks throw a dlib::fatal_error, which writes to
// std::cerr and can call abort(). Compiled code in an R package may do
// neither, so those checks throw an ordinary exception here instead, which
// cpp11 turns into an R error.
#pragma once
#pragma GCC system_header

#include <cstddef>
#include <cstring>
#include <ios>
#include <stdexcept>
#include <string>

#if defined(_LIBCPP_VERSION) && _LIBCPP_VERSION >= 190000
namespace std {
template <>
struct char_traits<unsigned int> {
  using char_type = unsigned int;
  using int_type = unsigned long;
  using off_type = streamoff;
  using pos_type = streampos;
  using state_type = mbstate_t;
  static void assign(char_type& a, const char_type& b) noexcept { a = b; }
  static constexpr bool eq(char_type a, char_type b) noexcept { return a == b; }
  static constexpr bool lt(char_type a, char_type b) noexcept { return a < b; }
  static int compare(const char_type* a, const char_type* b, size_t n) {
    for (size_t i = 0; i < n; ++i) {
      if (a[i] < b[i]) return -1;
      if (b[i] < a[i]) return 1;
    }
    return 0;
  }
  static size_t length(const char_type* s) {
    size_t n = 0;
    while (s[n] != 0) ++n;
    return n;
  }
  static const char_type* find(const char_type* s, size_t n,
                               const char_type& a) {
    for (size_t i = 0; i < n; ++i) {
      if (s[i] == a) return s + i;
    }
    return nullptr;
  }
  static char_type* move(char_type* d, const char_type* s, size_t n) {
    return static_cast<char_type*>(std::memmove(d, s, n * sizeof(char_type)));
  }
  static char_type* copy(char_type* d, const char_type* s, size_t n) {
    return static_cast<char_type*>(std::memcpy(d, s, n * sizeof(char_type)));
  }
  static char_type* assign(char_type* s, size_t n, char_type a) {
    for (size_t i = 0; i < n; ++i) s[i] = a;
    return s;
  }
  static constexpr int_type not_eof(int_type c) noexcept {
    return c == eof() ? 0 : c;
  }
  static constexpr char_type to_char_type(int_type c) noexcept {
    return char_type(c);
  }
  static constexpr int_type to_int_type(char_type c) noexcept {
    return int_type(c);
  }
  static constexpr bool eq_int_type(int_type a, int_type b) noexcept {
    return a == b;
  }
  static constexpr int_type eof() noexcept { return int_type(-1); }
};
}  // namespace std
#endif

#include <dlib/assert.h>

#undef DLIBM_CASSERT
#define DLIBM_CASSERT(_exp, _message)                                      \
  {                                                                        \
    if (!(_exp)) {                                                         \
      throw std::invalid_argument(std::string("dlib's check failed: ") +   \
                                  #_exp);                                  \
    }                                                                      \
  }

#include <dlib/optimization/optimization.h>
