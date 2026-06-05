#pragma once

// Small arithmetic expression evaluator used by the numeric input widgets.
// Supports +,-,*,/,^,(,), unary minus, constants pi/e and common math
// functions (sqrt, sin, cos, tan, log/ln, log10, exp, abs, floor, ceil).

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

namespace expr {

struct Parser {
  const char* p;
  void skip() { while (*p == ' ' || *p == '\t') ++p; }

  double primary() {
    skip();
    if (*p == '-') { ++p; return -primary(); }
    if (*p == '+') { ++p; return  primary(); }
    if (*p == '(') {
      ++p;
      double v = expr_val();
      skip(); if (*p == ')') ++p;
      return v;
    }
    if (std::isalpha((unsigned char)*p)) {
      const char* s = p;
      while (std::isalnum((unsigned char)*p) || *p == '_') ++p;
      std::string name(s, p);
      skip();
      if (*p == '(') {
        ++p;
        double a = expr_val();
        skip(); if (*p == ')') ++p;
        if (name=="sqrt")            return std::sqrt(a);
        if (name=="sin")             return std::sin(a);
        if (name=="cos")             return std::cos(a);
        if (name=="tan")             return std::tan(a);
        if (name=="log"||name=="ln") return std::log(a);
        if (name=="log10")           return std::log10(a);
        if (name=="exp")             return std::exp(a);
        if (name=="abs")             return std::abs(a);
        if (name=="floor")           return std::floor(a);
        if (name=="ceil")            return std::ceil(a);
        return std::numeric_limits<double>::quiet_NaN();
      }
      if (name=="pi"||name=="PI") return 3.14159265358979323846;
      if (name=="e" ||name=="E")  return 2.71828182845904523536;
      return std::numeric_limits<double>::quiet_NaN();
    }
    char* end;
    double v = std::strtod(p, &end);
    if (end == p) return std::numeric_limits<double>::quiet_NaN();
    p = end;
    return v;
  }

  double power() {
    double b = primary();
    skip();
    if (*p == '^') { ++p; return std::pow(b, power()); }
    return b;
  }

  double term() {
    double v = power();
    for (;;) {
      skip();
      if      (*p == '*') { ++p; v *= power(); }
      else if (*p == '/') { ++p; double d = power(); v = d != 0.0 ? v/d : std::numeric_limits<double>::quiet_NaN(); }
      else break;
    }
    return v;
  }

  double expr_val() {
    double v = term();
    for (;;) {
      skip();
      if      (*p == '+') { ++p; v += term(); }
      else if (*p == '-') { ++p; v -= term(); }
      else break;
    }
    return v;
  }
};

inline double eval(const char* s) {
  if (!s || !*s) return 0.0;
  Parser pr{s};
  return pr.expr_val();
}

}  // namespace expr
