#pragma once

#include <string>
#include <vector>

namespace mphys {

struct Field {
  std::string name;
  std::vector<double> values;

  Field() = default;
  Field(std::string name, int n_cells, double init_value = 0.0)
      : name(std::move(name)), values(n_cells, init_value) {}

  int NCells() const { return static_cast<int>(values.size()); }

  double& operator[](int i) { return values[i]; }
  double operator[](int i) const { return values[i]; }

  // Element-wise arithmetic
  Field operator+(const Field& rhs) const {
    Field result(name, NCells());
    for (int i = 0; i < NCells(); ++i) result[i] = values[i] + rhs[i];
    return result;
  }

  Field operator-(const Field& rhs) const {
    Field result(name, NCells());
    for (int i = 0; i < NCells(); ++i) result[i] = values[i] - rhs[i];
    return result;
  }

  Field operator*(double scalar) const {
    Field result(name, NCells());
    for (int i = 0; i < NCells(); ++i) result[i] = values[i] * scalar;
    return result;
  }

  Field operator-() const {
    Field result(name, NCells());
    for (int i = 0; i < NCells(); ++i) result[i] = -values[i];
    return result;
  }

  Field& operator+=(const Field& rhs) {
    for (int i = 0; i < NCells(); ++i) values[i] += rhs[i];
    return *this;
  }

  Field& operator-=(const Field& rhs) {
    for (int i = 0; i < NCells(); ++i) values[i] -= rhs[i];
    return *this;
  }
};

inline Field operator*(double scalar, const Field& f) { return f * scalar; }

}  // namespace mphys
