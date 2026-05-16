#pragma once

#include <vector>

#include "mphys/field.hpp"

namespace mphys {

struct Snapshot {
  double t = 0.0;
  std::vector<Field> fields;
  std::vector<double> algebraics;
};

struct SimResult {
  std::vector<Snapshot> snapshots;

  void Record(double t, const std::vector<Field>& fields,
              const std::vector<double>& algebraics) {
    snapshots.push_back({t, fields, algebraics});
  }
};

}  // namespace mphys
