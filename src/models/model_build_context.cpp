#include "mphys/models/model_info.hpp"

#include <algorithm>
#include <utility>

namespace mphys {

ModelBuildContext::ModelBuildContext(const Mesh1D& mesh, StateVector& sv,
                                     std::vector<int> cell_domain,
                                     ParamLookup param_lookup, BcLookup bc_lookup)
    : mesh_(mesh),
      sv_(sv),
      cell_domain_(std::move(cell_domain)),
      param_lookup_(std::move(param_lookup)),
      bc_lookup_(std::move(bc_lookup)) {}

int ModelBuildContext::n_domains() const {
  if (cell_domain_.empty()) return 0;
  return *std::max_element(cell_domain_.begin(), cell_domain_.end()) + 1;
}

double ModelBuildContext::DomainParam(int domain, const std::string& key) const {
  return param_lookup_(domain, key);
}

std::vector<double> ModelBuildContext::CellParam(const std::string& key) const {
  std::vector<double> v(mesh_.n_cells);
  for (int i = 0; i < mesh_.n_cells; ++i)
    v[i] = param_lookup_(cell_domain_[i], key);
  return v;
}

Field ModelBuildContext::FaceParam(const std::string& key,
                                   FaceInterp interp) const {
  int n = mesh_.n_cells;
  Field face(key + "_face", n + 1);
  if (n == 0) return face;

  face[0] = param_lookup_(cell_domain_[0], key);
  face[n] = param_lookup_(cell_domain_[n - 1], key);
  for (int i = 1; i < n; ++i) {
    double vL   = param_lookup_(cell_domain_[i - 1], key);
    double vR   = param_lookup_(cell_domain_[i], key);
    double dx_L = mesh_.dx[i - 1];
    double dx_R = mesh_.dx[i];
    switch (interp) {
      case FaceInterp::kHarmonic:
        // Resistance-weighted harmonic mean: v = (dxL+dxR)/(dxL/vL + dxR/vR)
        face[i] = (dx_L + dx_R) / (dx_L / vL + dx_R / vR);
        break;
      case FaceInterp::kArithmetic:
        face[i] = (vR * dx_L + vL * dx_R) / (dx_L + dx_R);
        break;
      case FaceInterp::kNone:
        face[i] = 0.5 * (vL + vR);
        break;
    }
  }
  return face;
}

BcChoice ModelBuildContext::Bc(const std::string& slot) const {
  return bc_lookup_(slot);
}

}  // namespace mphys
