#pragma once

#include <string>
#include <vector>

#include "mphys/materials/electrodes.hpp"
#include "mphys/materials/electrolytes.hpp"
#include "mphys/materials/material.hpp"

// ============================================================
// Materials database — registry and lookup.
//
// Built-in materials are identified by strongly-typed enums (ElectrodeId,
// ElectrolyteId) rather than strings: editor autocomplete lists every available
// material and a typo is a compile error instead of a runtime throw.  Each
// material still carries a human-readable `name` for display in a UI.
//
//   const auto& pos = Database::Electrode(ElectrodeId::kNmc811Chen2020);
//
// For building pickers, ElectrodeIds()/ElectrolyteIds() enumerate the catalogue
// (optionally filtered by domain) and ElectrodeName()/ElectrolyteName() give the
// display label for an id.
// ============================================================

namespace mphys::materials {

// Identifiers for the built-in electrode active materials.
enum class ElectrodeId {
  kGraphiteChen2020,
  kNmc811Chen2020,
  kLfpPrada2013,
  kNmc532Ecker2015,
};

// Identifiers for the built-in electrolytes.
enum class ElectrolyteId {
  kLipf6EcEmcNyman2008,
  kLipf6EcDmcMarquis2019,
};

class Database {
 public:
  // Resolve an identifier to its material. Built once, returned by reference.
  static const ElectrodeMaterial& Electrode(ElectrodeId id) {
    switch (id) {
      case ElectrodeId::kGraphiteChen2020: {
        static const ElectrodeMaterial m = Graphite_Chen2020();
        return m;
      }
      case ElectrodeId::kNmc811Chen2020: {
        static const ElectrodeMaterial m = Nmc811_Chen2020();
        return m;
      }
      case ElectrodeId::kLfpPrada2013: {
        static const ElectrodeMaterial m = Lfp_Prada2013();
        return m;
      }
      case ElectrodeId::kNmc532Ecker2015: {
        static const ElectrodeMaterial m = Nmc532_Ecker2015();
        return m;
      }
    }
    static const ElectrodeMaterial kUnknown;  // unreachable for valid enums
    return kUnknown;
  }

  static const ElectrolyteMaterial& Electrolyte(ElectrolyteId id) {
    switch (id) {
      case ElectrolyteId::kLipf6EcEmcNyman2008: {
        static const ElectrolyteMaterial m = Lipf6_EcEmc_Nyman2008();
        return m;
      }
      case ElectrolyteId::kLipf6EcDmcMarquis2019: {
        static const ElectrolyteMaterial m = Lipf6_EcDmc_Marquis2019();
        return m;
      }
    }
    static const ElectrolyteMaterial kUnknown;  // unreachable for valid enums
    return kUnknown;
  }

  // Catalogue of identifiers, optionally filtered to electrodes whose typical
  // domain matches `filter` (Domain::kAny returns everything).  Use this to
  // populate UI pickers — pair each id with ElectrodeName(id) for the label.
  static std::vector<ElectrodeId> ElectrodeIds(Domain filter = Domain::kAny) {
    static const ElectrodeId kAll[] = {
        ElectrodeId::kGraphiteChen2020,
        ElectrodeId::kNmc811Chen2020,
        ElectrodeId::kLfpPrada2013,
        ElectrodeId::kNmc532Ecker2015,
    };
    std::vector<ElectrodeId> ids;
    for (ElectrodeId id : kAll) {
      const Domain d = Electrode(id).typical_domain;
      if (filter == Domain::kAny || d == Domain::kAny || d == filter)
        ids.push_back(id);
    }
    return ids;
  }

  static std::vector<ElectrolyteId> ElectrolyteIds() {
    return {ElectrolyteId::kLipf6EcEmcNyman2008,
            ElectrolyteId::kLipf6EcDmcMarquis2019};
  }

  // Display labels.
  static const std::string& ElectrodeName(ElectrodeId id) {
    return Electrode(id).name;
  }
  static const std::string& ElectrolyteName(ElectrolyteId id) {
    return Electrolyte(id).name;
  }

  static std::vector<std::string> ElectrodeNames(Domain filter = Domain::kAny) {
    std::vector<std::string> names;
    for (ElectrodeId id : ElectrodeIds(filter)) names.push_back(ElectrodeName(id));
    return names;
  }
  static std::vector<std::string> ElectrolyteNames() {
    std::vector<std::string> names;
    for (ElectrolyteId id : ElectrolyteIds()) names.push_back(ElectrolyteName(id));
    return names;
  }
};

}  // namespace mphys::materials
