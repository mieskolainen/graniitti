// Nuclear identity and normalized density model
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MNucleus.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace gra::nuclear {
namespace {

// Decode one candidate nuclear code without throwing
bool TryDecode(const int pdg, NucleusID &id) {
  const std::int64_t signed_code = static_cast<std::int64_t>(pdg);
  const std::int64_t code        = signed_code < 0 ? -signed_code : signed_code;
  if (code / 100000000 != 10) { return false; }

  id.pdg    = pdg;
  id.anti   = signed_code < 0;
  id.isomer = static_cast<unsigned int>(code % 10);
  id.a      = static_cast<unsigned int>((code / 10) % 1000);
  id.z      = static_cast<unsigned int>((code / 10000) % 1000);
  id.lambda = static_cast<unsigned int>((code / 10000000) % 10);
  return id.a > 0 && id.z <= id.a && id.lambda <= id.a - id.z;
}

}  // namespace

// Compute whether a PDG code has a valid 10LZZZAAAI nuclear encoding
bool IsNuclearPDG(const int pdg) {
  NucleusID id;
  return TryDecode(pdg, id);
}

// Decode and validate one 10LZZZAAAI nuclear PDG code
NucleusID DecodeNuclearPDG(const int pdg) {
  NucleusID id;
  if (!TryDecode(pdg, id)) { throw std::invalid_argument("DecodeNuclearPDG: invalid 10LZZZAAAI nuclear PDG code"); }
  return id;
}

// Encode one validated 10LZZZAAAI nuclear PDG identity
int EncodeNuclearPDG(const unsigned int a, const unsigned int z, const unsigned int lambda, const unsigned int isomer,
                     const bool anti) {
  if (a == 0 || a > 999 || z > 999 || z > a || lambda > 9 || lambda > a - z || isomer > 9) {
    throw std::invalid_argument("EncodeNuclearPDG: invalid nuclear identity");
  }
  const std::int64_t code = 1000000000LL + static_cast<std::int64_t>(lambda) * 10000000LL +
                            static_cast<std::int64_t>(z) * 10000LL + static_cast<std::int64_t>(a) * 10LL + isomer;
  if (code > std::numeric_limits<int>::max()) {
    throw std::overflow_error("EncodeNuclearPDG: code exceeds signed int");
  }
  return anti ? -static_cast<int>(code) : static_cast<int>(code);
}

// Compute tuned spherical density parameters for one encoded nucleus
// R = r_0 A^(1/3) when no isotope specific radius is given
NucleusParam DefaultNucleusParam(const int pdg, const double mass, const GeometryParam &geometry) {
  const NucleusID id = DecodeNuclearPDG(pdg);
  if (!std::isfinite(mass) || !(mass > 0.0)) {
    throw std::invalid_argument("DefaultNucleusParam: mass must be finite and positive");
  }

  NucleusParam param;
  param.pdg  = pdg;
  param.mass = mass;
  if (!std::isfinite(geometry.radius_scale) || !std::isfinite(geometry.skin) || !std::isfinite(geometry.tail_skin) ||
      !(geometry.radius_scale > 0.0) || !(geometry.skin > 0.0) || !(geometry.tail_skin > 0.0) || geometry.nodes < 32 ||
      geometry.nodes > 4096) {
    throw std::invalid_argument("DefaultNucleusParam: invalid tuned nuclear geometry");
  }

  double charge_radius = geometry.radius_scale * std::cbrt(static_cast<double>(id.a));
  double charge_skin   = geometry.skin;
  double matter_radius = charge_radius;
  double matter_skin   = geometry.skin;
  for (const auto &shape : geometry.nucleus) {
    if (shape.a == id.a && shape.z == id.z) {
      charge_radius = shape.charge[0];
      charge_skin   = shape.charge[1];
      matter_radius = shape.matter[0];
      matter_skin   = shape.matter[1];
      break;
    }
  }
  const auto density = [&geometry](const double radius, const double skin) {
    if (!std::isfinite(radius) || !std::isfinite(skin) || !(radius > 0.0) || !(skin > 0.0)) {
      throw std::invalid_argument("DefaultNucleusParam: invalid Fermi-density parameters");
    }
    return DensityParam{radius,
                        skin,
                        radius + geometry.tail_skin * skin,
                        geometry.nodes,
                        geometry.cdf_nodes,
                        geometry.form_q_max,
                        geometry.form_nodes,
                        geometry.form_abs_tol};
  };
  param.charge = density(charge_radius, charge_skin);
  param.matter = density(matter_radius, matter_skin);
  return param;
}

// Construct one nucleus from explicit physical inputs
MNucleus::MNucleus(const NucleusParam &param)
    : param_(param), id_(DecodeNuclearPDG(param.pdg)), charge_(param.charge), matter_(param.matter) {
  if (!std::isfinite(param_.mass) || !(param_.mass > 0.0)) {
    throw std::invalid_argument("MNucleus: mass must be finite and positive");
  }
  const double mass_per_nucleon = param_.mass / static_cast<double>(id_.a);
  if (!(mass_per_nucleon > 0.5) || !(mass_per_nucleon < 1.5)) {
    throw std::invalid_argument("MNucleus: mass per nucleon is outside the physical range");
  }
  param_.charge = charge_.Param();
  param_.matter = matter_.Param();
}

// Compute the neutron number for an ordinary nucleus
// N = A - Z - Lambda
unsigned int MNucleus::N() const { return id_.a - id_.z - id_.lambda; }

// Compute the signed electric charge in positron units
// Q/e = (-1)^anti Z
int MNucleus::Charge() const {
  const int z = static_cast<int>(id_.z);
  return id_.anti ? -z : z;
}

}  // namespace gra::nuclear
