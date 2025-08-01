// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <random>

namespace HepMC3 {
class GenEvent;
}

namespace ActsExamples {

/// @brief Generator interface for vertices and particles
struct ParticlesGenerator {
  /// @brief Virtual destructor required
  virtual ~ParticlesGenerator() = default;
  /// @brief Generate vertices and particles for one interaction
  /// @note This method cannot be `const` because the Pythia8 generator
  ///       uses the Pythia8 interfaces, which is non-const
  ///
  /// @param rng Shared random number generator instance
  /// @return The vertex and particle containers
  virtual std::shared_ptr<HepMC3::GenEvent> operator()(RandomEngine& rng) = 0;
};

/// Generate particles from uniform parameter distributions.
///
/// Generates a single vertex with the given number of tracks. The track
/// direction is drawn from a uniform distribution on the unit sphere (within
/// the given limits). Its absolute momentum is drawn from a uniform
/// distribution. Position and time are always set to zero.
class ParametricParticleGenerator : public ParticlesGenerator {
 public:
  struct Config {
    /// Low, high (exclusive) for the transverse direction angle.
    double phiMin = -std::numbers::pi;
    double phiMax = std::numbers::pi;
    /// Low, high (inclusive) for  the longitudinal direction angle.
    ///
    /// This intentionally uses theta instead of eta so it can represent the
    /// full direction space with finite values.
    ///
    /// @note This is the standard generation, for detector performance
    /// classification, where a flat distribution in eta can be useful,
    /// this can be set by the etaUniform flag;
    ///
    double thetaMin = std::numeric_limits<double>::min();
    double thetaMax = std::numbers::pi - std::numeric_limits<double>::epsilon();
    bool etaUniform = false;
    /// Low, high (exclusive) for absolute/transverse momentum.
    double pMin = 1 * Acts::UnitConstants::GeV;
    double pMax = 10 * Acts::UnitConstants::GeV;
    /// Indicate if the momentum refers to transverse momentum.
    bool pTransverse = false;
    /// Indicate if the momentum should be uniformly distributed in log space.
    bool pLogUniform = false;
    /// (Absolute) PDG particle number to identify the particle type.
    Acts::PdgParticle pdg = Acts::PdgParticle::eMuon;
    /// Randomize the charge and flip the PDG particle number sign accordingly.
    bool randomizeCharge = false;
    /// Number of particles.
    std::size_t numParticles = 1;

    /// Overrides particle charge.
    std::optional<double> charge;
    /// Overrides particle mass.
    std::optional<double> mass;
  };

  explicit ParametricParticleGenerator(const Config& cfg);

  /// Generate a single primary vertex with the given number of particles.
  std::shared_ptr<HepMC3::GenEvent> operator()(RandomEngine& rng) override;

 private:
  using UniformIndex = std::uniform_int_distribution<std::uint8_t>;
  using UniformReal = std::uniform_real_distribution<double>;

  Config m_cfg;

  // will be automatically set from PDG data tables
  double m_mass{};

  // (anti-)particle choice is one random draw but defines two properties
  std::array<Acts::PdgParticle, 2> m_pdgChoices{};

  UniformIndex m_particleTypeChoice;
  UniformReal m_phiDist;
  std::function<std::pair<double, double>(RandomEngine& rng)> m_sinCosThetaDist;
  std::function<double(RandomEngine& rng)> m_somePDist;
};

}  // namespace ActsExamples
