// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <vector>

namespace ActsExamples {

/// The simulated path of one particle: points along it, thinned to a minimum
/// distance, from its production to where the simulation stopped it.
struct SimParticleTrajectory {
  /// The particle
  SimBarcode particleId;
  /// Four-positions (x, y, z, t) of the points
  std::vector<Acts::Vector4> positions;
  /// Absolute momentum at each point
  std::vector<double> absMomenta;
};

/// The trajectories of all simulated particles of an event, ordered by
/// particle id
using SimParticleTrajectoryContainer = std::vector<SimParticleTrajectory>;

}  // namespace ActsExamples
