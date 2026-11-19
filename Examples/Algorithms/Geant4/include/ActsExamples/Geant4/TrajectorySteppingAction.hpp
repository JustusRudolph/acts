// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"

#include <cstddef>
#include <memory>

#include <G4UserSteppingAction.hh>

class G4Step;

namespace ActsExamples::Geant4 {

/// Records the path of every simulated particle, in all volumes, into the
/// event store (EventStore::trajectories). Only reads the step, so it can run
/// next to the SensitiveSteppingAction without changing the hits.
///
/// To keep the output small, a point is only stored once it is at least
/// minStepLength away from the previous one; the production point and the
/// point where the particle stops (decay, interaction, kill, world exit) are
/// always stored.
class TrajectorySteppingAction : public G4UserSteppingAction {
 public:
  struct Config {
    /// event store
    std::shared_ptr<EventStore> eventStore;
    /// minimum distance between two stored points
    double minStepLength = 1 * Acts::UnitConstants::mm;
    /// at most this many points per particle (the end point is always kept)
    std::size_t maxPoints = 100000;
    /// record charged particles
    bool charged = true;
    /// record neutral particles
    bool neutral = true;
  };

  /// Construct the stepping action
  ///
  /// @param cfg the configuration struct
  /// @param logger the ACTS logging instance
  explicit TrajectorySteppingAction(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger =
          Acts::getDefaultLogger("TrajectorySteppingAction",
                                 Acts::Logging::INFO));

  /// Store the step's end point if it is far enough from the previous one, or
  /// if the particle stops there
  /// @param stepPtr is the Geant4 step of the particle
  void UserSteppingAction(const G4Step* stepPtr) override;

 private:
  EventStore& eventStore() const { return *m_cfg.eventStore; }

  const Acts::Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples::Geant4
