// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/TrajectorySteppingAction.hpp"

#include "ActsExamples/Geant4/AlgebraConverters.hpp"
#include "ActsExamples/Geant4/UnitConversion.hpp"

#include <cassert>
#include <cmath>
#include <utility>

#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>

namespace ActsExamples::Geant4 {

TrajectorySteppingAction::TrajectorySteppingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void TrajectorySteppingAction::UserSteppingAction(const G4Step* stepPtr) {
  assert(stepPtr != nullptr);
  const G4Step& step = *stepPtr;
  assert(step.GetTrack() != nullptr);
  const G4Track& track = *step.GetTrack();

  const double absCharge =
      std::abs(track.GetParticleDefinition()->GetPDGCharge());
  if ((absCharge > 0. && !m_cfg.charged) ||
      (absCharge == 0. && !m_cfg.neutral)) {
    return;
  }

  // only particles that made it into the particle record have a barcode
  const auto trackIt = eventStore().trackIdMapping.find(track.GetTrackID());
  if (trackIt == eventStore().trackIdMapping.end()) {
    return;
  }

  const G4StepPoint& preStepPoint = *step.GetPreStepPoint();
  const G4StepPoint& postStepPoint = *step.GetPostStepPoint();

  SimParticleTrajectory& trajectory =
      eventStore().trajectories[trackIt->second];
  auto store = [&trajectory](const G4StepPoint& point) {
    trajectory.positions.push_back(
        convertPosition(point.GetPosition(), point.GetGlobalTime()));
    trajectory.absMomenta.push_back(point.GetMomentum().mag() *
                                    convertEnergyToActs);
  };

  // the first step: its start is the production point
  if (trajectory.positions.empty()) {
    trajectory.particleId = trackIt->second;
    store(preStepPoint);
  }

  // the track status is set by the processes (and the ParticleKillAction
  // before this action) when the particle stops at the end of this step
  const bool isLastStep = track.GetTrackStatus() != fAlive ||
                          postStepPoint.GetStepStatus() == fWorldBoundary;
  const Acts::Vector3 position =
      convertPosition(postStepPoint.GetPosition());
  const bool isFarEnough =
      (position - trajectory.positions.back().head<3>()).norm() >=
      m_cfg.minStepLength;
  const bool hasRoom = trajectory.positions.size() + 1 < m_cfg.maxPoints;
  if (isLastStep || (isFarEnough && hasRoom)) {
    store(postStepPoint);
  }
}

}  // namespace ActsExamples::Geant4
