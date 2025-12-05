// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"

#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <ostream>
#include <stdexcept>
#include <utility>
#include <format>

namespace ActsExamples {

TruthTrackFinder::TruthTrackFinder(const Config& config,
                                   Acts::Logging::Level level)
    : IAlgorithm("TruthTrackFinder", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
}

ProcessCode TruthTrackFinder::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);

  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());

  std::string str_to_print = "";
  if (particles.size() != 1) {
    str_to_print = std::format("TRUTHTRACKFINDER: (my) event {} Created {} prototracks", ctx.eventNumber, tracks.size());
  }
  ACTS_VERBOSE("Create prototracks for " << particles.size() << " particles");
  for (const auto& particle : particles) {
    str_to_print += std::format("\tParticle {} has {} hits.", particle.particleId().value(), particle.numberOfHits());
    // find the corresponding hits for this particle
    const auto& measurements =
        makeRange(particleMeasurementsMap.equal_range(particle.particleId()));
    ACTS_VERBOSE(" - Prototrack from " << measurements.size()
                                       << " measurements");
    // fill hit indices to create the proto track
    ProtoTrack track;
    track.reserve(measurements.size());
    for (const auto& measurement : measurements) {
      track.emplace_back(measurement.second);
    }
    // add proto track to the output collection
    tracks.emplace_back(std::move(track));
  }
  // printf("%s\n", str_to_print.c_str());
  if (particles.size() != tracks.size()) {
    printf("TRUTHTRACKFINDER: MISMATCH in event %lu: Created %zu prototracks for %zu particles.\n",
      ctx.eventNumber, tracks.size(), particles.size());
  } else {
    printf("TRUTHTRACKFINDER: Created %zu prototracks for %zu particles in event %lu.\n",
      tracks.size(), particles.size(), ctx.eventNumber);
  }

  m_outputProtoTracks(ctx, std::move(tracks));
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
