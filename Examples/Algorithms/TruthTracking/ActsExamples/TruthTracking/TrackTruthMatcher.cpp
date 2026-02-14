// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TrackTruthMatcher.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace ActsExamples {

TrackTruthMatcher::TrackTruthMatcher(const Config& config,
                                     Acts::Logging::Level level)
    : IAlgorithm("TrackTruthMatcher", level), m_cfg(config) {
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input measurement particles map");
  }
  if (m_cfg.outputTrackParticleMatching.empty()) {
    throw std::invalid_argument("Missing output track particles matching");
  }
  if (m_cfg.outputParticleTrackMatching.empty()) {
    throw std::invalid_argument("Missing output particle track matching");
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_outputTrackParticleMatching.initialize(m_cfg.outputTrackParticleMatching);
  m_outputParticleTrackMatching.initialize(m_cfg.outputParticleTrackMatching);
}

ProcessCode TrackTruthMatcher::execute(const AlgorithmContext& ctx) const {
  // Read input tracks
  const auto& tracks = m_inputTracks(ctx);

  // Read truth input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& hitParticlesMap = m_inputMeasurementParticlesMap(ctx);

  TrackParticleMatching trackParticleMatching;
  ParticleTrackMatching particleTrackMatching;

  // TODO this may be computed in a separate algorithm
  // TODO can we wire this through?
  std::map<SimBarcode, std::size_t> particleTruthHitCount;
  for (const auto& [_, pid] : hitParticlesMap) {
    particleTruthHitCount[pid]++;
  }

  // For each particle within a track, how many hits did it contribute
  std::vector<ParticleHitCount> particleHitCounts;
  std::set<unsigned> interesting_events = {16, 374, 429, 813, 927};
  if ( true || (std::find(interesting_events.begin(), interesting_events.end(),
                  ctx.eventNumber) != interesting_events.end() ) ) {
    ACTS_INFO("event " << ctx.eventNumber << " has "
                        << particles.size() << " particles and "
                        << tracks.size() << " tracks");
  }
  for (const auto& track : tracks) {
    // Get the majority truth particle to this track
    identifyContributingParticles(hitParticlesMap, track, particleHitCounts);
    if (particleHitCounts.empty()) {
      ACTS_DEBUG(
          "No truth particle associated with this trajectory with tip index = "
          << track.tipIndex());
      continue;
    }

    // Get the majority particleId and majority particle counts
    // Note that the majority particle might not be in the truth seeds
    // collection
    ActsFatras::Barcode majorityParticleId =
        particleHitCounts.front().particleId;
    std::size_t nMajorityHits = particleHitCounts.front().hitCount;

    if (!particles.contains(majorityParticleId)) {
      ACTS_VERBOSE(
          "The majority particle is not in the input particle collection, "
          "majorityParticleId = "
          << majorityParticleId);
      continue;
    }

    // Check if the trajectory is matched with truth.
    // If not, it will be classified as 'fake'
    const bool recoMatched =
        static_cast<double>(nMajorityHits) / track.nMeasurements() >=
        m_cfg.matchingRatio;
    const double trackTruthRatio = static_cast<double>(nMajorityHits) /
                                   particleTruthHitCount.at(majorityParticleId);
    const bool truthMatched = trackTruthRatio >= m_cfg.matchingRatio;
    const auto trackWithWeight = std::make_pair(track.index(), trackTruthRatio);

    if ((!m_cfg.doubleMatching && recoMatched) ||
        (m_cfg.doubleMatching && recoMatched && truthMatched)) {
      auto& trackParticleMatch = trackParticleMatching[track.index()] = {
          TrackMatchClassification::Matched, majorityParticleId,
          particleHitCounts};

      auto& particleTrackMatch = particleTrackMatching[majorityParticleId];
      if (!particleTrackMatch.track) {
        particleTrackMatch.track = trackWithWeight;
      } else {
        // we already have a track associated with this particle and have to
        // resolve the ambiguity.
        // we will use the track with more hits and smaller chi2
        const auto& otherTrack =
            tracks.getTrack(particleTrackMatch.track.value().first);
        if (otherTrack.nMeasurements() < track.nMeasurements() ||
            otherTrack.chi2() > track.chi2()) {
          trackParticleMatching[otherTrack.index()].classification =
              TrackMatchClassification::Duplicate;
          // previous now goes to duplicate since it has a worse match
          particleTrackMatch.duplicateIdxs.push_back(
              particleTrackMatch.track.value());
          // assign the new track as the main matched track
          particleTrackMatch.track = trackWithWeight;
        } else {
          trackParticleMatch.classification =
              TrackMatchClassification::Duplicate;
          particleTrackMatch.duplicateIdxs.push_back(trackWithWeight);
        }

        ++particleTrackMatch.duplicates;
      }
    } else {  // not matched, i.e. fake track
      ACTS_DEBUG("Track " << track.tipIndex() << " in event " << ctx.eventNumber
                 << " NOT MATCHED to particle " << majorityParticleId.particle()
                 << " with " << particleTruthHitCount.at(majorityParticleId) 
                 << " hits in event " << ctx.eventNumber << ". Out of " 
                 << track.nMeasurements() << " track hits, " 
                 << nMajorityHits << " were right.");
      trackParticleMatching[track.index()] = {TrackMatchClassification::Fake,
                                              std::nullopt, particleHitCounts};

      auto& particleTrackMatch = particleTrackMatching[majorityParticleId];
      particleTrackMatch.fakeIdxs.push_back(trackWithWeight);
      ++particleTrackMatch.fakes;
    }
  }
  // Sort duplicate and fakes by decreasing weight
  for (auto& [_, entry] : particleTrackMatching) {
    std::sort(entry.duplicateIdxs.begin(), entry.duplicateIdxs.end(),
              [](const TrackIndexWithWeight& t1, const TrackIndexWithWeight& t2) {
                return t1.second > t2.second;  // sort by weight
              });
    std::sort(entry.fakeIdxs.begin(), entry.fakeIdxs.end(),
              [](const TrackIndexWithWeight& t1, const TrackIndexWithWeight& t2) {
                return t1.second > t2.second;
              });
  }

  m_outputTrackParticleMatching(ctx, std::move(trackParticleMatching));
  m_outputParticleTrackMatching(ctx, std::move(particleTrackMatching));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
