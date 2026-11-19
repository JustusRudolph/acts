// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticleTrajectory.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class TFile;
class TTree;

namespace ActsExamples {

/// Write out the simulated particle trajectories as a TTree.
///
/// Each entry is one particle: the event number, the particle barcode (as in
/// the RootParticleWriter) and the points along its path as vectors.
///
/// Safe to use from multiple writer threads. To avoid thread-safety issues,
/// the writer must be the sole owner of the underlying file. Thus, the
/// output file pointer can not be given from the outside.
class RootSimParticleTrajectoryWriter final
    : public WriterT<SimParticleTrajectoryContainer> {
 public:
  struct Config {
    /// Input trajectory collection to write.
    std::string inputTrajectories;
    /// Path to the output file.
    std::string filePath;
    /// Output file access mode.
    std::string fileMode = "RECREATE";
    /// Name of the tree within the output file.
    std::string treeName = "trajectories";
  };

  /// Construct the trajectory writer.
  ///
  /// @param config is the configuration object
  /// @param level is the logging level
  RootSimParticleTrajectoryWriter(const Config& config,
                                  Acts::Logging::Level level);

  /// Ensure underlying file is closed.
  ~RootSimParticleTrajectoryWriter() override;

  /// End-of-run hook
  ProcessCode finalize() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 protected:
  /// Type-specific write implementation.
  ///
  /// @param[in] ctx is the algorithm context
  /// @param[in] trajectories are the trajectories to be written
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const SimParticleTrajectoryContainer& trajectories) override;

 private:
  Config m_cfg;
  std::mutex m_writeMutex;
  TFile* m_outputFile = nullptr;
  TTree* m_outputTree = nullptr;
  /// Event identifier.
  std::uint32_t m_eventId = 0;
  /// Decoded barcode components.
  std::uint32_t m_vertexPrimary = 0;
  std::uint32_t m_vertexSecondary = 0;
  std::uint32_t m_particle = 0;
  std::uint32_t m_generation = 0;
  std::uint32_t m_subParticle = 0;
  /// Points along the path: position in mm, time in ns, |p| in GeV.
  std::vector<float> m_x, m_y, m_z, m_t, m_p;
};

}  // namespace ActsExamples
