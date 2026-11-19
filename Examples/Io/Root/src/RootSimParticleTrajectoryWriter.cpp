// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootSimParticleTrajectoryWriter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <ios>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

namespace ActsExamples {

RootSimParticleTrajectoryWriter::RootSimParticleTrajectoryWriter(
    const RootSimParticleTrajectoryWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.inputTrajectories, "RootSimParticleTrajectoryWriter",
              level),
      m_cfg(config) {
  // inputTrajectories is already checked by base constructor
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing file path");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath + "'");
  }
  m_outputFile->cd();
  m_outputTree = new TTree(m_cfg.treeName.c_str(), m_cfg.treeName.c_str());
  if (m_outputTree == nullptr) {
    throw std::bad_alloc();
  }

  m_outputTree->Branch("event_id", &m_eventId);
  m_outputTree->Branch("vertex_primary", &m_vertexPrimary);
  m_outputTree->Branch("vertex_secondary", &m_vertexSecondary);
  m_outputTree->Branch("particle", &m_particle);
  m_outputTree->Branch("generation", &m_generation);
  m_outputTree->Branch("sub_particle", &m_subParticle);
  m_outputTree->Branch("step_x", &m_x);
  m_outputTree->Branch("step_y", &m_y);
  m_outputTree->Branch("step_z", &m_z);
  m_outputTree->Branch("step_t", &m_t);
  m_outputTree->Branch("step_p", &m_p);
}

RootSimParticleTrajectoryWriter::~RootSimParticleTrajectoryWriter() {
  if (m_outputFile != nullptr) {
    m_outputFile->Close();
  }
}

ProcessCode RootSimParticleTrajectoryWriter::finalize() {
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();

  ACTS_VERBOSE("Wrote particle trajectories to tree '"
               << m_cfg.treeName << "' in '" << m_cfg.filePath << "'");

  return ProcessCode::SUCCESS;
}

ProcessCode RootSimParticleTrajectoryWriter::writeT(
    const AlgorithmContext& ctx,
    const SimParticleTrajectoryContainer& trajectories) {
  // ensure exclusive access to tree/file while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_eventId = ctx.eventNumber;
  for (const auto& trajectory : trajectories) {
    const auto& barcode = trajectory.particleId;
    m_vertexPrimary = barcode.vertexPrimary();
    m_vertexSecondary = barcode.vertexSecondary();
    m_particle = barcode.particle();
    m_generation = barcode.generation();
    m_subParticle = barcode.subParticle();

    for (auto* column : {&m_x, &m_y, &m_z, &m_t, &m_p}) {
      column->clear();
      column->reserve(trajectory.positions.size());
    }
    for (std::size_t i = 0; i < trajectory.positions.size(); ++i) {
      const auto& position = trajectory.positions[i];
      m_x.push_back(static_cast<float>(position.x() / Acts::UnitConstants::mm));
      m_y.push_back(static_cast<float>(position.y() / Acts::UnitConstants::mm));
      m_z.push_back(static_cast<float>(position.z() / Acts::UnitConstants::mm));
      m_t.push_back(static_cast<float>(position.w() / Acts::UnitConstants::ns));
      m_p.push_back(static_cast<float>(trajectory.absMomenta[i] /
                                       Acts::UnitConstants::GeV));
    }
    m_outputTree->Fill();
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
