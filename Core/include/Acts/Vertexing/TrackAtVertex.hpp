// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

/// @class TrackAtVertex
///
/// @brief Defines a track at vertex object
///
/// @tparam input_track_t Track object type

template <typename input_track_t>
struct TrackAtVertex {
  /// Deleted default constructor
  TrackAtVertex() = delete;

  /// @brief Parameterized constructor
  ///
  /// @param chi2perTrack Chi2 of track
  /// @param paramsAtVertex Fitted perigee parameter
  /// @param originalParams Original perigee parameter
  TrackAtVertex(double chi2perTrack, const BoundParameters& paramsAtVertex,
                const input_track_t& originalTrack)
      : m_chi2Track(chi2perTrack),
        m_fittedParams(paramsAtVertex),
        m_originalTrack(originalTrack) {}

  /// Chi2 of track
  double m_chi2Track;

  /// Fitted perigee
  BoundParameters m_fittedParams;

  /// Original input track
  input_track_t m_originalTrack;
};

}  // namespace Acts