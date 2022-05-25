// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp" 

namespace ActsExamples {

  class MikadoSelector {
  public:
    MikadoSelector() = default;

    uint32_t acccept(Acts::MultiTrajectory::TrackStateProxy trackState) const {
      const auto& sourceLink =
        static_cast<const IndexSourceLink&>(trackState.uncalibrated());

      // std::cout << "measurement unique index: " << sourceLink.index() << std::endl;
      return sourceLink.index();      
    }
  };

} // namespace