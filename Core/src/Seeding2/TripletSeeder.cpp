// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/TripletSeeder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"
#include "Acts/Seeding2/TripletSeedFinder.hpp"

#include <Eigen/Dense>

namespace Acts {

namespace {

template <typename DoubletCollections>
void createAndFilterTriplets(TripletSeeder::Cache& cache,
                             const TripletSeedFinder& tripletFinder,
                             const ITripletSeedFilter& filter,
                             const SpacePointContainer2& spacePoints,
                             DoubletCollections bottomDoublets,
                             const ConstSpacePointProxy2& spM,
                             DoubletCollections topDoublets,
                             unsigned evNo = 36) {
  // if (evNo == 36) {
  //   std::cout << "TRIPLET SEEDER: Event " << evNo
  //             << " createAndFilterTriplets called with "
  //             << bottomDoublets.size() << " bottom doublets and "
  //             << topDoublets.size() << " top doublets" << std::endl;
  // }
  for (auto bottomDoublet : bottomDoublets) {
    if (topDoublets.empty()) {
      break;
    }

    cache.tripletTopCandidates.clear();
    tripletFinder.createTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                             topDoublets,
                                             cache.tripletTopCandidates, evNo);
      // if (evNo == 36) {
      //   std::cout << "TRIPLET SEEDER: Event " << evNo
      //             << " Created " << cache.tripletTopCandidates.size()
      //             << " triplet top candidates" << std::endl;
      // }

    filter.filterTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      cache.tripletTopCandidates);
    // if (evNo == 36) {
    //   std::cout << "TRIPLET SEEDER: Event " << evNo << " After filtering: "
    //             << cache.tripletTopCandidates.size()
    //             << " triplet top candidates remain" << std::endl;
    // }
  }
}

template <typename SpacePointCollections>
void createSeedsFromGroupsImpl(
    const Logger& logger, TripletSeeder::Cache& cache,
    const DoubletSeedFinder& bottomFinder, const DoubletSeedFinder& topFinder,
    const TripletSeedFinder& tripletFinder, const ITripletSeedFilter& filter,
    const SpacePointContainer2& spacePoints,
    SpacePointCollections& bottomSpGroups,
    const ConstSpacePointProxy2& middleSp, SpacePointCollections& topSpGroups,
    SeedContainer2& outputSeeds, unsigned evNo=36) {
  MiddleSpInfo middleSpInfo = DoubletSeedFinder::computeMiddleSpInfo(middleSp);

  // create middle-top doublets
  cache.topDoublets.clear();
  for (auto& topSpGroup : topSpGroups) {
    // if (evNo == 36) {
    //   ACTS_INFO("TRIPLET SEEDER: Event " << evNo
    //               << " creating top doublets for middle candidate indexed "
    //               << middleSp.index() << " with coordinates (" << middleSp.xy()[0]
    //               << ", " << middleSp.xy()[1] << ", " << middleSp.zr()[0]
    //               << ") and radius " << middleSp.zr()[1]);
    // }
    topFinder.createDoublets(middleSp, middleSpInfo, topSpGroup,
                             cache.topDoublets, evNo);
  }

  // no top SP found -> cannot form any triplet
  // if (evNo == 36) {
  //   ACTS_INFO("TRIPLET SEEDER: Event " << evNo << " Created "
  //               << cache.topDoublets.size() << " top doublets for middle SP "
  //               << middleSp.index());
  // }
  if (cache.topDoublets.empty()) {
    // if (evNo == 36) {
    //   ACTS_INFO("Event " << evNo
    //               << " No compatible Tops for middle candidate indexed "
    //               << middleSp.index() << ", returning");
    // }
    ACTS_VERBOSE("No compatible Tops, returning");
    return;
  }

  if (!filter.sufficientTopDoublets(spacePoints, middleSp, cache.topDoublets)) {
    // if (evNo == 36) {
    //   ACTS_INFO("TRIPLET SEEDER: Event " << evNo
    //               << " Insufficient top doublets for middle SP "
    //               << middleSp.index());
    // }
    return;
  }

  // create middle-bottom doublets
  cache.bottomDoublets.clear();
  for (auto& bottomSpGroup : bottomSpGroups) {
    // if (evNo == 36) {
    //   ACTS_INFO("TRIPLET SEEDER: Event " << evNo
    //               << " creating bottom doublets for middle candidate indexed "
    //               << middleSp.index() << " with coordinates (" << middleSp.xy()[0]
    //               << ", " << middleSp.xy()[1] << ", " << middleSp.zr()[0]
    //               << ") and radius " << middleSp.zr()[1]);
    // }
    bottomFinder.createDoublets(middleSp, middleSpInfo, bottomSpGroup,
                                cache.bottomDoublets, evNo);
  }

  // no bottom SP found -> cannot form any triplet
  // if (evNo == 36) {
  //   ACTS_INFO("TRIPLET SEEDER: Event " << evNo << " Created "
  //               << cache.bottomDoublets.size()
  //               << " bottom doublets for middle SP " << middleSp.index());
  // }
  if (cache.bottomDoublets.empty()) {
    // if (evNo == 36) {
    //   ACTS_INFO("Event " << evNo
    //               << " No compatible Bottoms for middle candidate indexed "
    //               << middleSp.index() << ", returning");
    // }
    ACTS_VERBOSE("No compatible Bottoms, returning");
    return;
  }

  // if (evNo == 36) {
  //   ACTS_INFO("TRIPLET SEEDER: Event " << evNo << " Combining "
  //               << cache.bottomDoublets.size() << " bottom doublets with "
  //               << cache.topDoublets.size() << " top doublets for middle SP "
  //               << middleSp.index());
  // }
  ACTS_VERBOSE("Candidates: " << cache.bottomDoublets.size() << " bottoms and "
                              << cache.topDoublets.size()
                              << " tops for middle candidate indexed "
                              << middleSp.index());

  // combine doublets to triplets
  if (tripletFinder.config().sortedByCotTheta) {
    cache.bottomDoublets.sortByCotTheta({0, cache.bottomDoublets.size()},
                                        cache.sortedBottoms);
    cache.topDoublets.sortByCotTheta({0, cache.topDoublets.size()},
                                     cache.sortedTops);

    createAndFilterTriplets(cache, tripletFinder, filter, spacePoints,
                            cache.bottomDoublets.subset(cache.sortedBottoms),
                            middleSp,
                            cache.topDoublets.subset(cache.sortedTops), evNo);
  } else {
    createAndFilterTriplets(cache, tripletFinder, filter, spacePoints,
                            cache.bottomDoublets.range(), middleSp,
                            cache.topDoublets.range(), evNo);
  }

  // std::size_t seedsBefore = outputSeeds.size();
  filter.filterTripletsMiddleFixed(spacePoints, outputSeeds);
  // if (evNo == 36) {
  //   std::cout << "TRIPLET SEEDER: Event " << evNo
  //             << " After filterTripletsMiddleFixed: added "
  //             << (outputSeeds.size() - seedsBefore) << " seeds (total now: "
  //             << outputSeeds.size() << ")" << std::endl;
  // }
}

}  // namespace

TripletSeeder::TripletSeeder(std::unique_ptr<const Logger> logger_)
    : m_logger(std::move(logger_)) {
  if (m_logger == nullptr) {
    throw std::invalid_argument("TripletSeeder: logger cannot be null");
  }
}

void TripletSeeder::createSeedsFromGroup(
    Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const ITripletSeedFilter& filter, const SpacePointContainer2& spacePoints,
    SpacePointContainer2::ConstSubset& bottomSps,
    const ConstSpacePointProxy2& middleSp,
    SpacePointContainer2::ConstSubset& topSps,
    SeedContainer2& outputSeeds, unsigned evNo) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");

  std::array<SpacePointContainer2::ConstSubset, 1> bottomSpGroups{bottomSps};
  std::array<SpacePointContainer2::ConstSubset, 1> topSpGroups{topSps};

  createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                            tripletFinder, filter, spacePoints, bottomSpGroups,
                            middleSp, topSpGroups, outputSeeds, evNo);
}

void TripletSeeder::createSeedsFromGroups(
    Cache& cache, const DoubletSeedFinder& bottomFinder,
    const DoubletSeedFinder& topFinder, const TripletSeedFinder& tripletFinder,
    const ITripletSeedFilter& filter, const SpacePointContainer2& spacePoints,
    const std::span<SpacePointContainer2::ConstRange>& bottomSpGroups,
    const SpacePointContainer2::ConstRange& middleSpGroup,
    const std::span<SpacePointContainer2::ConstRange>& topSpGroups,
    const std::pair<float, float>& radiusRangeForMiddle,
    SeedContainer2& outputSeeds, unsigned evNo) const {
  assert((bottomFinder.config().spacePointsSortedByRadius ==
          topFinder.config().spacePointsSortedByRadius) &&
         "Inconsistent space point sorting");
  const bool spacePointsSortedByRadius =
      bottomFinder.config().spacePointsSortedByRadius;

  if (middleSpGroup.empty()) {
    return;
  }

  if (spacePointsSortedByRadius) {
    // Initialize initial offsets for bottom and top space points with binary
    // search. This requires at least one middle space point to be present which
    // is already checked above.
    const ConstSpacePointProxy2 firstMiddleSp = middleSpGroup.front();
    const float firstMiddleSpR = firstMiddleSp.zr()[1];

    for (auto& bottomSpGroup : bottomSpGroups) {
      // Find the first bottom space point that is within the deltaRMax of the
      // first middle space point.
      const auto low = std::ranges::lower_bound(
          bottomSpGroup, firstMiddleSpR - bottomFinder.config().deltaRMax, {},
          [&](const ConstSpacePointProxy2& sp) { return sp.zr()[1]; });
      bottomSpGroup = bottomSpGroup.subrange(low - bottomSpGroup.begin());
    }

    for (auto& topSpGroup : topSpGroups) {
      // Find the first top space point that is within the deltaRMin of the
      // first middle space point.
      const auto low = std::ranges::lower_bound(
          topSpGroup, firstMiddleSpR + topFinder.config().deltaRMin, {},
          [&](const ConstSpacePointProxy2& sp) { return sp.zr()[1]; });
      topSpGroup = topSpGroup.subrange(low - topSpGroup.begin());
    }
  }

  for (ConstSpacePointProxy2 spM : middleSpGroup) {
    const float rM = spM.zr()[1];

    if (spacePointsSortedByRadius) {
      // check if spM is outside our radial region of interest
      if (rM < radiusRangeForMiddle.first) {
        continue;
      }
      if (rM > radiusRangeForMiddle.second) {
        // break because SPs are sorted in r
        break;
      }
    }

    createSeedsFromGroupsImpl(*m_logger, cache, bottomFinder, topFinder,
                              tripletFinder, filter, spacePoints,
                              bottomSpGroups, spM, topSpGroups, outputSeeds, evNo);
  }
}

}  // namespace Acts
