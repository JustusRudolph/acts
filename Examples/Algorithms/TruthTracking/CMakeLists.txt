add_library(
  ActsExamplesTruthTracking SHARED
  ActsExamples/TruthTracking/ParticleSelector.cpp
  ActsExamples/TruthTracking/ParticleSmearing.cpp
  ActsExamples/TruthTracking/ParticleSmearingOptions.cpp
  ActsExamples/TruthTracking/TrackSelector.cpp
  ActsExamples/TruthTracking/TruthSeedSelector.cpp
  ActsExamples/TruthTracking/TruthTrackFinder.cpp
  ActsExamples/TruthTracking/TruthVertexFinder.cpp
)
target_include_directories(
  ActsExamplesTruthTracking
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(
  ActsExamplesTruthTracking
  PUBLIC ActsCore ActsExamplesFramework Boost::program_options)

install(
  TARGETS ActsExamplesTruthTracking
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
