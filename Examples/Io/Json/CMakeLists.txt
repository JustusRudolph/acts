add_library(
  ActsExamplesIoJson SHARED
  src/JsonGeometryList.cpp
  src/JsonMaterialWriter.cpp
  src/JsonSurfacesWriter.cpp
  src/JsonOptionsWriter.cpp
  src/JsonDigitizationConfig.cpp)
target_include_directories(
  ActsExamplesIoJson
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoJson
  PUBLIC ActsCore ActsPluginJson ActsExamplesDigitization ActsExamplesFramework ActsExamplesMaterialMapping Boost::program_options)

install(
  TARGETS ActsExamplesIoJson
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
