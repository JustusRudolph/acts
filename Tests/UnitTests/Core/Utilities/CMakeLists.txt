add_unittest(AnnealingUtility AnnealingUtilityTests.cpp)
add_unittest(Axes AxesTests.cpp)
add_unittest(BFieldMapUtils BFieldMapUtilsTests.cpp)
add_unittest(BinAdjustment BinAdjustmentTests.cpp)
add_unittest(BinAdjustmentVolume BinAdjustmentVolumeTests.cpp)
add_unittest(BinningData BinningDataTests.cpp)
add_unittest(BinUtility BinUtilityTests.cpp)
add_unittest(BoundingBox BoundingBoxTest.cpp)
add_unittest(Extendable ExtendableTests.cpp)
add_unittest(FiniteStateMachine FiniteStateMachineTests.cpp)
add_unittest(Frustum FrustumTest.cpp)
add_unittest(Grid GridTests.cpp)
add_unittest(Helpers HelpersTests.cpp)
add_unittest(Interpolation InterpolationTests.cpp)
add_unittest(Intersection IntersectionTests.cpp)
add_unittest(KDTree KDTreeTests.cpp)
add_unittest(Logger LoggerTests.cpp)
add_unittest(MaterialMapUtils MaterialMapUtilsTests.cpp)
add_unittest(MPL MPLTests.cpp)
add_unittest(MultiIndex MultiIndexTests.cpp)
add_unittest(Periodic PeriodicTests.cpp)
add_unittest(Range1D Range1DTests.cpp)
add_unittest(RangeXD RangeXDTests.cpp)
add_unittest(Ray RayTest.cpp)
add_unittest(RealQuadraticEquation RealQuadraticEquationTests.cpp)
add_unittest(Result ResultTests.cpp)
add_unittest(Subspace SubspaceTests.cpp)
add_unittest(TypeTraits TypeTraitsTest.cpp)
add_unittest(UnitVectors UnitVectorsTests.cpp)
add_unittest(Delegate DelegateTests.cpp)
if (ACTS_BUILD_CUDA_FEATURES)  
  add_unittest(Cuda CudaTests.cu)
  add_unittest(CudaMostSimplified CudaMostSimplifiedTests.cu)
endif()
