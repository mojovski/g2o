ADD_EXECUTABLE(sba_depth_expmap
  sba_depth_expmap.cpp
)
INCLUDE_DIRECTORIES(${CSPARSE_INCLUDE_DIR})

TARGET_LINK_LIBRARIES(sba_depth_expmap core types_sba types_slam3d solver_csparse)


