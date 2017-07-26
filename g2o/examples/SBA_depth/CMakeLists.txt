ADD_EXECUTABLE(sba_depth
  sba_depth.cpp
)
INCLUDE_DIRECTORIES(${CSPARSE_INCLUDE_DIR})

TARGET_LINK_LIBRARIES(sba_depth core types_sba types_slam3d solver_csparse)


