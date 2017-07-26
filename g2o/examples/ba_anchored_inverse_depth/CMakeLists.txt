ADD_EXECUTABLE(ba_anchored_inverse_depth_demo
  ba_anchored_inverse_depth_demo.cpp
)

#INCLUDE_DIRECTORIES(${CHOLMOD_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CSPARSE_INCLUDE_DIR})

TARGET_LINK_LIBRARIES(ba_anchored_inverse_depth_demo core types_sba solver_csparse)


