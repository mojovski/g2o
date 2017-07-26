# Description

This example demonstrates how to implement single pixel depth measurements
into the SBA g2o framework.
The scene considers a case, where a camera (mono) is flying around
and a 6dof-trajectory has been precomputed.
(next version 7dof??)

This builds a graph network with all the inter-pose (relative) odometries
and adds marker observations, which have been acquired from
scaled markers. This are modelled as a Edge between the vertex
VertexSBAPointXYZ (unknoen 3D point)
and the pose
VertexSE3Expmap or VertexSE3 (from slam3d).

The edge, the observation between VertexSE3 and VertexSBAPointXYZ is defined here as a custom class.


