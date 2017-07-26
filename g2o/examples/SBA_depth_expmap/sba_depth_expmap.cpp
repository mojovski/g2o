// g2o - General Graph Optimization
// Copyright (C) 2012 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Eigen/StdVector>
#include <iostream>
#include <stdint.h>

#include <unordered_set>

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include <g2o/solvers/eigen/linear_solver_eigen.h>

#include "g2o/types/sba/types_six_dof_expmap.h"
//#include "g2o/math_groups/se3quat.h"
#include "g2o/solvers/structure_only/structure_only_solver.h"
#include <cmath>
#include <g2o/types/slam3d/edge_se3_pointxyz_depth.h>
#include <g2o/types/slam3d/edge_se3expmap_pointxyz_depth.h>

#include <g2o/types/slam3d/se3quat.h>
#include <g2o/types/slam3d/edge_se3.h>
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include <g2o/types/slam3d/isometry3d_mappings.h>


using namespace Eigen;
using namespace std;
#ifndef PI
	#define PI 3.14159265
#endif

class Sample{
public:
  static int uniform(int from, int to);
  static double uniform();
  static double gaussian(double sigma);
};

static double uniform_rand(double lowerBndr, double upperBndr){
  return lowerBndr + ((double) std::rand() / (RAND_MAX + 1.0)) * (upperBndr - lowerBndr);
}

static double gauss_rand(double mean, double sigma){
  double x, y, r2;
  do {
    x = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
    y = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0.0);
  return mean + sigma * y * std::sqrt(-2.0 * log(r2) / r2);
}

int Sample::uniform(int from, int to){
  return static_cast<int>(uniform_rand(from, to));
}

double Sample::uniform(){
  return uniform_rand(0., 1.);
}

double Sample::gaussian(double sigma){
  return gauss_rand(0., sigma);
}

Vector2d project2d(const Vector3d& v){
  Vector2d res;
  res(0) = v(0)/v(2);
  res(1) = v(1)/v(2);
  return res;
}

Vector3d unproject2d(const Vector2d& v){
  Vector3d res;
  res(0) = v(0);
  res(1) = v(1);
  res(2) = 1;
  return res;
}

inline Vector3d invert_depth(const Vector3d & x){
  return unproject2d(x.head<2>())/x[2];
}

void print_edge_errors(g2o::SparseOptimizer& optimizer, int num_odometry_edges)
{
	/*
	Print all the edges errors
	*/
	int edge_id = 0;
	int n= optimizer.edges().size();
	for (g2o::HyperGraph::EdgeSet::iterator eit = optimizer.edges().begin(); eit != optimizer.edges().end(); eit++)
	{

		if (dynamic_cast<g2o::EdgeSE3*>(*eit))
		{
			g2o::EdgeSE3Expmap* e = dynamic_cast<g2o::EdgeSE3Expmap*>(*eit);
			e->computeError();
			double err = e->error().norm();
			std::cout << "Error Odometry Edge id: " << edge_id << ", error: " << err << "\n";
			edge_id += 1;
		}
		else
		{
			g2o::EdgeSE3PointXYZDepth* e = dynamic_cast<g2o::EdgeSE3PointXYZDepth*>(*eit);
			e->computeError();
			double err = e->error().norm();
			std::cout << "Error Depth Edge id: " << edge_id << ", error: " << err << "\n";
			edge_id += 1;
		}
	}
}

void print_abs_pose_errors(g2o::SparseOptimizer& optimizer, vector<g2o::SE3Quat, aligned_allocator<g2o::SE3Quat> > & true_poses);

void print_point_errors(g2o::SparseOptimizer& optimizer, vector<Vector3d>& true_markers, std::map<int, std::vector<int> >& marker_poses_links, int num_pose_vertices);

int main(int argc, const char* argv[]){
  //if (argc<2)
  //{
    cout << endl;
    cout << "Please type: " << endl;
    cout << "sba_depth [PIXEL_NOISE] [OUTLIER RATIO] [ROBUST_KERNEL] [SCHUR-TRICK]" << endl;
    cout << endl;
    cout << "PIXEL_NOISE: noise in image space (E.g.: 1)" << endl;
    cout << "OUTLIER_RATIO: probability of spuroius observation  (default: 0.0)" << endl;
    cout << "ROBUST_KERNEL: use robust kernel (0 or 1; default: 0==false)" << endl;
    cout << "SCHUR-TRICK: Use Schur-complement trick (0 or 1; default: 0==true)" << endl;
    cout << endl;
    cout << "Note, if OUTLIER_RATIO is above 0, ROBUST_KERNEL should be set to 1==true." << endl;
    cout << endl;
    //exit(0);
 // }

  double PIXEL_NOISE = 0.0; //atof(argv[1]);

  double OUTLIER_RATIO = 0.0;
  if (argc>2){
    OUTLIER_RATIO = atof(argv[2]);
  }

  bool ROBUST_KERNEL = false;
  if (argc>3){
    ROBUST_KERNEL = atoi(argv[3]) != 0;
  }

  bool SCHUR_TRICK = true;
  if (argc>4){
    SCHUR_TRICK = atoi(argv[4]) != 0;
  }

  cout << "PIXEL_NOISE: " <<  PIXEL_NOISE << endl;
  cout << "OUTLIER_RATIO: " << OUTLIER_RATIO<<  endl;
  cout << "ROBUST_KERNEL: " << ROBUST_KERNEL << endl;
  cout << "SCHUR-TRICK: " << SCHUR_TRICK << endl;

  g2o::SparseOptimizer optimizer;
  //SparseOptimizer optimizer;
  //g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType> *linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

  optimizer.setVerbose(true);
  g2o::OptimizationAlgorithmLevenberg * solver;

  if (SCHUR_TRICK){
    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

	//linearSolver = new g2o::LinearSolverCSparse<g2o::BlockSolver_6_3::PoseMatrixType>();
	cerr << "Using CSPARSE" << endl;

    g2o::BlockSolver_6_3 * solver_ptr
        = new g2o::BlockSolver_6_3(linearSolver);
    solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
  } else {
    g2o::BlockSolverX::LinearSolverType * linearSolver;
    linearSolver
        = new g2o::LinearSolverEigen<g2o
        ::BlockSolverX::PoseMatrixType>();
    g2o::BlockSolverX * solver_ptr
        = new g2o::BlockSolverX(linearSolver);
    solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
  }

  optimizer.setAlgorithm(solver);

  //====================================================
  //  Setup the scene
  vector<Vector3d> true_markers;
  Vector3d m1(0,0,1);
  Vector3d m2(1,0,1);
  Vector3d m3(1,1,1);
  true_markers.push_back(m1);
  true_markers.push_back(m2);
  true_markers.push_back(m3);

  /* 
  Poses:
  The camera moves on circle on the z plane
  Create the Gt poses and add them as vertices VertexSE3Expmap
  */
  vector<g2o::SE3Quat, aligned_allocator<g2o::SE3Quat> > true_poses;
  vector<g2o::SE3Quat, aligned_allocator<g2o::SE3Quat> > noisy_poses;

  int n=20;
  double angle_step=360.0/(double)n/180.0*PI;
  int overall_vertex_counter=0;
  for (int i=0; i<n; i++)
  {
    Vector3d t(sin(angle_step*i), //+Sample::gaussian(0.1), 
      1.0-cos(angle_step*i), //+Sample::gaussian(0.1), 
		-1); // +Sample::gaussian(0.1));

    Eigen::Quaterniond q(1.0, 0, 0, 0);
    g2o::SE3Quat pose(q, t);
    true_poses.push_back(pose);

	g2o::SE3Quat noisy_pose(pose);
	if (i>0)
		noisy_pose.setTranslation(pose.translation() + g2o::Vector3D(Sample::gaussian(0.1), Sample::gaussian(0.1), Sample::gaussian(0.1)));

	noisy_poses.push_back(noisy_pose);


	g2o::SE3Quat est_pose = noisy_pose;
	//est_pose.setTranslation(est_pose.translation() + g2o::Vector3D(Sample::gaussian(0.1), Sample::gaussian(0.1), Sample::gaussian(0.1)));

    //create a vertex representing the pose
    g2o::VertexSE3Expmap * v_se3
        = new g2o::VertexSE3Expmap();
    v_se3->setId(overall_vertex_counter);
    v_se3->setEstimate(est_pose);
	if (i == 0)
		v_se3->setFixed(true);
    optimizer.addVertex(v_se3);
    overall_vertex_counter++;

  }
  int num_pose_vertices=n;


  //assign the links between the poses and the markers?
  std::map<int,std::vector<int> > markers_poses_links;
  markers_poses_links[0].push_back(0); //0.th marker is visible in poses 0,1,2
  markers_poses_links[0].push_back(1);
  markers_poses_links[0].push_back(2);

  markers_poses_links[1].push_back(4);
  markers_poses_links[1].push_back(7);
  markers_poses_links[1].push_back(10);

  markers_poses_links[2].push_back(13);
  markers_poses_links[2].push_back(14);
  markers_poses_links[2].push_back(n-1);
  

  //-------- setup the pose covariance (information matrix)

  g2o::EdgeSE3Expmap::InformationType information_pose;
  information_pose.setIdentity();
  information_pose=information_pose;
  cout << "information of each pose: \n" << information_pose << endl;


  //setup the camera properties
  double focal_length= 100.;
  Vector2d principal_point(320., 240.);

  g2o::ParameterCamera * cam_params = new g2o::ParameterCamera();
  cam_params->setKcam(focal_length, focal_length, principal_point[0], principal_point[1]);
  //you may also set the depth sensor offset via cam_params->setOffset(Isometry3D)
  cam_params->setId(0);

  //add the parameter to the optimizer, so the values can be adjusted??
  if (!optimizer.addParameter(cam_params)){
    assert(false);
  }

  /*----------------------------------
  Add the odometry edges
  */
  //noisy_poses.push_back(true_poses[0]);
  for (int i = 1; i < num_pose_vertices; i++)
  {
	  g2o::SE3Quat pose1 = noisy_poses[i - 1];
	  g2o::SE3Quat pose2 = noisy_poses[i];
	  g2o::SE3Quat rel_pose = pose1.inverse()*pose2;
	  //add noise to the odometry
	  //rel_pose.setTranslation(rel_pose.translation()+ g2o::Vector3D(Sample::gaussian(0.1), Sample::gaussian(0.1), Sample::gaussian(0.1)));
	  //noisy_poses.push_back(*(noisy_poses.end() - 1)*rel_pose);


	  g2o::EdgeSE3Expmap* odometry = new g2o::EdgeSE3Expmap;
	  odometry->vertices()[0] = optimizer.vertex(i - 1);
	  odometry->vertices()[1] = optimizer.vertex(i);

	  int id1 = optimizer.vertex(i - 1)->id();
	  int id2 = optimizer.vertex(i)->id();
	  if (id1 != (i - 1) || (id2 != i))
	  {
		  std::cerr << "Error: The Vertex IDS do not match!";
		  throw std::runtime_error("The vertex ids do not match when addeing odometry");
	  }

	  
	  odometry->setMeasurement(rel_pose);
	  g2o::EdgeSE3Expmap::InformationType rel_pose_info_mat;

	  rel_pose_info_mat.setIdentity();
	  //rel_pose_info_mat *= 1000.0;
	  odometry->setInformation(rel_pose_info_mat);
	  optimizer.addEdge(odometry);

  }
  int num_odometry_edges = optimizer.edges().size();

  


  /*---------------------------
  Add the markers as vertices and the depth edges
  */
  for (size_t i=0; i<true_markers.size(); ++i)
  {
    g2o::VertexSBAPointXYZ * v_marker = new g2o::VertexSBAPointXYZ();
    v_marker->setId(overall_vertex_counter);
    if (SCHUR_TRICK){
      v_marker->setMarginalized(true);
    }
    

    //Add noise to the 3d point estimation of the marker
	Vector3d point_w = true_markers.at(i);
	/*
                + Vector3d(Sample::gaussian(0.1),
                           Sample::gaussian(0.1),
                           Sample::gaussian(0.1)); 
						   */
	//Set
	int pose_idx0 = markers_poses_links[i][0];
	g2o::SE3Quat T_cam_in_world = noisy_poses[pose_idx0];
	g2o::SE3Quat T_world_in_cam = T_cam_in_world.inverse();
	
	Vector3d init_estimate = point_w + 
		Vector3d(Sample::gaussian(0.05),
			Sample::gaussian(0.05),
			Sample::gaussian(0.12));

	v_marker->setId(num_pose_vertices+i);
    v_marker->setEstimate(init_estimate);
    optimizer.addVertex(v_marker);

    // Increment the vertex counter
    overall_vertex_counter++;

    // Add the edges for all poses connected to this marker.
    for (size_t j=0; j<markers_poses_links[i].size(); j++)
    {
		Vector3d measurement_marker_w = true_markers.at(i);
		/*
			+ Vector3d(Sample::gaussian(0.15),
				Sample::gaussian(0.15),
				Sample::gaussian(0.15));
				*/
      int pose_idx=markers_poses_links[i][j];
      g2o::SE3Quat T_cam_in_world = true_poses[pose_idx];
      g2o::SE3Quat T_world_in_cam=T_cam_in_world.inverse();

      //Vector3d marker_in_cam_frame = T_world_in_cam*measurement_marker_w;
	  //project the marker and use the result as measurement estimation
	  Eigen::Transform<double, 3, Eigen::ColMajor> w2i;
	  w2i.matrix().setIdentity();
	  w2i.matrix().topLeftCorner<3, 4>() = cam_params->Kcam() * T_world_in_cam.to_homogeneous_matrix().topLeftCorner<3, 4>();
	  Vector3d marker_in_image_uvd = w2i*measurement_marker_w;
	  marker_in_image_uvd.head<2>() = marker_in_image_uvd.head<2>() / marker_in_image_uvd(2);

	 /* marker_in_image_uvd(0) += Sample::gaussian(0.3);
	  marker_in_image_uvd(1) += Sample::gaussian(0.3);
	  marker_in_image_uvd(2) += Sample::gaussian(0.3);
	  */

	  //Vector3d marker_in_image=cam_params->
      //marker_in_cam_frame);
      //check if the point is visible.
      //std::cout << "The marker " << i << " in pose " << pose_idx << " is:\n" << marker_in_cam_frame;

      //Add also a corresponding observation

      g2o::EdgeSE3PointXYZDepth * e = new g2o::EdgeSE3PointXYZDepth();
      //e->resize(2);
      e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertices().find(pose_idx)->second)); //0: pose 
      e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(v_marker)); 

      e->setMeasurement(marker_in_image_uvd);//the marker measurement in image (u,v,depth)
      e->information() = Matrix3d::Identity()*100;
	  e->information()(2, 2) = 1000;

      if (ROBUST_KERNEL) {
        g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
        e->setRobustKernel(rk);
      }
      e->setParameterId(0, 0); //grants access to the cam_param object??
      optimizer.addEdge(e);

    }
  }


  


  /*
  Compute the error before optimization
  */
  print_point_errors(optimizer, true_markers, markers_poses_links, num_pose_vertices);
  print_edge_errors(optimizer, num_odometry_edges);
  print_abs_pose_errors(optimizer, true_poses);


  /*
  Optimize
  */

  optimizer.initializeOptimization();
  optimizer.setVerbose(true);
  cout << endl;

  optimizer.optimize(30);


  /*Compute the error after optimization
  */
  print_point_errors(optimizer, true_markers, markers_poses_links, num_pose_vertices);
  print_edge_errors(optimizer, num_odometry_edges);
  print_abs_pose_errors(optimizer, true_poses);

  std::cout << "done\n";

}


void print_abs_pose_errors(g2o::SparseOptimizer& optimizer, vector<g2o::SE3Quat, aligned_allocator<g2o::SE3Quat> > & true_poses)
{

	std::cout << "============ Errors between GT poses and the vertex estimates================\n";
	double errors = 0;
	int num = 0;
	for (g2o::HyperGraph::VertexIDMap::iterator v_it_poses = optimizer.vertices().begin(); v_it_poses != optimizer.vertices().end(); v_it_poses++)
	{
		
		if (dynamic_cast<g2o::VertexSE3Expmap*>(v_it_poses->second))
		{
			g2o::VertexSE3Expmap * v_pose = dynamic_cast< g2o::VertexSE3Expmap * > (v_it_poses->second);
			g2o::Isometry3D est_pose = v_pose->estimate();
			g2o::SE3Quat true_pose = true_poses[v_pose->id()];

			g2o::Isometry3D delta = est_pose.inverse()*true_pose;
			Eigen::Vector3d dt=delta.translation();
			g2o::Vector6d _error = g2o::internal::toVectorMQT(delta);

			double err = _error.norm();
			std::cout << "Error Pose Vertex, id: " << v_pose->id() << ", error: " << err << ", dt: " << dt.norm() << "\n";
			errors += err;
			num += 1;
		}
		
	}
	double mean_error = errors / (double)num;
	std::cout << "Mean pose error: " << mean_error << "\n";

}

void print_point_errors(g2o::SparseOptimizer& optimizer, vector<Vector3d>& true_markers, std::map<int, std::vector<int> >& markers_poses_links, int num_pose_vertices)
{
	std::cout << "\n======================Marker Errors in World Frame============================\n";
	double sum_diff2 = 0;
	int m_num = 0;
	for (size_t i = 0; i < true_markers.size(); ++i)
	{
		g2o::VertexSBAPointXYZ * v_marker = dynamic_cast<g2o::VertexSBAPointXYZ *> (optimizer.vertices()[num_pose_vertices + i]);
		Vector3d diff = v_marker->estimate() - true_markers[i];
		std::cout << "Error Marker " << i << ": " << diff.norm() << "\n";

		sum_diff2 += diff.norm();
		m_num += 1;

	}
	cout << "Mean Marker 3d error: " << sum_diff2 / m_num << endl;
}