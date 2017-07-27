#ifndef SBA_DEPTH_H
#define SBA_DEPTH_H

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

#include <g2o/types/slam3d/se3quat.h>
#include <g2o/types/slam3d/edge_se3.h>
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include <g2o/types/slam3d/isometry3d_mappings.h>
#include "CSVReader.h"
#include <string>
#include <fstream>
namespace g2o {
	using namespace Eigen;
	
	/*This class enables to use slam odometry with detected markers in images
	and to use the markers for trajectory optimization.
	The markers need also to provide the depth between the camera and the marker.
	*/
	class SBADepth{
	public:
		class Options{
		public:
			Options():SCHUR_TRICK(true),ROBUST_KERNEL(false){}
			bool SCHUR_TRICK;
			bool ROBUST_KERNEL;
		};

		class SE3QuatStamped{
		public:
			typedef Eigen::Matrix<double, 3, 3, Eigen::ColMajor> InformationType;
			SE3Quat se3;
			uint64_t time;
			InformationType information; //inverse of covariance

		};
		class MarkerObservation{
		public:
			uint64_t id; //unique id. use same ID of obervation of the same physical marker
			Vector3d xyz_estimate_in_cam_frame; //marker 3d coordinates in the camera frame
			uint64_t pose_time;
			Vector3d uvd; //<!- Projection on the undistorted image u,v,depth

			
		};

		Options options;

		SBADepth():overall_vertex_counter(0),optimization_initialized(false){}

		void setCamParams(float fx, float fy, float cx, float cy)
		{
			cam_params.setKcam(fx, fy, cx, cy);
			//you may also set the depth sensor offset via cam_params->setOffset(Isometry3D)
			cam_params.setId(0);
			optimizer.clearParameters();
			if (!optimizer.addParameter(&cam_params)){
				assert(false);
			}
		}

		void addPose(SE3QuatStamped& pose)
		{
			poses.push_back(pose);
			time_pose_map[pose.time]=pose;
		}

		void addMarker(MarkerObservation& m){marker_obs.push_back(m);}



		void setupOptimization()
		{
			if (options.SCHUR_TRICK){
			    g2o::BlockSolver_6_3::LinearSolverType * linearSolver;

			    linearSolver = new g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>();

				//linearSolver = new g2o::LinearSolverCSparse<g2o::BlockSolver_6_3::PoseMatrixType>();
			    g2o::BlockSolver_6_3 * solver_ptr = new g2o::BlockSolver_6_3(linearSolver);
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

			addPosesAsVertices();
			addMarkersAsVertices();
			addMarkerEdges();
			optimization_initialized=true;
		}

		void optimize(bool verbose, int num_iterations)
		{
			if (!optimization_initialized)
				throw std::runtime_error("You need to call SBA_Depth::setupOptimization before calling optimize()");

			optimizer.initializeOptimization();
			optimizer.setVerbose(verbose);
			optimizer.optimize(num_iterations);
		}

		/*@brief expects a csv file with these columns:
		#pose_time,qw,qx,qy,qz,tx,ty,tz
		*/
		void readPosesFromFile(std::string& file_path)
		{
			CSVReader reader;
		    reader.read(file_path);

		    CSVReader::LinesType lines=reader.lines;
		    for (size_t i=0; i<reader.lines.size(); i++)
		    {
		        SE3QuatStamped pose;
		        pose.se3.setRotation(Quaterniond(std::stod(reader.lines[i].elements[1]),
		        	std::stod(reader.lines[i].elements[2]),
		        	std::stod(reader.lines[i].elements[3]),
		        	std::stod(reader.lines[i].elements[4])));
		        pose.se3.setTranslation(Vector3d(std::stod(reader.lines[i].elements[5]),
		        	std::stod(reader.lines[i].elements[6]),
		        	std::stod(reader.lines[i].elements[7])));
		        pose.time=std::stol(reader.lines[i].elements[0]);
		        poses.push_back(pose);
		    }
		};

		/*@brief expects a csv file with these columns:
		#id,pose_time,u,v,depth
		*/
		void readMarkersFromFile(std::string& file_path)
		{
			CSVReader reader;
		    reader.read(file_path);

		    CSVReader::LinesType lines=reader.lines;
		    for (size_t i=0; i<reader.lines.size(); i++)
		    {
		        MarkerObservation mo;
		        mo.id=std::stol(reader.lines[i].elements[0]);
		        mo.uvd[0]=std::stod(reader.lines[i].elements[1]);
		        mo.uvd[1]=std::stod(reader.lines[i].elements[2]);
		        mo.uvd[2]=std::stod(reader.lines[i].elements[3]);
		        marker_obs.push_back(mo);
		    }
		};

		/*@brief
		Writes all the vertex-pose estimates to a csv file.
		The format: Isometry3D --> (x, y, z, qx, qy, qz, qw) for a single row.
		*/
		void exportVertexPosesToCSV(std::string dest_filepath)
		{
			std::ofstream csv_file;
			csv_file.open (dest_filepath);
			csv_file << "#time,tx,ty,tz,qx,qy,qz,qw\n"; //covariances are not considered yet., , cov(9, col-major)\n";
			IOFormat csvFmt(6, 0, ",", ";", "", "");

			for (g2o::HyperGraph::VertexIDMap::iterator v_it_poses = optimizer.vertices().begin(); v_it_poses != optimizer.vertices().end(); v_it_poses++)
			{
				
				if (dynamic_cast<g2o::VertexSE3*>(v_it_poses->second))
				{
					g2o::VertexSE3 * v_pose = dynamic_cast< g2o::VertexSE3 * > (v_it_poses->second);
					g2o::Isometry3D est_pose = v_pose->estimate();
					//Isometry3D --> (x, y, z, qx, qy, qz, qw)
					g2o::Vector7d est_pose_vec = g2o::internal::toVectorQT(est_pose);
					csv_file << poses[v_pose->id()].time << "," << est_pose_vec.format(csvFmt) << "\n";
				}
			}
			csv_file.close();
		}


	protected:
		bool optimization_initialized;
		std::map<uint64_t, SE3QuatStamped> time_pose_map; //pose time to pose
		std::map<uint64_t, g2o::VertexSE3*> vertex_pose_map; //pose time to vertex pointer
		std::map<uint64_t, g2o::VertexSBAPointXYZ*> vertex_marker_map; //MarkerObservation.id to its vertex pointer

		std::vector<SE3QuatStamped> poses;
		std::vector<MarkerObservation> marker_obs;
		g2o::ParameterCamera cam_params;
		int overall_vertex_counter;

		g2o::SparseOptimizer optimizer;
		g2o::OptimizationAlgorithmLevenberg * solver;

		/*@brief
		Returns all marker observations for the corresponding unique marker id
		*/
		std::vector<MarkerObservation> markers_by_id(uint64_t marker_id)
		{
			std::vector<MarkerObservation> res;
			for (size_t i=0; i<marker_obs.size(); i++)
			{
				if (marker_obs[i].id=marker_id)
					res.push_back(marker_obs[i]);
			}
			return res;
		}

		/*@brief
		Searches for all markers referring to the same marker_id and estimates the
		mean value of all the observations of the marker in the world frame.
		*/
		Vector3d meanMarkerPositionInWorld(uint64_t marker_id)
		{
			float n=0;
			Vector3d mean_in_w(0,0,0);
			for (size_t i=0;i<marker_obs.size(); i++)
			{
				if (marker_obs[i].id==marker_id)
				{
					mean_in_w+=time_pose_map[marker_obs[i].pose_time].se3*marker_obs[i].xyz_estimate_in_cam_frame;
					n+=1.0;
				}
			}
			if (n==0)
			{
				throw std::runtime_error("No marker observations with marker id"+std::to_string(marker_id)+" could be found.");
			}

			mean_in_w=mean_in_w/n;
			return mean_in_w;
		}

		/*
		@brief
		Returns all poses from where a marker with the given unique id has been observed.
		*/
		std::vector<SE3QuatStamped> getAllPosesForMarkerId(uint64_t marker_id)
		{
			std::vector<SE3QuatStamped> res;
			for (size_t i=0; i<marker_obs.size(); i++)
			{
				if (marker_obs[i].id==marker_id)
				{
					res.push_back(time_pose_map[marker_obs[i].pose_time]);
				}
			}
			return res;
		}

		std::vector<MarkerObservation> getAllObservationsForMarkerId(uint64_t marker_id)
		{
			std::vector<MarkerObservation> res;
			for (size_t i=0; i<marker_obs.size(); i++)
			{
				if (marker_obs[i].id==marker_id)
				{
					res.push_back(marker_obs[i]);
				}
			}
			return res;
		}

		void addPosesAsVertices()
		{
			for (size_t i=0; i<poses.size(); i++)
			{
				g2o::VertexSE3 * v_se3 = new g2o::VertexSE3();
				v_se3->setId(i);
				v_se3->setEstimate(poses[i].se3);
				overall_vertex_counter+=1;
				vertex_pose_map[poses[i].time]=v_se3;
			}
		}

		void addOdometryEdges()
		{
			for (size_t i = 1; i < poses.size(); i++)
			{

				SE3Quat pose1=poses[i-1].se3;
				SE3Quat pose2=poses[i].se3;
				SE3Quat rel_pose=pose1.inverse()*pose2;

				g2o::EdgeSE3* odometry = new g2o::EdgeSE3;
				odometry->vertices()[0] = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex_pose_map[poses[i-1].time]); //0: pose optimizer.vertex(vertex_pose_map[poses[i-1].time]);
				odometry->vertices()[1] = dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex_pose_map[poses[i].time]);

				//TODO: Compute the covariance/information of the relative pose properly!
				g2o::EdgeSE3::InformationType information_delta_pose;
				information_delta_pose.setIdentity();
				//TODO:_ What comes first? Rotation or translation???
				information_delta_pose.matrix().block(3,3,3,3)*=0.1; //translation?
				information_delta_pose.matrix().block(0,0,3,3)*=90000; //rotation?

				odometry->setMeasurement(rel_pose);
				//information_delta_pose *= 1000.0;
				odometry->setInformation(information_delta_pose);
				optimizer.addEdge(odometry);

			}
		}

		/*@brief
		Returns a list of MarkerObservations, where a physical marker id
		is present only once
		*/
		std::map<uint64_t/*id*/, MarkerObservation> getUniqueMarkers()
		{
			std::map<uint64_t, MarkerObservation> uniqe_markers;
			for (std::vector<MarkerObservation>::iterator it=marker_obs.begin();
				it!=marker_obs.end(); it++)
			{
				if (uniqe_markers.find(it->id)==uniqe_markers.end())
				{
					uniqe_markers[it->id]=*it;
				}
			}
			return uniqe_markers;

		}

		void addMarkersAsVertices()
		{
			//collect uniqe markers and create a vertex for each of them
			std::map<uint64_t/*id*/, MarkerObservation> uniqe_markers=getUniqueMarkers();
			for (std::map<uint64_t/*id*/, MarkerObservation>::iterator it=uniqe_markers.begin(); 
				it!=uniqe_markers.end();
				it++)
			{
				g2o::VertexSBAPointXYZ * v_marker = new g2o::VertexSBAPointXYZ();
				v_marker->setId(overall_vertex_counter);
				if (options.SCHUR_TRICK){
					v_marker->setMarginalized(true);
				}
				Vector3d point_w=meanMarkerPositionInWorld(it->second.id);
				v_marker->setEstimate(point_w);
				optimizer.addVertex(v_marker);

				vertex_marker_map[it->second.id]=v_marker;
				overall_vertex_counter++;
			}
		}

		void addMarkerEdges()
		{
			
			for (std::vector<MarkerObservation>::iterator it=marker_obs.begin();
				it!=marker_obs.end(); it++)
			{
				g2o::EdgeSE3PointXYZDepth * e = new g2o::EdgeSE3PointXYZDepth();
				e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex_pose_map[it->pose_time])); //0: pose 
				e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex_marker_map[it->id]));
				e->setMeasurement(it->uvd);//the marker measurement in image (u,v,depth)

				//TODO: Compute the covariances (important)
				e->information() = Matrix3d::Identity()*10;
				e->information()(2, 2) = 3;

				if (options.ROBUST_KERNEL) {
					g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
					e->setRobustKernel(rk);
				}
				e->setParameterId(0, 0); //grants access to the cam_param object
				optimizer.addEdge(e);

			}
		}



	}; //end class SBADepth
} //end namespace

#endif