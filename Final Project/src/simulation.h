
#ifndef SIM_H
#define SIM_H

#include "LB_consts.h"
#include <memory>

#ifdef _MPI_
#include <mpi.h>
#include "helper.h"
#endif

namespace LBM {
	
	#define FLUID 		1
	#define NOSLIP 		2
	#define MOVINGWALL 	3
	#define FREESLIP    4
	#define INFLOW      5
	#define OUTFLOW     6
	#define OBS			8

	#define NEIGHBORS 		19
	#define FACES 			6
	#define SHAREDVEL 		5

	class SimDomain {

		public:
			SimDomain( int *size, double inTau, int *flagField );
			
			void stream();
			void collide();
			void swapFields();
			void treatBoundaries();
#ifdef _MPI_
			void exchangePressureGhostCells( int axis, int sendtag, int recvtag, MPI_Comm comm, MPI_Status status );

			void visualize( int t, int rank, int widthRank, int depthRank, int heightRank );
#elif _SEQ_
			void visualize( int t );
#endif
			

		protected:
			int xlen, ylen, zlen, xylen;
			int size[DIM] 	= {0};
			std::unique_ptr<double[]> collisionField;
			std::unique_ptr<double[]> streamField;
			int *FLAG;
			// VISUALIZATION:
			std::unique_ptr<double[]> visDensity;
			std::unique_ptr<double[]> visVelocity;
			double wallVelocity[DIM] = { 1.0, 0.0, 0.0 };
			double tau;
			// Outflow/Inflow:
			double p_ref_out = 1.0;
			double p_ref_in = 1.005;
			double inflowVelocity[3] = { 0.0, 0.05, 0.0 };

			void computeDensity( int cellidx, double *density );
			void computeVelocity( int cellidx, double *density, double *velocity );
			void computeFEQ( double *velocity, double density, double *f_eq );
			void computePostCollisionDistributions( int cellidx, double *f_eq );
			
			// Boundary Condition Functions:
			void noSlip( int cellidx, int *neighbors );
			void freeSlip( int cellidx, int x, int y, int z );
			void inflow( int cellidx, int *neighbors  );
			void outflow( int cellidx, int *neighbors );
			void movingWall( int cellidx, int *neighbors );

			// Helper function:
			void getNeighbors( int x, int y, int z, int *neighbors );
	};
}

#endif