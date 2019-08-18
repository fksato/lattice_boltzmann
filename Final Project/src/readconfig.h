
#ifndef CONFIG_H
#define CONFIG_H

namespace LBM {
#ifdef _SEQ_
	void readConfig( const char* configFile, char* geometry
						, int *x_cells, int *y_cells, int *z_cells
						, int *time_steps, double *tau, int *vis_t 
					);
	void buildFlag( const char* geometry, int *FLAG, int *size );

#elif _MPI_
	void readConfig( const char* configFile, char* geometry
						, int *x_cells, int *y_cells, int *z_cells
						, int *time_steps, double *tau, int *vis_t
						, int *x_proc, int *y_proc, int *z_proc );

	void buildFlag( const char* geometry, int *FLAG, int *size, int x_proc, int y_proc, int z_proc, int *scatter_size );
#endif	
}

#endif