
#include "readconfig.h"
#include "LB_consts.h"
#include "helper.h"
#include <cstdio>

#ifdef _MPI_
#include <memory>
#endif

#ifdef _SEQ_

void LBM::readConfig( const char* configFile, char* geometry
						, int *x_cells, int *y_cells, int *z_cells
						, int *time_steps, double *tau, int *vis_t ) 
{
	READ_STRING( configFile, geometry );
	//
	READ_INT( configFile, *x_cells );
	READ_INT( configFile, *y_cells );
	READ_INT( configFile, *z_cells );
	//
	READ_INT( configFile, *time_steps );
	//
	READ_DOUBLE( configFile, *tau );
	//
	READ_INT( configFile, *vis_t );
}

void LBM::buildFlag( const char* filename, int *FLAG, int *size ) {

	FILE *input = NULL;
	char line[1024];
	
	int xlen = ( size[0] + 2 );
	int ylen = ( size[1] + 2 );

	if ( ( input=fopen( filename,"rb" ) ) ==0 )
    {
       char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", filename );
	   ERROR( szBuff );
    }

    /* skip the comments */
    do
    fgets( line, sizeof line, input );
    while( *line=='#' );
    
    int x, y, z, flag, idx;
    while( fscanf( input, "%i %i %i %i\n", &x, &y, &z, &flag ) == 4 ) {
    	idx = z * xlen * ylen + y * xlen + x;
    	FLAG[ idx ] = flag;
    }
    /* close file */
    fclose(input);
}

#elif _MPI_

void LBM::readConfig( const char* configFile, char* geometry
						, int *x_cells, int *y_cells, int *z_cells
						, int *time_steps, double *tau, int *vis_t
						, int *x_proc, int *y_proc, int *z_proc ) 
{
	READ_STRING( configFile, geometry );
	//
	READ_INT( configFile, *x_proc );
	READ_INT( configFile, *y_proc );
	READ_INT( configFile, *z_proc );
	//
	READ_INT( configFile, *x_cells );
	READ_INT( configFile, *y_cells );
	READ_INT( configFile, *z_cells );
	//
	READ_INT( configFile, *time_steps );
	//
	READ_DOUBLE( configFile, *tau );
	//
	READ_INT( configFile, *vis_t );
}

void LBM::buildFlag( const char* filename, int *FLAG, int *size
					, int x_proc, int y_proc, int z_proc, int *scatter_size ) {

	FILE *input = NULL;
	char line[1024];
	
	int xlen = ( size[0] + 2 );
	int ylen = ( size[1] + 2 );
	int zlen = ( size[2] + 2 );
	std::unique_ptr<int[]>temp( new int[ xlen * ylen * zlen] );

	if ( ( input=fopen( filename,"rb" ) ) ==0 )
    {
       char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", filename );
	   ERROR( szBuff );
    }

    /* skip the comments */
    do
    fgets( line, sizeof line, input );
    while( *line=='#' );
    
    int x, y, z, flag, idx;
    while( fscanf( input, "%i %i %i %i\n", &x, &y, &z, &flag ) == 4 ) {
    	idx = z * xlen * ylen + y * xlen + x;
    	temp[ idx ] = flag;
    }

    int sizes[DIM];
    int global_cnt = 0;
    
    int local_size_x = size[0] / x_proc;
    int local_size_y = size[1] / y_proc;
    int local_size_z = size[2] / z_proc;

    int local_xlen;
    int local_ylen;
    int local_zlen;
    
    for( int i = 0; i < x_proc; i++ ) {
		for( int j = 0; j < y_proc; j++ ) {
			for( int k = 0; k < z_proc; k++ ) {

		    	idx = k *(x_proc*y_proc) + j * x_proc + i;

		    	sizes[0] = scatter_size[ idx + 0 ];
		    	sizes[1] = scatter_size[ idx + 1 ];
		    	sizes[2] = scatter_size[ idx + 2 ];

		    	local_xlen = sizes[0] + 2;
		    	local_ylen = sizes[1] + 2;
		    	local_zlen = sizes[2] + 2;

    			for( int height = 0; height < local_zlen; height++ ) {
    				for( int depth = 0; depth < local_ylen; depth++ ) {
						for( int width = 0; width < local_xlen; width++ ) {

    						FLAG[ global_cnt + height * ( local_xlen * local_ylen ) + ( depth * local_xlen ) + width ]
    								= temp[ ( k * local_size_z + height * xlen * ylen ) 
    									+ ( j * local_size_y + depth * xlen ) 
    									+ ( i * local_size_x + width ) ];
    					}
    				}
    			}
    			global_cnt += local_xlen * local_ylen * local_zlen;
    		}
    	}
    }
    /* close file */
    fclose(input);
}
#endif
