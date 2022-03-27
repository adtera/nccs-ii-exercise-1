 # include <stdio.h>
#include <stdlib.h>
# include <iostream>
# include <mpi.h>
# include <cmath>

int main(int argc, char ** argv)
{
    int i,j,a,b,c,sum,squares ,master, worker1,worker2,worker3, myrank,numprocs,tag,buffer;
    int inbuffer[4];
    double answer;


    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Status status;
    // Debug
    std::cout << "numprocs = " << numprocs << std::endl;
    
  
    // Build Cartesian Grid
    std::cout << "Rank=" << myrank << std::endl;
    int ndims = 2;
    int dims[2] = {0,0};
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Dims_create(numprocs, ndims, dims);
    std::cout << "dims=(" << dims[0] << "," << dims[1] <<  ")" << std::endl;
    MPI_Comm GRID_COMM;
    int bcs[2] = {0,0};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &GRID_COMM);
    int coords[2] = {};
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Cart_coords(GRID_COMM, myrank, ndims, coords);
    std::cout << "coords=(" << coords[0] << "," << coords[1] << ")" << std::endl;
    
    
    
    // create subdomain size dynamically
    
    int resolution = 20;
    double stepsize = 1/((float)resolution-1);
    double master_domain[resolution][resolution] {};
    for (i=0; i < resolution; i++) {
        master_domain[i][resolution - 1] = sin(2*M_PI*i*stepsize) * sinh(2*M_PI);
        }
    
    for (i=0;i<resolution;i++){
            for (j=0;j<resolution;j++){
                std::cout << myrank << "master value (" << i << " , " << j << ") = " << master_domain[i][j] << std::endl;
            }
        }
    
    
    int width, width_x, width_y, last_width, last_width_x, last_width_y;
    int modulo, modulo_x, modulo_y;
//    width = (resolution - 2)/numprocs;
    width_x = (resolution - 2)/dims[0];
    width_y = (resolution - 2)/dims[1];
    
//    modulo = resolution % numprocs;
    modulo_x = resolution % dims[0];
    modulo_y = resolution % dims[1];
    
//    last_width = width + modulo;
    last_width_x = width_x + modulo_x;
    last_width_y = width_y + modulo_y;
    
//    int subdomain[resolution][width+2];
    double subdomain[width_x+2][width_y+2];
//    int last_subdomain[resolution][last_width+2];
    double x_last_subdomain[last_width_x+2][width_y+2];
    double y_last_subdomain[width_x+2][last_width_y+2];
    double xy_last_subdomain[last_width_x+2][last_width_y+2];
    

    // subdomain
    
    for (int coord_x = 0; coord_x < dims[0]; coord_x++){
        for (int coord_y = 0; coord_y < dims[1]; coord_y++){
        
    //
            if (coord_x + 1 != dims[0] && coord_y + 1 != dims[1]) {
                for (i=0;i<width_x+2;i++){
                    for (j=0;j<width_y+2;j++){
                        subdomain[i][j] = master_domain[width_x*coord_x + i][width_y*coord_y + j];
                        
                        std::cout << myrank << "sub value (" << i << " , " << j << ") = " << subdomain[i][j] << std::endl;
                        }
                }
            // x_last
            } else if (coord_x + 1 == dims[0] && (coord_y + 1) != dims[1]) {
                for (i=0;i<last_width_x+2;i++){
                    for (j=0;j<width_y+2;j++){
                        x_last_subdomain[i][j] = master_domain[width_x*coord_x + i][width_y*coord_y + j];
                        
                        std::cout << myrank << "sub value (" << i << " , " << j << ") = " << x_last_subdomain[i][j] << std::endl;
                        }
                }
            //y_last
            } else if (coord_x + 1 != dims[0] && coord_y + 1 == dims[1]) {
                for (i=0;i<width_x+2;i++){
                    for (j=0;j<last_width_y+2;j++){
                        y_last_subdomain[i][j] = master_domain[width_x*coord_x + i][width_y*coord_y + j];
                        
                        std::cout << myrank << "sub value (" << i << " , " << j << ") = " << y_last_subdomain[i][j] << std::endl;
                        }
                }
            // xy_last
            } else  if (coord_x + 1 == dims[0] && coord_y + 1 == dims[1]) {
        //    } else {
                for (i=0;i<last_width_x+2;i++){
                    for (j=0;j<last_width_y+2;j++){
                        xy_last_subdomain[i][j] = master_domain[width_x*coord_x + i][width_y*coord_y + j];
                        
                        std::cout << myrank << "sub value (" << i << " , " << j << ") = " << xy_last_subdomain[i][j] << std::endl;
                        }
                }
            }
            MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank)
            MPI_SEND
            
        }
    }
    //
    // for loop start
    // MPI BARRIER
    // MPI SOLVE
    // MPI_SEND
    // MPI_RECV
    // end
    
    
/*
    if (myrank+1 != numprocs) {
        for (i=0;i<resolution;i++){
            for (j=0;j<width+2;j++){
                subdomain[i][j] = master_domain[i][width*myrank+j];
                }
        }
    } else {
        for (i=0;i<resolution;i++){
            for (j=0;j<last_width+2;j++){
                last_subdomain[i][j] = master_domain[i][width*myrank+j];
                }
        }
    }

    
    
    for (i=0;i<resolution;i++){
            for (j=0;j<width+2;j++){
                std::cout << myrank << "sub value (" << i << " , " << j << ") = " << subdomain[i][j] << std::endl;
            }
        }
 */
        
        
        
    // create subdomain
    int subdomains[numprocs][5][5] {};
    for (int n=0; n<numprocs;n++){
        for (i=1;i<4;i++){
            for (j=1;j<4;j++){
                subdomains[n][i][j] = 1;
            }
        }
    }
    // Initialised subdomain with ghost layer on boundary, subdomain is now matrix of 1s surrounded by 0s, 
    //Here subdomain of process 1 is shown 
    for (i=0;i<5;i++){
        for (j=0;j<5;j++){
            std::cout << "value (" << i << " , " << j << ") = " << subdomains[myrank][i][j] << std::endl;
        }
    }







    // Create Directions framework
    enum DIRECTIONS {UP,DOWN,RIGHT,LEFT};
    char* neighbours_names[4]{"up","down", "right", "left"};
    int neighbours_ranks[4];

    // Let every rank know of his neigbhour processes and write them to neighbours_ranks
    std::cout << "I am process number: " << myrank << " and I have coords(" << coords[0] << "," << coords[1] << ") " << std::endl;
    MPI_Cart_shift(GRID_COMM, 1,1, &neighbours_ranks[DOWN], &neighbours_ranks[UP]);
    MPI_Cart_shift(GRID_COMM, 0,1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);

    // Print out which process is neughbor with whom
    for(int i = 0; i < 4; i++)
    {
        if(neighbours_ranks[i] == MPI_PROC_NULL)
            printf("[MPI process %d] I have no %s neighbour.\n", myrank, neighbours_names[i]);
        else
            printf("[MPI process %d] I have a %s neighbour: process %d.\n", myrank, neighbours_names[i], neighbours_ranks[i]);
    }
    
    // Pairwise: Send data to DOWN neighbour - Receive data from UP neighbour 
    // Retrieve data, that needs to be sent DOWN, which is the penultimate lowest line of a matrix
    // Matrix is 3x3, but with ghost layers (and BC) it is 5x5
    int row_to_send[3];
    for (i=0;i<3;i++){
        row_to_send[i] = subdomains[myrank][3][i+1];
        std::cout << row_to_send[i] << std::endl;
    }
    

    MPI_Send(&row_to_send,1,MPI_INT,neighbours_ranks[DOWN],1,GRID_COMM);
    std::cout << "Process number: " << myrank << " I sent to DOWN Process number " << neighbours_ranks[DOWN]  << std::endl;
    MPI_Recv(&inbuffer,1,MPI_INT,neighbours_ranks[UP],1,GRID_COMM,&status);
    std::cout << "Process number: " << myrank << " I received from UP Process number " << neighbours_ranks[UP]  << std::endl;
    
    // Write received from UP Neigbour in own upper ghost layer, that is the most upper layer
    for (i=1;i<4;i++){
        subdomains[myrank][0][i] = row_to_send[i-1];
    }


    // Print out new domain
    for (i=0;i<5;i++){
        for (j=0;j<5;j++){
            std::cout << "value (" << i << " , " << j << ") = " << subdomains[myrank][i][j] << std::endl;
        }
    }

return 0;
}

/*
// Function for sending and receiving in arbitrary directions
void sendReceiveGHOST(DIRECTIONS direction_send_to, DIRECTIONS direction_recv_from){
    int row_to_send[3];
    for (i=0;i<3;i++){
        // Specify which row has to be sent, this is different for every direction
        if direction_send_to == DOWN{
        row_to_send[i] = subdomains[myrank][3][i+1];
        std::cout << row_to_send[i] << std::endl; 
        }

        if direction_send_to == UP{
        row_to_send[i] = subdomains[myrank][1][i+1];
        std::cout << row_to_send[i] << std::endl; 
        }

        if direction_send_to == LEFT{
        row_to_send[i] = subdomains[myrank][i+1][1];
        std::cout << row_to_send[i] << std::endl; 
        }

        if direction_send_to == RIGHT{
        row_to_send[i] = subdomains[myrank][i+1][1];
        std::cout << row_to_send[i] << std::endl; 
        }

        }
    MPI_Send(&row_to_send,1,MPI_INT,neighbours_ranks[direction_send_to],1,GRID_COMM);
    std::cout << "Process number: " << myrank << " I sent to DOWN Process number " << neighbours_ranks[DOWN]  << std::endl;
    MPI_Recv(&inbuffer,1,MPI_INT,neighbours_ranks[direction_recv,from],1,GRID_COMM,&status);
    std::cout << "Process number: " << myrank << " I received from UP Process number " << neighbours_ranks[UP]  << std::endl;
}
*/

