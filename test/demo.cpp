#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define MAX_ROWS 1000 // upper limits #define MAX_COLS 1000 // upper limits #define MIN(a, b) ((a) > (b) ? (b) : (a)) #define DONE MAX_ROWS+1
int main(int argc, char **argv) {
    double** A;
    double* b, * c, * buffer;
    double ans;
    int myid, master, numprocs, i, j, numsent, sender, done, anstype, row, rows, cols;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    A = (double **)calloc(MAX_ROWS, sizeof(double *));
    for (i=0; i<MAX_ROWS; i++)
        A[i] = (double *)calloc(MAX_COLS, sizeof(double));
    b = (double *)calloc(MAX_COLS, sizeof(double));
    c = (double *)calloc(MAX_ROWS, sizeof(double));
    buffer = (double *)calloc(MAX_COLS, sizeof(double));
    master = 0;
    rows = 2; // for debugging ...
    cols = 2; // for debugging ...
    if (myid == master) {
        /* master code: next slides */
        /* Initialize A and b (arbitrary) */
        A[0][0] = 1.0; A[0][1] = 2.0; A[1][0] = 1.0; A[1][1] = 2.0; b[0] = 1.0; b[1] = 1.0;
        numsent = 0;
        /* Send b to each worker process */
        MPI_Bcast(b, cols, MPI_DOUBLE, master, MPI_COMM_WORLD); /* Send a row to each worker process; tag with row number */
        for (i = 0; i < MIN(numprocs - 1, rows); i++) {
            MPI_Send(&A[i][0], cols, MPI_DOUBLE, i+1, i, MPI_COMM_WORLD);
            numsent++;
            for (i = 0; i < rows; i++) {
                MPI_Recv(&ans, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                sender = status.MPI_SOURCE; /* row is tag value */
                anstype = status.MPI_TAG; c[anstype] = ans;
                /* send another row */
                if (numsent < rows) {
                    MPI_Send(&A[numsent][0], cols, MPI_DOUBLE, sender, numsent, MPI_COMM_WORLD);
                    numsent++;
                } else {
                    /* Tell sender that there is no more work */
                    MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, DONE, MPI_COMM_WORLD);
                }
            }
        }
    } else { 
        /* worker code: next slides */
        MPI_Bcast(b, cols, MPI_DOUBLE, master, MPI_COMM_WORLD); /* Skip if more processes than work */
        done = myid > rows;
        while (!done) {
            MPI_Recv(buffer, cols, MPI_DOUBLE, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            done = status.MPI_TAG == DONE; if (!done) {
                row = status.MPI_TAG; ans = 0.0;
                for (i = 0; i < cols; i++) {
                    ans += buffer[i] * b[i];
                }
                MPI_Send(&ans, 1, MPI_DOUBLE, master, row, MPI_COMM_WORLD);
            }
        }
    }
// result c vector finished: considering the test values, all element values should be 3.0 for (i=0; i<MAX_ROWS; i++)
    for (i=0; i<MAX_ROWS; i++)
        free(A[i]);
    free(A);
    free(b);
    free(c);
    free(buffer);
    MPI_Finalize();
    return 0;
    }