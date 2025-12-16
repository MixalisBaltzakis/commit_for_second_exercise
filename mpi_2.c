#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>


typedef struct {
    int n;          // number of nodes
    long long nnz;        // number of nonzeros (edges * 2 for undirected)
    int local_n;    // number of local process nodes
    long long local_nnz;   // number of local process nonzeros
    long long *row_ptr;   // column pointers (size n+1)
    int *col_ind;   // row indices (size nnz)
    double *values; // all ones (binary)
} CSCMatrix;

void print_csc(const CSCMatrix *A) {
    printf("Binary undirected graph: %d nodes, %d edges (nnz=%d)\n",
           A->n, A->nnz / 2, A->nnz);
    printf("\nrow_ptr = [ ");
    for (int j = 0; j <= A->n; j++)
        printf("%d ", A->row_ptr[j]);
    printf("]\n");

    printf("col_ind = [ ");
    for (int k = 0; k < A->nnz; k++)
        printf("%d ", A->col_ind[k]);
    printf("]\n\n");
}

CSCMatrix *load_problem_parallel(int my_rank, int procnum, const char *dataset_name) {
    
    int n;
    long long nnz;
    char filepath[512];

    // 1. Διάβασμα Metadata
    if (my_rank == 0) {
        sprintf(filepath, "%s/meta.txt", dataset_name);
        
        FILE *f = fopen(filepath, "r");
        if (!f) { printf("Error opening meta.txt\n"); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, 1); }
        fscanf(f, "%d\n%lld", &n, &nnz);
        fclose(f);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    // Last Process takes Remainder
   
    int base_chunk = n / procnum;   // Το κομμάτι που παίρνουν όλοι
    int remainder  = n % procnum;   // Οτι περισσευει το παιρνει ο τελευταιος
    
    int local_n = base_chunk;

    if (my_rank == procnum - 1) {
        local_n += remainder;
    }
    
    // Επειδή ΟΛΟΙ οι προηγούμενοι από εμένα πήραν σίγουρα `base_chunk`,
    int my_start_node = my_rank * base_chunk;

    CSCMatrix *A = (CSCMatrix *)malloc(sizeof(CSCMatrix));
    A->n = n;
    A->nnz = nnz;
    A->local_n = local_n;
    A->row_ptr = (long long *)malloc((local_n + 1) * sizeof(long long));

    // 3. MPI I/O για το row_ptr.bin
    sprintf(filepath, "%s/row_ptr.bin", dataset_name);
    
    MPI_File fh;
    MPI_Status status;
    
    int err = MPI_File_open(MPI_COMM_WORLD, filepath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err) { printf("Error opening %s\n", filepath); MPI_Abort(MPI_COMM_WORLD, 1); }

    MPI_Offset offset = (MPI_Offset)my_start_node * sizeof(long long);  
    
    // Διαβάζουμε local_n + 1 στοιχεία
    MPI_File_read_at(fh, offset, A->row_ptr, local_n + 1, MPI_LONG_LONG, &status);
    MPI_File_close(&fh);

    // 4. Υπολογισμός local NNZ (Edges)
    long long start = A->row_ptr[0]; 
    long long end   = A->row_ptr[local_n];
    A->local_nnz = end - start;
    A->col_ind = (int *)malloc(A->local_nnz * sizeof(int));

    // 5. MPI I/O για το col_ind.bin
    sprintf(filepath, "%s/col_ind.bin", dataset_name);
    
    err = MPI_File_open(MPI_COMM_WORLD, filepath, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (err) { printf("Error opening %s\n", filepath); MPI_Abort(MPI_COMM_WORLD, 1); }

    offset = (MPI_Offset)start * sizeof(int);
    
    MPI_File_read_at(fh, offset, A->col_ind, A->local_nnz, MPI_INT, &status);
    MPI_File_close(&fh);

    // 6. Normalization
    long long global_start_offset = A->row_ptr[0];
    for (int i = 0; i <= local_n; i++) {
        A->row_ptr[i] -= global_start_offset;
    }

    return A;
}

void free_csc(CSCMatrix *A) {
    free(A->row_ptr);
    free(A->col_ind);
    free(A);
}

int vector_min(const int *v, int n, int min_val) {
    for (int i = 0; i < n; i++) {
        if (v[i] < min_val)
            min_val = v[i];
    }
    return min_val;
}

void print_vector(int *v, int n){
    printf("[");
    for (int i = 0; i < n; i++) {
        printf("%d ", v[i]);  
    }
    printf("]\n");
}

static int cmp_int(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    if (ia < ib) return -1;
    if (ia > ib) return 1;
    return 0;
}

size_t sort_unique_int(int *v, size_t n) {
    if (n == 0) return 0;

    // 1) sort
    qsort(v, n, sizeof(int), cmp_int);

    size_t w = 1;              // write index
    for (size_t r = 1; r < n; ++r) {   // read index
        if (v[r] != v[w - 1]) {
            v[w] = v[r];
            ++w;
        }
    }

    return w;   // number of unique elements
}

double wall_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char *argv[]) {
    int my_rank, procnum;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procnum);
    
    const char *dataset_name = argv[1];
    double t0=wall_time();
    
    
    
    //load the graph and print the biggest number of neighbors  that a node has
    CSCMatrix *G = load_problem_parallel(my_rank, procnum, dataset_name);
    
    if (G == NULL) {
        if (my_rank == 0) printf("Failed to load graph. Exiting.\n");
        MPI_Finalize();
        return 1;
    }

    const int local_node_number = G->local_n;
    const int node_number = G->n;
    double t1 = wall_time();
    
    if(my_rank == 0) {         
        printf("Time to load the .mat file is %.6f\n", t1-t0);       
    }
    
    int *recvcounts = malloc(procnum * sizeof(int));
    int *rdispls    = malloc(procnum * sizeof(int));
    int offset = 0;
    for(int i = 0 ; i<procnum ; i++) {
        if(i == (procnum-1)){
            recvcounts[i] = (node_number/procnum) + (node_number % procnum);
            rdispls[i]    = offset;
        }
        else  {
            recvcounts[i] = node_number/procnum;
            rdispls[i]    = offset;
        } 
        offset += recvcounts[i];

    }  
    
    //initialize labels
    int *labels = malloc(node_number * sizeof(int));    
    int *labels_check = malloc(local_node_number * sizeof(int));
    for(int i =0; i<node_number; i++ ) {
        labels[i] = i;
    }
    
    int *my_initial_labels = &labels[rdispls[my_rank]];
    memcpy(labels_check, my_initial_labels, local_node_number * sizeof(int));
    

    int exit = 0;
    int global_exit = 0;
    int checking_after_iterations = 0;
    int communicating_after_iterations = 0;

    #pragma omp parallel 
    {       
        #pragma omp single
        {
            if(my_rank == 0) printf("We have %d nodes, %d processes and %d threads\n", node_number, procnum, omp_get_num_threads());
        }

    while(true) {
        
        #pragma omp single
        {    
            communicating_after_iterations++;
            if(communicating_after_iterations == 10){
                communicating_after_iterations = 0;
                MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,  labels, recvcounts, rdispls, MPI_INT, MPI_COMM_WORLD);
            }   
            
        }

        // using ckeck labels to check if finished
        #pragma omp for    
        for (int j = 0; j < local_node_number; j++) {
            // j is the local node  while workers_node is the global one
            long long start = G->row_ptr[j];
            long long end   = G->row_ptr[j+1];
            int workers_node = rdispls[my_rank] + j;                         
            
            if (start == end) { 
                //labels_check[j] = labels[workers_node];
                continue;
            }
            if (labels[workers_node] == 0 ) {
                //labels_check[j] = 0;
                continue;
            }
                                                                                                                        
            int min_val = labels[workers_node];
            for (long long i = start; i < end; ++i) {
                
                int nei = G->col_ind[i];
                int nei_label;
        
                nei_label = labels[nei];
                
                if (nei_label == 0) {
                    min_val = 0;
                    break;
                }
                if (nei_label < min_val) {
                    min_val = nei_label;
                }
            }                                
                //labels_check[j] = labels[workers_node];
            
                labels[workers_node] = min_val;

        }    
                    
        // Check if finished i have to change this for sure because it checks only if one process has ended
        bool changed = false;

        #pragma omp single      
        {    
            checking_after_iterations++;
            if(checking_after_iterations == 20) {
                        
                checking_after_iterations = 0;
                for (int j = 0; j < local_node_number; ++j) {
                    int workers_node = rdispls[my_rank] + j; 
                    if (labels_check[j] != labels[workers_node]) {
                        changed = true;
                        break;
                    }
                }     
                
                if(changed == false) exit = 1;
                else  exit = 0;                         

                MPI_Allreduce(&exit, &global_exit, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                
                // check labels update for next check
                int *my_current_labels = &labels[rdispls[my_rank]];
                memcpy(labels_check, my_current_labels, local_node_number * sizeof(int));

            }
        } 
        
        if(global_exit == procnum) {      //this ensures that every process has finished    
            break;
        } else {
            exit = 0;
            
        }       
    }  
}   
 
    double t2 = wall_time();
    if(my_rank == 0) {
    int unique_labels = sort_unique_int(labels, node_number);
    print_vector(labels, unique_labels);
    //if(n_strongly_connect_componenents == unique_labels) printf("CONGRATS WE ARE CORRECT\nWe have %d unique numbers in labels same as the #SCC\n", unique_labels);
    //else printf("Something went wrong\n");
    printf("We have %d unique numbers in labels same as the #SCC\n", unique_labels);
    printf("Time to propagate is %.6f\n", t2-t1);
    printf("nnz = %lld\n\n",G->nnz);
    }
    
    free(recvcounts);
    free(rdispls); 
    free(labels);
    free_csc(G);
    free(labels_check);
    MPI_Finalize();
    return 0;
}