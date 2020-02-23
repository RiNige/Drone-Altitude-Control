# include <iostream>
# include <mpi.h>
# include "LidDrivenCavity.h"

using namespace std;

int main(int argc, char *argv[]){


    //------------------------------------- MPI configuretion-----------------------------------------


    int rank = 0, size = 0, root = 0; // specify rank 0 as root processor
    int err = MPI_Init(&argc, &argv);
    if(err != MPI_SUCCESS){
        cout << "There is an error when intialising MPI" << endl;
        return -1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    try{
        double* data_doub = new double[5]; // store double type input for each process
        int* data_int = new int[8];        // store int type input for each process


        // ----------------read input from user and check validity on root process-----------------------


        if (rank == root){
            solver->Initialise();
            int check = solver->validityCheck(size);

            if (check == 0){
                throw std::logic_error("dt must smaller than re*dx*dy/4");
            }
            else if(check == 1){
                throw std::logic_error("Partitons in x and y are not compatible with the number of process");
            }
            else if(check == 2){
                throw std::logic_error("Partitons in x and y cannot be 1");
            }
            else{
                cout << "Input Parameters pass the validation check" << endl;
            }

            // extract user input from class and broadcast to other processes
            solver->OutputVal(data_doub,data_int);
        }

        // all processes wait until root process finish reading data
        MPI_Barrier(MPI_COMM_WORLD);


        // ----------------Boradcast input to all processed and intialise matrix-----------------------


        MPI_Bcast(data_doub,5,MPI_DOUBLE,root,MPI_COMM_WORLD);
        MPI_Bcast(data_int,8,MPI_INT,root,MPI_COMM_WORLD);
        MPI_Comm comm;
        int dim[2] = {data_int[5], data_int[4]};
        int period[2] = {0,0}, coord[2];

        // Topology Optimisation and rearrange the metrix
        MPI_Cart_create(MPI_COMM_WORLD,2,dim,period,1,&comm);
        MPI_Cart_coords(comm, rank, 2, coord);
        //cout << "This is rank" << rank << endl;
        //cout << "Coordinates after toplogy: " << endl;
        //cout << coord[0] << ' '<< coord[1] << endl;
        

        // Read data broadcasted from root
        solver->InputVal(data_doub,data_int,rank);

        //check accuracy of individual rank dimension, then create matrix
        solver->valueCheck();
        solver->CreateMatrix();
        MPI_Barrier(comm);


        // --------------------Run the solver in parallel-----------------------------------------------


        solver->Integrate(coord,rank,comm);


        // --------Collect vorticity and stream matrix from each process and output in one dat file------


        // extract two matrices on each process and gather 
        //solver->gather(comm,rank,size,coord);
        //MPI_Barrier(comm);
        solver->GetObj(coord,rank);


        //-------------------------Return memory--------------------------------------------------------


        solver->deallocate();
        solver->~LidDrivenCavity();

        delete[] data_doub;
        delete[] data_int;
    }
    catch(const std::logic_error & e){
        cout << "An error has occured due to following reason: " << e.what() << endl;
        cout << "Please terminate the program" << endl;
    }
    delete solver;
    MPI_Finalize();
	return 0;
}