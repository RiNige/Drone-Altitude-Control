# include <iostream>
# include <mpi.h>
#include "LidDrivenCavity.h"

using namespace std;

int main(int argc, char *argv[]){
    // MPI configuretion
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
        int* data_int = new int[6];        // store int type input for each process

        // read input from user and check validity on root process
        if (rank == root){
            solver->Initialise();
            int check = solver->validityCheck(size);
            if (check == 0){
                throw std::logic_error("dt must smaller than re*dx*dy/4");
            }
            else if(check == 1){
                throw std::logic_error("Partitons in x and y are not compatible with the number of process");
            }
            else{
                cout << "Input Parameters pass the validation check" << endl;
            }
            // extract user input from class and broadcast to other processes
            solver->OutputVal(data_doub,data_int);
        }
        MPI_Bcast(data_doub,5,MPI_DOUBLE,root,MPI_COMM_WORLD);
        MPI_Bcast(data_int,6,MPI_INT,root,MPI_COMM_WORLD);
        solver->InputVal(data_doub,data_int,rank);
        cout << "This is rank: " << rank << endl;
        solver->valueCheck();
        
        // Create separated matrix in each process
        solver->CreateMatrix();

        // Run the solver
        //solver->Integrate(rank);
        //solver->GetObj();

        // return memory
        solver->deallocate();
        solver->~LidDrivenCavity(); 
        delete[] data_doub;
        delete[] data_int;
        delete solver;
    }
    catch(const std::logic_error & e){
        cout << "An error has occured due to following reason: " << e.what() << endl;
    }
    MPI_Finalize();
	return 0;
}