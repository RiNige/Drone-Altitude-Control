# include <iostream>
using namespace std;

#include "LidDrivenCavity.h"

int main(int argc, char **argv)
{
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    int check = 0;
    try{
        // read input from user and check validity 
        solver->Initialise();
        check = solver->validityCheck();
        if (check == 0){
            throw std::logic_error("dt must smaller than re*dx*dy/4");
        }
        else{
            cout << "Input Parameters pass the validation check" << endl;
        }
    
        // Run the solver
        solver->Integrate();
        solver->GetObj();

        // return the heap memory
        solver->deallocate();
        //solver->~LidDrivenCavity();    // check why destructor has been called twice
        delete solver;
    }
    catch(const std::logic_error & e){
        cout << "Error: " << e.what() << endl;
    }
	return 0;
}