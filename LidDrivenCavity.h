#pragma once

#include <string>
using namespace std;

class LidDrivenCavity{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    /*void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);*/

    double** Initialise();
    void Integrate();
    inline int validityCheck();
    void currentOmegaBC();
    void deallocate(double** ptr);

    // Add any other public functions

private:
    double* v = nullptr;
    double* s = nullptr;

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    double dx;
    double dy;
};


LidDrivenCavity::LidDrivenCavity(){
    dt = 0;
    T = 0;
    Nx = 0;
    Ny = 0;
    Lx = 0;
    Ly = 0;
    Re = 0;
    dx = 0;
    dy = 0;
}

LidDrivenCavity::~LidDrivenCavity(){

    delete this;
}



double** LidDrivenCavity::Initialise(){
    std::cout << "Please input required parameters from the command line" << std::endl;
    cin >> Lx >> Ly >> Nx >> Ny >> dt >> T >>Re;
    dx = Lx/(Nx - 1);
    dy = Ly/(Ny - 1);
    std::cout << "User input completed, validating..." << std::endl;
    double **cavity = new double*[Ny];
    for(int i = 0; i < Nx; i++){
        cavity[i] = new double [Nx];
    }
    for(int i = 0; i < Nx; i++){
        for(int j = 0; j < Ny; j++){
            cavity[j][i] = 0;
        }
    }
    return cavity;

}

inline int LidDrivenCavity :: validityCheck(){
    if (dt < Re*dx*dy/4) return 1;
    else return 0;
}

void LidDrivenCavity::currentOmegaBC(){

}

void LidDrivenCavity::Integrate(){
    std::cout << "Not finished yet" << std::endl;
}

void LidDrivenCavity::deallocate(double** ptr){
    for(int i = 0; i < Nx; i++){
        delete[] ptr[i];
    }
    delete[] ptr;
}

