#pragma once

# include <string>
# include <fstream>
# include <iomanip>
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

    void Initialise();
    void Integrate();
    inline int validityCheck();
    void currentOmegaBC();
    void GetObj();
    void deallocate();

private:
    double** v = nullptr;        // vorticity
    double** s = nullptr;        // stream function

    double dt;                   // Time step
    double T;                    // Total time
    int    Nx;                   // Number of grid in x
    int    Ny;                   // Number of grid in y
    double Lx;                   // Length of x
    double Ly;                   // Length of y
    double Re;                   // Reynolds number
    double dx;                   // Distance between each grid in x 
    double dy;                   // Distance between each grid in y
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
    cout << "A constructor has been called" << endl;
}

LidDrivenCavity::~LidDrivenCavity(){
    cout << "A destructor has been called" << endl;
};



void LidDrivenCavity::Initialise(){
    std::cout << "Please input required parameters from the command line" << std::endl;
    cin >> Lx >> Ly >> Nx >> Ny >> dt >> T >> Re;
    dx = Lx/(Nx - 1);
    dy = Ly/(Ny - 1);
    std::cout << "User input completed, validating..." << std::endl;
    v = new double*[Ny];
    for(int i = 0; i < Nx; i++){           //creating 2d vorticity array using heap
        v[i] = new double [Nx];
    }
    s = new double*[Ny];
    for(int i = 0; i < Nx; i++){           //creating 2d stream function array using heap
        s[i] = new double [Nx];
    }
    for(int i = 0; i < Nx; i++){            // intialise both arrays by initial condition(0)
        for(int j = 0; j < Ny; j++){
            v[i][j] = 0;
            s[i][j] = 0;
        }
    }

}

inline int LidDrivenCavity :: validityCheck(){   // function is made inline to reduce computation time
    if (dt < Re*dx*dy/4) return 1;
    else return 0;
}

void LidDrivenCavity::currentOmegaBC(){

}

void LidDrivenCavity::Integrate(){
    double time = 0;
    while(v && s){                 // check if v is nullptr
    cout << "v and s are not nullptrs" << endl;
    break;
    }
    std::cout << "Not finished yet" << std::endl;
}

void LidDrivenCavity::GetObj(){
    ofstream Vout("Stream.dat",ios::out|ios::trunc);
        Vout.precision(5);
        if(Vout.good()){                            // check if the file is ready to write
            for(int i = 0; i < Nx; i++){            // output the stream function
                for(int j = 0; j < Ny; j++){
                    Vout << setw(10) << fixed << left << s[i][j];
                }
                Vout << endl;
            }
        }
    Vout.close();
}

void LidDrivenCavity::deallocate(){
    for(int i = 0; i < Nx; i++){
        delete[] v[i];
        delete[] s[i];
    }
    delete[] v;
    delete[] s;
}

