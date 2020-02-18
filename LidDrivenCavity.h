# pragma once
# include <string>
# include <fstream>
# include <iomanip>
# include <cmath>
using namespace std;

//-----------------------Class Definition--------------------------------

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
        inline int validityCheck(int n);
        void OutputVal(double* ptr_double, int* ptr_int);
        void InputVal(double* ptr_double, int* ptr_int, int rank);
        void CreateMatrix();
        void Integrate(int rank);
        void currentOmegaBC();
        void currentOmegaInt();
        void nextOmegaInt();
        void GetObj();
        void deallocate();
        void valueCheck();
        friend class Poisson;
    private:
        double** v;        // vorticity
        double** s;        // stream function
        double dt;                   // Time step
        double T;                    // Total time
        int    Nx;                   // Number of grid in x
        int    Ny;                   // Number of grid in y
        int    Px;                   // Number of parallel in x
        int    Py;                   // Number of parallel in y
        double Lx;                   // Length of x
        double Ly;                   // Length of y
        double Re;                   // Reynolds number
        double dx;                   // Distance between each grid in x 
        double dy;                   // Distance between each grid in y
        double U;
};

class Poisson{
    public:
        Poisson();
        void GetVal(double**  ptr1, double**  ptr2);
    private:
        
};


//------------------Poisson class memeber function--------------------------
Poisson::Poisson(){}

void Poisson::GetVal(double**  ptr1, double**  ptr2){
    
}

//-----------------LidDrivenCavity Class member function---------------------

LidDrivenCavity::LidDrivenCavity(){
    double** s = nullptr;
    double** v = nullptr;
    dt = 0.0;
    T = 0.0;
    Nx = 0;
    Ny = 0;
    Px = 0;
    Py = 0;
    Lx = 0.0;
    Ly = 0.0;
    Re = 0.0;
    dx = 0.0;
    dy = 0.0;
    U = 1.0;
}

LidDrivenCavity::~LidDrivenCavity(){
};



void LidDrivenCavity::Initialise(){
    cout << "Please input required parameters from the command line" << endl;
    cin >> Lx >> Ly >> Nx >> Ny >> Px >> Py >> dt >> T >> Re;
    dx = Lx/(Nx - 1);
    dy = Ly/(Ny - 1);
    cout << "User input completed, validating..." << endl;
}

inline int LidDrivenCavity :: validityCheck(int n){      // function is made inline to reduce computation time
    if (dt < Re*dx*dy/4 ){
        if (Px*Py == n){
            return 2;
        }
        else return 1;
    }
    else return 0;
}

void LidDrivenCavity::OutputVal(double* ptr_double, int* ptr_int){
    int newX = Nx/Px, newY = Ny/Py, myMod_x = Nx%Px, myMod_y = Ny%Py;
    double data_double[5] = {dx,dy,dt,T,Re};
    int data_int[6] = {newX,newY,myMod_x,myMod_y,Px,Py};
    for(int i = 0; i < 6; i++){
        ptr_int[i] = data_int[i];
    }
    for(int i = 0; i < 5; i++){
        ptr_double[i] = data_double[i];
    }
}

void LidDrivenCavity::InputVal(double* ptr_double, int* ptr_int, int rank){
    dx = ptr_double[0];
    dy = ptr_double[1];
    dt = ptr_double[2];
    T  = ptr_double[3];
    Re = ptr_double[4];
    if (rank == 0){
        Nx = ptr_int[0]+ptr_int[2];
        Ny = ptr_int[1]+ptr_int[3];
    }
    else if(rank < ptr_int[4]) {
        Nx = ptr_int[0];
        Ny = ptr_int[1]+ptr_int[3];
    }
    else if(rank%ptr_int[4] == 0){
        Nx = ptr_int[0]+ptr_int[2];
        Ny = ptr_int[1];        
    }
    else{
        Nx = ptr_int[0];
        Ny = ptr_int[1];
    }
}

void LidDrivenCavity::CreateMatrix(){
    v = new double*[Ny];
    s = new double*[Ny];
    for(int i = 0; i < Ny; i++){                    //creating 2d vorticity array using heap
        v[i] = new double [Nx];
        s[i] = new double [Nx];
    }
    for(int i = 0; i < Ny; i++){                    // intialise both arrays by initial condition(0)
        for(int j = 0; j < Nx; j++){
            v[i][j] = 0;
            s[i][j] = 0;
        }
    }
}

void LidDrivenCavity::currentOmegaBC(){
    for(int i = 0; i < Nx; i++){
        v[0][i] = (s[0][i] - s[1][i])*2/pow(dy,2) - 2*U/dy; // top boundary
        v[Ny-1][i] = (s[Ny-1][i] - s[Ny-2][i])*2/pow(dy,2); // bottom boundary
    }

    for(int i = 1; i < Ny - 1; i++){
        v[i][0] = (s[i][0] - s[i][1])*2/pow(dx,2); // left
        v[i][Nx-1] = (s[i][Nx-1] - s[i][Nx-2])*2/pow(dx,2); // right
    }
}

void LidDrivenCavity::currentOmegaInt(){
    for(int i = 1; i < Nx - 1; i++){
        for(int j = 0; j < Ny - 1; j++){
            v[i][j] = (s[i][j+1] - 2*s[i][j] + s[i][j-1])/pow(dx,2) + (s[i-1][j] - 2*s[i][j] + s[i+1][j])/pow(dy,2);
        }
    }
}

void LidDrivenCavity:: nextOmegaInt(){
    double part1 = 0, part2 = 0, part3 = 0;
    for(int i = 1; i < Nx - 1; i++){
        for(int j = 0; j < Ny - 1; j++){
            part1 = (s[i-1][j] - s[i+1][j])/(2*dy)*(v[i][j+1] - v[i][j-1])/(2*dx);
            part2 = (s[i][j+1] - s[i][j-1])/(2*dx)*(v[i-1][j] - v[i+1][j])/(2*dy);
            part3 = (v[i][j+1] - 2*v[i][j] + v[i][j-1])/pow(dx,2) + (v[i-1][j] - 2*v[i][j] + v[i+1][j])/pow(dy,2);
            v[i][j] += dt*(1/Re*part3 + part2 - part1);
        }
    }    
}

void LidDrivenCavity::Integrate(int rank){
    double time = 0.0;
    while(time <= T){                                  // only integrate if both array has been given the value
        if(time == 0) cout << "Simulation running..." << endl;
        this->currentOmegaBC();
        this->currentOmegaInt();
        this->nextOmegaInt();
        cout << time << endl;
        time += dt;
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

            Vout << endl << endl;

            for(int i = 0; i < Nx; i++){            // output the vorticity function
                for(int j = 0; j < Ny; j++){
                    Vout << setw(10) << fixed << left << v[i][j];
                }
                Vout << endl;
            }
        }
    Vout.close();
}

void LidDrivenCavity::deallocate(){
    for(int i = 0; i < Ny; i++){
        delete[] v[i];
        delete[] s[i];
    }
    delete[] v;
    delete[] s;
}

void LidDrivenCavity::valueCheck(){
    cout << "Nx = " << Nx << endl << "Ny = " << Ny << endl;
}

