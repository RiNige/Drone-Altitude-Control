# pragma once
# include <string>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <stdlib.h>
using namespace std;

//-----------------------Class Definition--------------------------------

class LidDrivenCavity{
    public:
        LidDrivenCavity();
        ~LidDrivenCavity();
        void Initialise();
        inline int validityCheck(int n);
        void OutputVal(double* ptr_double, int* ptr_int);
        void InputVal(double* ptr_double, int* ptr_int, int rank);
        void CreateMatrix();
        void Integrate(int *coord, int rank,MPI_Comm comm);
        void currentOmegaBC(int *pos);
        void currentOmegaInt(int *pos,int* coord,MPI_Comm comm);
        void nextOmegaInt(int* pos, int* coord,MPI_Comm comm);
        void gather(MPI_Comm comm,int rank,int size,int* ptr);
        void GetObj(int* ptr, int rank);
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




// Empty constructor
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




// Destructor
LidDrivenCavity::~LidDrivenCavity(){
};




// Getting user input from command line
void LidDrivenCavity::Initialise(){
    cout << "Please input required parameters from the command line" << endl;
    cin >> Lx >> Ly >> Nx >> Ny >> Px >> Py >> dt >> T >> Re;
    dx = Lx/(Nx - 1);
    dy = Ly/(Ny - 1);
    cout << "User input completed, validating..." << endl;
}




// Check the correctness of user input
inline int LidDrivenCavity :: validityCheck(int n){      // function is made inline to reduce computation time
    if (dt < Re*dx*dy/4 ){
        if (Px*Py == n){
            if (Px != 1 && Py != 1) {
                return 3;
            }
            else return 2;
        }
        else return 1;
    }
    else return 0;
}




// Output private value
void LidDrivenCavity::OutputVal(double* ptr_double, int* ptr_int){
    int newX = Nx/Px, newY = Ny/Py, myMod_x = Nx%Px, myMod_y = Ny%Py;
    double data_double[5] = {dx,dy,dt,T,Re};
    int data_int[8] = {newX,newY,myMod_x,myMod_y,Px,Py,Nx,Ny};
    for(int i = 0; i < 8; i++){
        ptr_int[i] = data_int[i];
    }
    for(int i = 0; i < 5; i++){
        ptr_double[i] = data_double[i];
    }
}




// Input private value
void LidDrivenCavity::InputVal(double* ptr_double, int* ptr_int, int rank){
    dx = ptr_double[0];
    dy = ptr_double[1];
    dt = ptr_double[2];
    T  = ptr_double[3];
    Re = ptr_double[4];
    Px = ptr_int[4];
    Py = ptr_int[5];
    /* Numbering process in the entire array(e.g. Px=3, Py = 3)
    ---------------------------------------
    |            |            |           |
    |      0     |      1     |      2    |
    |            |            |           |
    ---------------------------------------
    |            |            |           |
    |      3     |      4     |      5    |
    |            |            |           |
    ---------------------------------------
    |            |            |           |
    |      6     |      7     |      8    |
    |            |            |           |
    ---------------------------------------
    */
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




// Create empty vorticity and stream matrix
void LidDrivenCavity::CreateMatrix(){
    srand(time(NULL));
    v = new double*[Ny];
    s = new double*[Ny];
    for(int i = 0; i < Ny; i++){                    //creating 2d vorticity array using heap
        v[i] = new double [Nx];
        s[i] = new double [Nx];
    }
    for(int i = 0; i < Ny; i++){                    // intialise both arrays by initial condition(0)
        for(int j = 0; j < Nx; j++){
            v[i][j] = 0;
            s[i][j] = double (rand())/RAND_MAX;
        }
    }
}




// running the algorithm
void LidDrivenCavity::Integrate(int* coord, int rank,MPI_Comm comm){
    double time = 0.0; 
    int Pos[4] = {0};                               // [top,bot,left,right]           
    if (coord[0] == 0 && coord[1] == 0){                // [1,0,1,0]
        Pos[0] = 1;
        Pos[2] = 1;
    }          
    else if (coord[0] == 0 && coord[1] == (Px - 1)) {     //[1,0,0,1]
        Pos[0] = 1;
        Pos[3] = 1;        
    }
    else if (coord[0] == (Py - 1) && coord[1] == 0){      //[0,1,1,0]
        Pos[1] = 1;
        Pos[2] = 1; 
    }
    else if (coord[0] == (Py - 1) && coord[1] == (Px - 1)){ //[0,1,0,1]
        Pos[1] = 1;
        Pos[3] = 1;    
    }
    else if (coord[0] == 0){                          //[1,0,0,0]
        Pos[0] = 1;
    }
    else if (coord[0] == (Py - 1)){                     //[0,1,0,0]
        Pos[1] = 1;
    }
    else if(coord[1] == 0){                           //[0,0,1,0]
        Pos[2] = 1;
    }
    else if(coord[1] == (Px - 1)){                      //[]0,0,0,1]
        Pos[3] = 1;                                 
    }
    while(time <= T){                                  
        if(time == 0 && rank == 0) cout << "Simulation running..." << endl;
        this->currentOmegaBC(Pos);
        this->currentOmegaInt(Pos,coord,comm);
        this->nextOmegaInt(Pos,coord,comm);
        //if(rank == 0) cout << time << endl;
        time += dt;
    }
    if (rank == 0) cout << "Not finished yet" << endl;
}




// Update boundary vorticity 
void LidDrivenCavity::currentOmegaBC(int* pos){
    // process includes top boundary
    if(pos[0]){                                        
        for(int i = 0; i < Nx; i++){
            v[0][i] = (s[0][i] - s[1][i])*2/pow(dy,2) - 2*U/dy; // update top boundary
        }
    }
    // process includes bottom boundary
    if(pos[1]){
        for(int i = 0; i < Nx; i++){
            v[Ny-1][i] = (s[Ny-1][i] - s[Ny-2][i])*2/pow(dy,2); // update bottom boundary
        }                
    }
    // process includes left boundary
    if(pos[2]){
        for(int i = 0; i < Ny; i++){
            v[i][0] = (s[i][0] - s[i][1])*2/pow(dx,2); // update left boundary
        }
    }
    // process includes right boundary
    if(pos[3]){
        for(int i = 0; i < Ny; i++){
            v[i][Nx-1] = (s[i][Nx-1] - s[i][Nx-2])*2/pow(dx,2); // update right boundary
        }
    }

}




                                    //([top,bot,left,right],[row,colum])
void LidDrivenCavity::currentOmegaInt(int* pos, int* coord,MPI_Comm comm){

    // creating maximum six buffer for top ,bot,left,right
    double* buff_top_recv = new double[Nx];
    double* buff_bot_recv = new double[Nx];
    double* buff_left_send = new double[Ny];
    double* buff_left_recv = new double[Ny];
    double* buff_right_send = new double[Ny];
    double* buff_right_recv = new double[Ny];

    // assign buffer value if matrix boundaries are connected
    if(pos[0] == 0){                                                        //top
        MPI_Send(&(s[0][0]),Nx,MPI_DOUBLE,(coord[0]-1)*Px+coord[1],0,comm);
        MPI_Recv(&(buff_top_recv[0]),Nx,MPI_DOUBLE,(coord[0]-1)*Px+coord[1],1,comm,MPI_STATUS_IGNORE);
    }
    if(pos[1] == 0){                                                        //bot
        MPI_Send(&(s[0][0]),Nx,MPI_DOUBLE,(coord[0]+1)*Px+coord[1],1,comm);
        MPI_Recv(&(buff_bot_recv[0]),Nx,MPI_DOUBLE,(coord[0]+1)*Px+coord[1],0,comm,MPI_STATUS_IGNORE);   
    }
    if(pos[2] == 0){                                                        //left
        for(int i = 0; i < Ny;i++){
            buff_left_send[i] = s[i][0];
        }
        MPI_Send(&(buff_left_send[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]-1,2,comm);
        MPI_Recv(&(buff_left_recv[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]-1,3,comm,MPI_STATUS_IGNORE);   
    }    
    if(pos[3] == 0){                                                        //right
        for(int i = 0; i < Ny;i++){
            buff_right_send[i] = s[i][0];
        }
        MPI_Send(&(buff_right_send[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]+1,3,comm);
        MPI_Recv(&(buff_right_recv[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]+1,2,comm,MPI_STATUS_IGNORE);    
    }  
    MPI_Barrier(comm);

    // update inter value and leave four matrix boundary untouched
    for(int i = 1; i < (Ny - 1);i++){
        for(int j = 1; j < (Nx - 1); j++){
            v[i][j] = -((s[i][j+1] - 2*s[i][j] + s[i][j-1])/pow(dx,2) + (s[i-1][j] - 2*s[i][j] + s[i+1][j])/pow(dy,2));
        }
    }

    //update four outlayer without the corner
    if(pos[0] == 0){ // if top row is not boundary
        for(int i = 1; i < (Nx - 1); i++){
            v[0][i] = -((s[0][i+1] - 2*s[0][i] + s[0][i-1])/pow(dx,2) + (buff_top_recv[i] - 2*s[0][i] + s[1][i])/pow(dy,2));
        }
    }
    if(pos[1] == 0){ // if bot row is not boundary
        for(int i = 1; i < (Nx - 1); i++){
            v[Ny-1][i] = -((s[Ny-1][i+1] - 2*s[Ny-1][i] + s[Ny-1][i-1])/pow(dx,2) + (s[Ny-2][i] - 2*s[Ny-1][i] + buff_bot_recv[i])/pow(dy,2));
        }
    }
    if(pos[2] == 0){ // if left column is not boundary
        for(int i = 1; i < (Ny - 1); i++){
            v[i][0] = -((s[i][1] - 2*s[i][0] + buff_left_recv[i])/pow(dx,2) + (s[i-1][0] - 2*s[i][0] + s[i+1][0])/pow(dy,2));
        }
    }
    if(pos[3] == 0){ // if right solumn is not boundary
        for(int i = 1; i < (Ny - 1); i++){
            v[i][Nx-1] = -((buff_right_recv[i] - 2*s[i][Nx-1] + s[i][Nx-2])/pow(dx,2) + (s[i-1][Nx-1] - 2*s[i][Nx-1] + s[i+1][Nx-1])/pow(dy,2));
        }
    }

    //update four corner
    if(pos[0] == 0 && pos[2] == 0){ // update top left corner
        v[0][0] = -((s[0][1] - 2*s[0][0] + buff_left_recv[0])/pow(dx,2) + (buff_top_recv[0] - 2*s[0][0] + s[1][0])/pow(dy,2));
    }
    if(pos[0] == 0 && pos[3] == 0){ // update top right corner
        v[0][Nx-1] = -((buff_right_recv[0] - 2*s[0][Nx-1] + s[0][Nx-2])/pow(dx,2) + (buff_top_recv[Nx-1] - 2*s[0][Nx-1] + s[1][Nx-1])/pow(dy,2));
    }
    if(pos[1] == 0 && pos[2] == 0){ // update bottom left corner
        v[Ny-1][0] = -((s[Ny-1][1] - 2*s[Ny-1][0] + buff_left_recv[Ny-1])/pow(dx,2) + (s[Ny-2][0] - 2*s[Ny-1][0] + buff_top_recv[0])/pow(dy,2));
    }
    if(pos[1] == 0 && pos[3] == 0){ // update bottom right corner
        v[Ny-1][Nx-1] = -((buff_right_recv[Ny-1] - 2*s[Ny-1][Nx-1] + s[Ny-1][Nx-2])/pow(dx,2) + (s[Ny-2][Nx-1] - 2*s[Ny-1][Nx-1] + buff_top_recv[Nx-1])/pow(dy,2));
    }

    // returning memory at the end of process
    delete[] buff_top_recv;
    delete[] buff_bot_recv;
    delete[] buff_left_send;
    delete[] buff_left_recv;
    delete[] buff_right_send;
    delete[] buff_right_recv;
}



                                    //([top,bot,left,right],[row,colum],mpi_comm)
void LidDrivenCavity:: nextOmegaInt(int* pos, int* coord, MPI_Comm comm){

    // creating maximum six buffer for top ,bot,left,right
    double* buff_top_recv = new double[2*Nx];
    double* buff_bot_recv = new double[2*Nx];
    double* buff_left_send = new double[2*Ny];
    double* buff_left_recv = new double[2*Ny];
    double* buff_right_send = new double[2*Ny];
    double* buff_right_recv = new double[2*Ny];

    // temp variable used in the code
    double part1 = 0, part2 = 0, part3 = 0;
    double ** v_temp = new double *[Ny];            // temped variable used to store next step omega value
    for(int i = 0; i < Ny; i++){     
        v_temp[i] = new double[Nx];
    }

    for(int i = 0; i < Ny; i++){
        for(int j = 0; j < Nx; j++){
            v_temp[i][j] = 0.0;
        }
    }

    // assign buffer value if matrix boundaries are connected
    if(pos[0] == 0){                                                        //top
        MPI_Send(&(s[0][0]),Nx,MPI_DOUBLE,(coord[0]-1)*Px+coord[1],0,comm);
        MPI_Send(&(v[0][0]),Nx,MPI_DOUBLE,(coord[0]-1)*Px+coord[1],1,comm);
        MPI_Recv(&(buff_top_recv[0]),Nx,MPI_DOUBLE,(coord[0]-1)*Px+coord[1],0,comm,MPI_STATUS_IGNORE);
        MPI_Recv(&(buff_top_recv[Nx]),Nx,MPI_DOUBLE,(coord[0]-1)*Px+coord[1],1,comm,MPI_STATUS_IGNORE);
    }
    if(pos[1] == 0){                                                        //bot
        MPI_Send(&(s[0][0]),Nx,MPI_DOUBLE,(coord[0]+1)*Px+coord[1],0,comm);
        MPI_Send(&(v[0][0]),Nx,MPI_DOUBLE,(coord[0]+1)*Px+coord[1],1,comm);
        MPI_Recv(&(buff_bot_recv[0]),Nx,MPI_DOUBLE,(coord[0]+1)*Px+coord[1],0,comm,MPI_STATUS_IGNORE);   
        MPI_Recv(&(buff_bot_recv[Nx]),Nx,MPI_DOUBLE,(coord[0]+1)*Px+coord[1],1,comm,MPI_STATUS_IGNORE); 
    }
    if(pos[2] == 0){                                                        //left
        for(int i = 0; i < 2*Ny;i++){
            if(i < Ny) buff_left_send[i] = s[i][0];
            else buff_left_send[i] = v[i-Ny][0];
        }
        MPI_Send(&(buff_left_send[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]-1,0,comm);
        MPI_Send(&(buff_left_send[Ny]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]-1,1,comm);
        MPI_Recv(&(buff_left_recv[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]-1,0,comm,MPI_STATUS_IGNORE);  
        MPI_Recv(&(buff_left_recv[Ny]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]-1,1,comm,MPI_STATUS_IGNORE);  
    }    
    if(pos[3] == 0){                                                        //right
        for(int i = 0; i < 2*Ny;i++){
            if(i < Ny) buff_right_send[i] = s[i][0];
            else buff_right_send[i] = v[i-Ny][0];
        }
        MPI_Send(&(buff_right_send[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]+1,0,comm);
        MPI_Send(&(buff_right_send[Ny]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]+1,1,comm);
        MPI_Recv(&(buff_right_recv[0]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]+1,0,comm,MPI_STATUS_IGNORE); 
        MPI_Recv(&(buff_right_recv[Ny]),Ny,MPI_DOUBLE,coord[0]*Px+coord[1]+1,1,comm,MPI_STATUS_IGNORE);    
    }  
    MPI_Barrier(comm);


    // running the update for inner blocks using local value
    for(int i = 1; i < Ny - 1; i++){
        for(int j = 1; j < Nx - 1; j++){
            part1 = (s[i-1][j] - s[i+1][j])/(2*dy)*(v[i][j+1] - v[i][j-1])/(2*dx);
            part2 = (s[i][j+1] - s[i][j-1])/(2*dx)*(v[i-1][j] - v[i+1][j])/(2*dy);
            part3 = (v[i][j+1] - 2*v[i][j] + v[i][j-1])/pow(dx,2) + (v[i-1][j] - 2*v[i][j] + v[i+1][j])/pow(dy,2);
            v_temp[i][j] = v[i][j] + dt*(1/Re*part3 + part2 - part1);
        }
    } 

    //update four boundary without the corner
    if(pos[0] == 0){ // if top row is not boundary
        for(int i = 1; i < (Nx - 1); i++){
            part1 = (buff_top_recv[i] - s[1][i])/(2*dy)*(v[0][i+1] - v[0][i-1])/(2*dx);
            part2 = (s[0][i+1] - s[0][i-1])/(2*dx)*(buff_top_recv[Nx+i] - v[1][i])/(2*dy);
            part3 = (v[0][i+1] - 2*v[0][i] + v[0][i-1])/pow(dx,2) + (buff_top_recv[Nx+i] - 2*v[0][i] + v[1][i])/pow(dy,2);
            v_temp[0][i] = v[0][i] + dt*(1/Re*part3 + part2 - part1);
        }
    }
    if(pos[1] == 0){ // if bot row is not boundary
        for(int i = 1; i < (Nx - 1); i++){
            part1 = (s[Ny-2][i] - buff_bot_recv[i])/(2*dy)*(v[Ny-1][i+1] - v[Ny-1][i-1])/(2*dx);
            part2 = (s[Ny-1][i+1] - s[Ny-1][i-1])/(2*dx)*(v[Ny-2][i] - buff_bot_recv[Nx+i])/(2*dy);
            part3 = (v[Ny-1][i+1] - 2*v[Ny-1][i] + v[Ny-1][i-1])/pow(dx,2) + (v[Ny-2][i] - 2*v[Ny-1][i] + buff_bot_recv[Nx+i])/pow(dy,2);
            v_temp[Ny-1][i] = v[Ny-1][i] + dt*(1/Re*part3 + part2 - part1);
        }
    }
    if(pos[2] == 0){ // if left column is not boundary
        for(int i = 1; i < (Ny - 1); i++){
            part1 = (s[i-1][0] - s[i+1][0])/(2*dy)*(v[i][1] - buff_left_recv[Ny+i])/(2*dx);
            part2 = (s[i][1] - buff_left_recv[i])/(2*dx)*(v[i-1][0] - v[i+1][0])/(2*dy);
            part3 = (v[i][1] - 2*v[i][0] + buff_left_recv[Ny+i])/pow(dx,2) + (v[i-1][0] - 2*v[i][0] + v[i+1][0])/pow(dy,2);
            v_temp[i][0] = v[i][0] + dt*(1/Re*part3 + part2 - part1);
        }
    }
    if(pos[3] == 0){ // if right solumn is not boundary
        for(int i = 1; i < (Ny - 1); i++){
            part1 = (s[i-1][Nx-1] - s[i+1][Nx-1])/(2*dy)*(buff_right_recv[Ny+i] - v[i][Nx-2])/(2*dx);
            part2 = (buff_right_recv[i] - s[i][Nx-2])/(2*dx)*(v[i-1][Nx-1] - v[i+1][Nx-1])/(2*dy);
            part3 = (buff_right_recv[Ny+i] - 2*v[i][Nx-1] + v[i][Nx-2])/pow(dx,2) + (v[i-1][Nx-1] - 2*v[i][Nx-1] + v[i+1][Nx-1])/pow(dy,2);
            v_temp[i][Nx-1] = v[i][Nx-1] + dt*(1/Re*part3 + part2 - part1);
        }
    } 

    //update four corner
    if(pos[0] == 0 && pos[2] == 0){ // update top left corner
            part1 = (buff_top_recv[0] - s[1][0])/(2*dy)*(v[0][1] - buff_left_recv[Ny])/(2*dx);
            part2 = (s[0][1] - buff_left_recv[0])/(2*dx)*(buff_top_recv[Nx] - v[1][0])/(2*dy);
            part3 = (v[0][1] - 2*v[0][0] + buff_left_recv[Ny])/pow(dx,2) + (buff_top_recv[Nx] - 2*v[0][0] + v[1][0])/pow(dy,2);
            v_temp[0][0] = v[0][0] + dt*(1/Re*part3 + part2 - part1);
    }
    if(pos[0] == 0 && pos[3] == 0){ // update top right corner
            part1 = (buff_top_recv[Nx-1] - s[1][Nx-1])/(2*dy)*(buff_right_recv[Ny] - v[0][Nx-2])/(2*dx);
            part2 = (buff_right_recv[0] - s[0][Nx-2])/(2*dx)*(buff_top_recv[2*Nx-1] - v[1][Nx-1])/(2*dy);
            part3 = (buff_right_recv[Ny] - 2*v[0][Nx-1] + v[0][Nx-2])/pow(dx,2) + (buff_top_recv[2*Nx-1] - 2*v[0][Nx-1] + v[1][Nx-1])/pow(dy,2);
            v_temp[0][Nx-1] = v[0][Nx-1] + dt*(1/Re*part3 + part2 - part1);
    }
    if(pos[1] == 0 && pos[2] == 0){ // update bottom left corner
            part1 = (s[Ny-2][0] - buff_bot_recv[0])/(2*dy)*(v[Ny-1][1] - buff_left_recv[2*Ny-1])/(2*dx);
            part2 = (s[Ny-1][1] - buff_left_recv[Ny-1])/(2*dx)*(v[Ny-2][0] - buff_bot_recv[Nx])/(2*dy);
            part3 = (v[Ny-1][1] - 2*v[Ny-1][0] + buff_left_recv[2*Ny-1])/pow(dx,2) + (v[Ny-2][0] - 2*v[Ny-1][0] + buff_bot_recv[Nx])/pow(dy,2);
            v_temp[Ny-1][0] = v[Ny-1][0] + dt*(1/Re*part3 + part2 - part1);
    }
    if(pos[1] == 0 && pos[3] == 0){ // update bottom right corner
            part1 = (s[Ny-2][Nx-1] - buff_bot_recv[Nx-1])/(2*dy)*(buff_right_recv[2*Ny-1] - v[Ny-1][Nx-2])/(2*dx);
            part2 = (buff_right_recv[Ny-1] - s[Ny-1][Nx-2])/(2*dx)*(v[Ny-2][Nx-1] - buff_bot_recv[2*Nx-1])/(2*dy);
            part3 = (buff_right_recv[2*Ny-1] - 2*v[Ny-1][Nx-1] + v[Ny-1][Nx-2])/pow(dx,2) + (v[Ny-2][Nx-1] - 2*v[Ny-1][Nx-1]+ buff_bot_recv[2*Nx-1])/pow(dy,2);
            v_temp[Ny-1][Nx-1] = v[Ny-1][Nx-1] + dt*(1/Re*part3 + part2 - part1);
    }

    // substitude v with v_temp
    int begin_x = 0, begin_y = 0,end_x = Nx, end_y = Ny;
    if(pos[0] == 1) begin_y = 1;
    if(pos[1] == 1) end_y = Ny - 1;
    if(pos[2] == 1) begin_x = 1;
    if(pos[3] == 1) end_x = Nx - 1;
    for(int i = begin_y; i < end_y;i++){
        for(int j = begin_x; j < end_x; j++){
            v[i][j]= v_temp[i][j];
        }
    }

    // returning memory at the end of process
    delete[] buff_top_recv;
    delete[] buff_bot_recv;
    delete[] buff_left_send;
    delete[] buff_left_recv;
    delete[] buff_right_send;
    delete[] buff_right_recv;   
}




/* void LidDrivenCavity::gather(MPI_Comm comm,int rank,int size,int* ptr){
    double v_final[161][161] = {1};
    double s_final[161][161] = {1};
    int Nx_tot = 0, Ny_tot = 0, Nx_sp = 0, Ny_sp = 0;
    if (rank == 0){
        Nx_sp = Nx;
        Ny_sp = Ny;
    }
    MPI_Bcast(&(Nx_sp),1,MPI_INT,0,comm);
    MPI_Bcast(&(Ny_sp),1,MPI_INT,0,comm);
    if(rank == 0){
        for(int i = 0; i < Ny; i++){
            for(int j = 0; j < Nx; j++){
                v_final[i][j] = v[i][j];
                s_final[i][j] = s[i][j];
            }
        }
        for(int i = 1; i <= size; i++){
            MPI_Recv(&(v_final[Ny_sp + (ptr[0]-1)*(Ny-1)][Nx_sp + (ptr[1]-1)*(Ny-1)]),Ny*Nx,MPI_DOUBLE,i,i*2-1,comm,MPI_STATUS_IGNORE);
            MPI_Recv(&(s_final[Ny_sp + (ptr[0]-1)*(Ny-1)][Nx_sp + (ptr[1]-1)*(Ny-1)]),Ny*Nx,MPI_DOUBLE,i,i*2,comm,MPI_STATUS_IGNORE);
        }
    }
    else{
        MPI_Send(&(v[0][0]),Ny*Nx,MPI_DOUBLE,0,rank*2-1,comm);
        MPI_Send(&(s[0][0]),Ny*Nx,MPI_DOUBLE,0,rank*2,comm);
    }
    MPI_Reduce(&(Nx),&(Nx_tot),1,MPI_INT,MPI_SUM,0,comm);
    MPI_Reduce(&(Nx),&(Ny_tot),1,MPI_INT,MPI_SUM,0,comm);
    Nx_tot /= Px;
    Ny_tot /= Py;
    if(rank == 0) cout <<"This is final Nx and Ny: " << Nx_tot << ' ' << Ny_tot << endl;
    if (rank == 0){
        cout <<"This is final Nx and Ny: " << Nx_tot << ' ' << Ny_tot << endl;
        ofstream Vout("data.dat",ios::out|ios::trunc);
            Vout.precision(5);
            if(Vout.good()){                            // check if the file is ready to write
                for(int i = 0; i < Ny_tot; i++){            // output the stream function
                    for(int j = 0; j < Nx_tot; j++){
                        Vout << setw(10) << fixed << left << s_final[i][j];
                    }
                    Vout << endl;
                }

                Vout << endl << endl;

                for(int i = 0; i < Ny_tot; i++){            // output the vorticity function
                    for(int j = 0; j < Nx_tot; j++){
                        Vout << setw(10) << fixed << left << v_final[i][j];
                    }
                    Vout << endl;
                }
            }
            Vout.close();
    }
} */




void LidDrivenCavity::GetObj(int* ptr, int rank){
    string row = to_string(ptr[0]);
    string colum = to_string(ptr[1]);
    string filename = "rank_" + to_string(rank) + "_coord_" + row + '_' + colum;
    ofstream Vout(filename.c_str(),ios::out|ios::trunc);
        Vout.precision(5);
        if(Vout.good()){                            // check if the file is ready to write
            for(int i = 0; i < Ny; i++){            // output the stream function
                for(int j = 0; j < Nx; j++){
                    Vout << setw(20) << fixed << left << s[i][j];
                }
                Vout << endl;
            }

            Vout << endl << endl;

            for(int i = 0; i < Ny; i++){            // output the vorticity function
                for(int j = 0; j < Nx; j++){
                    Vout << setw(20) << left << v[i][j];
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

