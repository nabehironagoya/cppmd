#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

// from string.cpp
string trim(string, char);
vector<string> split(string str, char del);

class Gro{
    static int cnt;
//    static double box[3];
    int resnr;
    string resname;
    string atomname;
    int num;
    double r[3];
public:
    Gro();
    ~Gro();
//    static void set_box(double, double, double);
//    static double *get_box();
    void set_params(int, string, string, int, double, double, double);
    int get_params();
    double get_r(int xyz){return this->r[xyz];}
    int get_resnr(){return this->resnr;}
    int get_num(){return this->num;}
    string get_atomname(){return this->atomname;}
    string get_resname(){return this->resname;}
};

//double Gro::box[] = {0.0, 0.0, 0.0};

//void Gro::set_box(double box_x, double box_y, double box_z){
//    box[0] = box_x;
//    box[1] = box_y;
//    box[2] = box_z;
//}

//double* Gro::get_box(){
//    cout << box[0] << endl;
//    return box;
//}

