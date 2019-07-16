#include <iostream>
#include <string>
#include <fstream>
using namespace std;

// from string.cpp
string trim(string, char);
vector<string> split(string str, char del);

class Gro{
    static int cnt;
//    static float box[3];
    int resnr;
    string resname;
    string atomname;
    int num;
    float r[3];
public:
    Gro();
    ~Gro();
//    static void set_box(float, float, float);
//    static float *get_box();
    void set_params(int, string, string, int, float, float, float);
    int get_params();
    float get_r(int xyz){return this->r[xyz];}
    int get_resnr(){return this->resnr;}
    int get_num(){return this->num;}
    string get_atomname(){return this->atomname;}
    string get_resname(){return this->resname;}
};

//float Gro::box[] = {0.0, 0.0, 0.0};

//void Gro::set_box(float box_x, float box_y, float box_z){
//    box[0] = box_x;
//    box[1] = box_y;
//    box[2] = box_z;
//}

//float* Gro::get_box(){
//    cout << box[0] << endl;
//    return box;
//}

