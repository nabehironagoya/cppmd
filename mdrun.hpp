#include <iostream>
#include <string>
using namespace std;


class Atom{
    static int n_atoms;

    int num;
    int resnr;
    string resname;
    string atomname;
    float mass;
    float epsilon;
    float sigma;

    float r[3];
    float v[3];
    float f[3];

public:
    static int get_n_atoms(){return n_atoms;}
    Atom();
    ~Atom(){n_atoms--;}

    void set_atomname(string atomname){this->atomname = atomname;}
    void set_mass(float mass){this->mass = mass;}
    void set_epsilon(float epsilon){this->epsilon = epsilon;}
    void set_sigma(float sigma){this->sigma = sigma;}

    string get_atomname(){return this->atomname;}
    float get_mass(){return this->mass;}
    float get_epsilon(){return this->epsilon;}
    float get_sigma(){return this->sigma;}
    
    float get_r(int xyz){return this->r[xyz];}
    float get_v(int xyz){return this->v[xyz];}
    float get_f(int xyz){return this->f[xyz];}
    void set_r(int xyz, float r){this->r[xyz] = r;}
    void set_v(int xyz, float v){this->v[xyz] = v;}
    void add_f(int xyz, float f){this->f[xyz] += f;}

    void integrate_r(float);
    void integrate_v(float);
    void reset_f();
};

