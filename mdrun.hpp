#include <iostream>
#include <string>
#include <cmath>
#include <vector>
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
    Atom(){
    n_atoms++;
    // initialize atom info
    for (int xyz=0; xyz<3; xyz++){
        this->r[xyz] = 0.0; this->v[xyz] = 0.0; this->f[xyz] = 0.0; 
    }
    this->num = n_atoms;
    this->resname = "None";
    this->resnr = 0;
    this->atomname = "None";
    this->mass = 1.0;
    this->epsilon = 0.5;
    this->sigma = 0.5;
    }

    ~Atom(){n_atoms--;}

    void set_resnr(int resnr){this->resnr = resnr;}
    void set_resname(string resname){this->resname = resname;}
    void set_atomname(string atomname){this->atomname = atomname;}
    void set_mass(float mass){this->mass = mass;}
    void set_epsilon(float epsilon){this->epsilon = epsilon;}
    void set_sigma(float sigma){this->sigma = sigma;}

    int get_num(){return this->num;}
    int get_resnr(int resnr){return this->resnr;}
    string get_resname(string resname){return this->resname;}
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

    void integrate_r(float dt){
        for (int xyz=0; xyz<3; xyz++){
            this->r[xyz] += this->v[xyz]*dt;
        }
    }

    void integrate_v(float dt){
        for (int xyz=0; xyz<3; xyz++){
        this->v[xyz] += 0.5*(this->f[xyz]*dt)/(this->mass);
        }
    }

    void reset_f(){
        for (int xyz=0; xyz<3; xyz++){
            this->f[xyz] = 0.0;
        }
    }
};

// for interaction between 2 particles
class Itp{
    int i, j;
    float eps_ij;
    float sig_ij;
    float d_ij;
    float p_ij;
    float f_ij[3];
    float r_ij[3];

public:
    Itp(){}
    ~Itp(){}
    
    void init(Atom *atom, int i, int j){
        this->i = i; 
        this->j = j;
        this->eps_ij = sqrt(atom[i].get_epsilon() * atom[j].get_epsilon());
        this->sig_ij = 0.5*(atom[i].get_sigma() + atom[j].get_epsilon());
    }

    float get_eps_ij(){return this->eps_ij;}
    float get_sig_ij(){return this->sig_ij;}
    void set_r_ij(Atom *atom, float *box){
        // box??
        this->d_ij = 0.0;
        for (int xyz=0; xyz<3; xyz++){
            this->r_ij[xyz] = atom[this->j].get_r(xyz) - atom[this->i].get_r(xyz);
            this->d_ij += pow(r_ij[xyz], 2);
        }
        this->d_ij = sqrt(this->d_ij);
    }

    void set_f_ij(Atom *atom, float *box){
        for (int xyz=0; xyz<3; xyz++){
            this->f_ij[xyz] = 4*(this->eps_ij/this->sig_ij)*(12*pow(this->sig_ij/this->d_ij, 13) - 6*pow(this->sig_ij/this->d_ij, 7))*(r_ij[xyz]/this->d_ij);
            atom[i].add_f(xyz, this->f_ij[xyz]);
            atom[j].add_f(xyz, -this->f_ij[xyz]);
        }
    }

    float set_p_ij(Atom *atom, float *box){
            this->p_ij = 4*this->eps_ij*(pow(this->sig_ij/this->d_ij, 12) - pow(this->sig_ij/this->d_ij, 6));
            return this->p_ij;
    }
};
