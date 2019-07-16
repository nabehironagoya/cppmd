#include "mdrun.hpp"
#include "read_gro.hpp"


// from read_gro.cpp
int get_n_atoms(string);
void read_gro(string, int, Gro*, float*);


int do_md(){
    
    string grofile = "em.gro";
    float dt = 0.001;
    int n_frames = 10000;
    int n_atoms = get_n_atoms(grofile);
    float box[3];

    Gro *gro = new Gro[n_atoms];
    Atom *atom = new Atom[n_atoms];

    read_gro(grofile, n_atoms, gro, box);

    cout << n_atoms << endl;
    for (int i=0; i<n_atoms;i++){
        gro[i].get_params();
    }
    cout << box[0] ;
    cout << gro[0].get_r(2);

    delete [] gro;


//    atom[0].add_f(0, -k*atom[0].get_r(0));
//    for (int frame=0; frame<n_frames; frame++){
//        atom[0].integrate_v(dt);
//        atom[0].integrate_r(dt);
//        atom[0].reset_f();
//        atom[0].add_f(0, -k*atom[0].get_r(0));
//        atom[0].integrate_v(dt);
//
//    }

    delete [] atom;
    return 0;

}

int main(){
    do_md();
}

int Atom::n_atoms = 0;

Atom::Atom(){
    n_atoms++;
    // initialize atom info
    for (int xyz=0; xyz<3; xyz++){
        this->r[xyz] = 0.0; this->v[xyz] = 0.0; this->f[xyz] = 0.0; } 
    this->mass = 0.0;
    this->epsilon = 0.0;
    this->sigma = 0.0;
    this->atomname = "None";
    }

inline void Atom::integrate_r(float dt){
    for (int xyz=0; xyz<3; xyz++){
        this->r[xyz] += this->v[xyz]*dt;
    }
}

inline void Atom::integrate_v(float dt){
    for (int xyz=0; xyz<3; xyz++){
        this->v[xyz] += 0.5*(this->f[xyz]*dt)/(this->mass);
    }
}

inline void Atom::reset_f(){
    for (int xyz=0; xyz<3; xyz++){
        this->f[xyz] = 0.0;
    }
}

