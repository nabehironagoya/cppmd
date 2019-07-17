#include <iostream>
#include <vector>
#include <cmath>
#include "mdrun.hpp"
using namespace std;
int Atom::n_atoms = 0;

// for interaction between 2 particles
class Itp{
    int i, j;
    float epsilon;
    float sigma;
    float distance;
public:
//    Itp(Atom *atom, int i, int j){
//        this->epsilon = sqrt(atom[i].get_epsilon() * atom[j].get_epsilon());
//        this->sigma = 0.5*(atom[i].get_sigma() + atom[j].get_epsilon());
//    }
    Itp(){}
    ~Itp(){}
    
    void init(Atom *atom, int i, int j){
        this->i = i; 
        this->j = j;
        this->epsilon = sqrt(atom[i].get_epsilon() * atom[j].get_epsilon());
        this->sigma = 0.5*(atom[i].get_sigma() + atom[j].get_epsilon());
    }

    float get_epsilon();
    float get_sigma();
    float get_distance(Atom *atom, float *box){
        // box??
        this->distance = 0.0;
        for (int xyz=0; xyz<3; xyz++){
            this->distance += pow(atom[this->j].get_r(xyz) - atom[this->i].get_r(xyz), 2);
        }
        this->distance = sqrt(this->distance);
        return this->distance;
    }
};

int main(){
    int n_atoms = 2;
    float box[3];
    Atom *atom = new Atom[n_atoms];
    vector<vector<Itp> > itp(n_atoms, vector<Itp>(n_atoms, Itp()));
    for (int i=0; i<n_atoms; i++){
        for (int j=0; j<i; j++){
            itp[i][j].init(atom, i, j);
        }
    }

    for (int xyz=0; xyz<3; xyz++){
        atom[0].set_r(xyz, 1.0);
        atom[1].set_r(xyz, 0.0);
    }

    cout << itp[1][0].get_distance(atom, box) << endl;
    return 0;
}

