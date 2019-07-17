#include "mdrun.hpp"
#include "read_gro.hpp"
#include <cmath>
using namespace std;
int Atom::n_atoms = 0;

// from read_gro.cpp
int get_n_atoms(string);
void read_gro(string, int, Gro*, float*);


int main(){
    
    string grofile = "sample.gro";
    float dt = 0.001;
    int n_frames = 10000;
    int n_atoms = get_n_atoms(grofile);
    float box[3];

    Gro *gro = new Gro[n_atoms];
    Atom *atom = new Atom[n_atoms];
    vector<vector<Itp> > itp(n_atoms, vector<Itp>(n_atoms, Itp()));

    read_gro(grofile, n_atoms, gro, box);

    // load atom info from <.gro>
    for (int i=0; i<n_atoms;i++){
        for (int xyz=0; xyz<3; xyz++){
            atom[i].set_r(xyz, gro[i].get_r(xyz));
        }

        atom[i].set_atomname(gro[i].get_atomname());
        atom[i].set_resname(gro[i].get_resname());
        atom[i].set_resnr(gro[i].get_resnr());
    }
    delete [] gro;

    // initialize Itp 
    for (int i=0; i<n_atoms; i++){
        for (int j=0; j<i; j++){
            itp[i][j].init(atom, i, j);
        }
    }

    /* start md */

    // calc force at 0 step    
    float pot = 0.0;
    float kin = 0.0;
    for (int i=0; i<n_atoms; i++){
        for (int j=0; j<i; j++){
            itp[i][j].set_r_ij(atom, box);
            itp[i][j].set_f_ij(atom, box);
            pot += itp[i][j].set_p_ij(atom, box);
        }

        for (int xyz=0; xyz<3; xyz++){
            kin += 0.5*atom[i].get_mass()*pow(atom[i].get_v(xyz), 2);
        }
    }

    for (int frame=0; frame<n_frames; frame++){
        cout << frame <<  " " << pot << " " << kin << " " <<  pot + kin  << endl;
        for (int i=0; i<n_atoms; i++){
            atom[i].integrate_v(dt);
            atom[i].integrate_r(dt);
            atom[i].reset_f();
        }

        // calc force at frame step
        pot = 0.0;
        kin = 0.0;
        for (int i=0; i<n_atoms; i++){
            for (int j=0; j<i; j++){
                itp[i][j].set_r_ij(atom, box);
                itp[i][j].set_f_ij(atom, box);
                pot += itp[i][j].set_p_ij(atom, box);
            }

            for (int xyz=0; xyz<3; xyz++){
                kin += 0.5*atom[i].get_mass()*pow(atom[i].get_v(xyz), 2);
            }
        }

        for (int i=0; i<n_atoms; i++){
            atom[i].integrate_v(dt);
        }
    }

    delete [] atom;
    return 0;
}


