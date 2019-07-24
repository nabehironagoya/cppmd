#include "mdrun.hpp"
#include "gro_util.hpp"
#include <cmath>
#include <fstream>

#ifdef OPENMP
#include <omp.h>
#endif /* OPENMP */

using namespace std;
int Atom::n_atoms = 0;

// from read_gro.cpp
int get_n_atoms(string);
void read_gro(string, int, Gro*, double*);
void write_gro(string, int, Gro*, double*);


int main(){
    
	//declarations
    string grofile = "sample.gro";
    double dt = 0.001;
    int n_frames = 10000;
    int n_atoms = get_n_atoms(grofile);
    double r_c = 2.0;
    double box[3];
    double pot = 0.0;
    double kin = 0.0;
	double T = 0.0;
	double T_ref = 100.0;
	double Q = 1.0;
	double P = 0.0;
	double vir = 0.0;
	double zeta = 0.0;
	double eta = 0.0;
	double H = 0.0;
	double s = 0.0;
	const double R = 0.00831; // kJ/mol K


	ofstream ofs_coord("coord.xyz");
	ofstream ofs_energy("energy.dat");

    Gro *gro = new Gro[n_atoms];
    Atom *atom = new Atom[n_atoms];
    vector<vector<Itp> > itp(n_atoms, vector<Itp>(n_atoms, Itp()));


    // load atom info from <.gro>
    read_gro(grofile, n_atoms, gro, box);
    for (int i=0; i<n_atoms;i++){
        for (int xyz=0; xyz<3; xyz++){
            atom[i].set_r(xyz, gro[i].get_r(xyz));
        }

        atom[i].set_atomname(gro[i].get_atomname());
        atom[i].set_resname(gro[i].get_resname());
        atom[i].set_resnr(gro[i].get_resnr());
    }

    // initialize Itp 
    for (int i=0; i<n_atoms; i++){
        for (int j=0; j<i; j++){
            itp[i][j].init(atom, i, j);
        }
    }

    /* start md */

    // calc force at 0 step    
    pot = 0.0;
    kin = 0.0;
    for (int i=0; i<n_atoms; i++){
        for (int j=0; j<i; j++){
            itp[i][j].set_r_ij(atom, box);
            itp[i][j].set_f_ij(atom, r_c);
            pot += itp[i][j].set_p_ij(atom, r_c);
        }
        kin += atom[i].get_kin();

    }
	
	// energy calc at 0 step
	T = (2.0*kin)/(3.0*n_atoms*R);
	vir = 0.0;
    for (int i=0; i<n_atoms; i++){
        for (int j=0; j<i; j++){
			for (int xyz=0; xyz<3 ; xyz++){
				vir += atom[i].get_r(xyz)*atom[i].get_f(xyz);
			}
		}
	}
	P = (2.0*kin + vir)/(3*box[0]*box[1]*box[2]);
	H = kin + pot + 0.5*Q*pow(zeta, 2) + 3*n_atoms*R*T_ref*eta;

    for (int frame=0; frame<n_frames; frame++){
        cout << frame <<" " <<  H  << endl;

        for (int i=0; i<n_atoms; i++){
			for (int xyz=0; xyz<3; xyz++){
				atom[i].set_v(xyz, atom[i].get_v(xyz)*exp(-0.5*zeta*dt));
				atom[i].set_v(xyz, atom[i].get_v(xyz) + atom[i].get_f(xyz)*0.5*dt/atom[i].get_mass());
				atom[i].set_r(xyz, atom[i].get_r(xyz) + atom[i].get_v(xyz)*dt);
				atom[i].set_r(xyz, mod(atom[i].get_r(xyz), box[xyz]));
			}
            atom[i].reset_f();
        }

		zeta = zeta + (2*kin - 3*n_atoms*R*T_ref)*(dt/Q);
		eta = eta + zeta*dt;
		s = exp(eta);

        // calc force at frame step
        pot = 0.0;

        for (int i=0; i<n_atoms; i++){
		#ifdef OPENMP
		#pragma omp parallel for
		#endif /* OPENMP */
            for (int j=0; j<i; j++){
                itp[i][j].set_r_ij(atom, box);
                itp[i][j].set_f_ij(atom, r_c);
                pot += itp[i][j].set_p_ij(atom, r_c);
            }
        }


        kin = 0.0;
        for (int i=0; i<n_atoms; i++){
			for (int xyz=0; xyz<3; xyz++){
				atom[i].set_v(xyz, atom[i].get_v(xyz) + atom[i].get_f(xyz)*0.5*dt/atom[i].get_mass());
				atom[i].set_v(xyz, atom[i].get_v(xyz)*exp(-0.5*zeta*dt));
			}
            kin += atom[i].get_kin();
        }

		// energy calc
		T = (2.0*kin)/(3.0*n_atoms*R);
		vir = 0.0;
        for (int i=0; i<n_atoms; i++){
            for (int j=0; j<i; j++){
				for (int xyz=0; xyz<3 ; xyz++){
					vir += atom[i].get_r(xyz)*atom[i].get_f(xyz);
				}
			}
		}
		P = (2.0*kin + vir)/(3*box[0]*box[1]*box[2]);
		H = kin + pot + 0.5*Q*pow(zeta, 2) + 3*n_atoms*R*T_ref*eta;

		// energy
		ofs_energy << dt*frame << " "
		           << H << " "
		           << pot << " "
		           << T << " "
		           << s << " "
		           << P << " "
		           << vir << " "
		           << " " << endl; 
		//vmd
		ofs_coord << n_atoms << endl;
		ofs_coord << dt*frame << endl;
		for (int i=0; i<n_atoms; i++){
				ofs_coord << atom[i].get_atomname() << " "
						<< atom[i].get_r(0) << " " 
						<< atom[i].get_r(1) << " " 
						<< atom[i].get_r(2) << endl;
		}


    }

    for (int i=0; i<n_atoms;i++){
		gro[i].set_params(atom[i].get_resnr(),
		               atom[i].get_resname(),
					   atom[i].get_atomname(),
					   atom[i].get_num(),
					   atom[i].get_r(0),
					   atom[i].get_r(1),
					   atom[i].get_r(2)
					  );
    }
	write_gro("out.gro", n_atoms, gro, box);

	// end of md


    delete [] gro;
    delete [] atom;
    return 0;
}


