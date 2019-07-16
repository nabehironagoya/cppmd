#include "read_gro.hpp"

int Gro::cnt = 0;
Gro::Gro(){
    cnt++;
    this->resnr = ' ';
    this->resname = ' ';
    this->atomname = ' ';
    this->num = cnt;
    this->r[0] = 0.0;
    this->r[1] = 0.0;
    this->r[2] = 0.0;
}

Gro::~Gro(){
    cnt--;
}

void Gro::set_params(int resnr, string resname, string atomname, int num, float x, float y, float z){
    this->resnr = resnr;
    this->resname = trim(resname, ' ');
    this->atomname = trim(atomname, ' ');
    this->num = num;
    this->r[0] = x;
    this->r[1] = y;
    this->r[2] = z;
}

int Gro::get_params(){
#ifdef DEBUG
    cout << resname << " " <<atomname << endl;
#endif /* DEBUG */
    return num;
}



int get_n_atoms(string grofile){
    string buf;
    ifstream ifs(grofile);
    getline(ifs, buf);
    getline(ifs, buf);
    int n_atoms = stoi(buf);
    return n_atoms;
}

void read_gro(string grofile, int n_atoms,  Gro *gro, float *box){
    string buf;
    ifstream ifs(grofile);
    getline(ifs, buf);
    getline(ifs, buf);

    for (int i=0; i<n_atoms;i++){
        getline(ifs, buf);
        gro[i].set_params(
                            stoi(buf.substr(0,5)),
                            buf.substr(5,5),
                            buf.substr(10,5),
                            stoi(buf.substr(15,5)),
                            stof(buf.substr(20,8)),
                            stof(buf.substr(28,8)),
                            stof(buf.substr(36,8))
                           );
    }

    getline(ifs, buf);
    box[0] = stof(buf.substr(0,10));
    box[1] = stof(buf.substr(10,10));
    box[2] = stof(buf.substr(20,10));

}

//int main(int argc, char **argv){
//    string gro = argv[1];
//    int n_atoms = get_n_atoms(gro);
//    float *box = new float[3];
//    Gro *gro = new Gro[n_atoms];
//
//    read_gro(gro,n_atoms, gro, box);
//
//    cout << n_atoms << endl;
//    for (int i=0; i<n_atoms;i++){
//        gro[i].get_params();
//    }
//    cout << box[0] ;
//    cout << gro[0].get_r(2);
//
//    delete [] gro;
//    delete [] box;
//    return 0;
//}
