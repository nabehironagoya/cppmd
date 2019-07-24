#include <iostream>
#include <string>
#include <vector>
using namespace std;


vector<string> split(string str, char del){
    vector<string> result;
    string temp_str;

    for (char c : str){
        if ( c == del ){
            result.push_back(temp_str);
            temp_str.clear();
        }else{
            temp_str += c;
        }
    }

    result.push_back(temp_str);
    return result;
}

string trim(string str, char del){
    while(str.find(del) < str.size()){
        str.erase(str.begin() + str.find(del));
    }
    return str;
}
