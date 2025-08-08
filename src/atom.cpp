/*
 *  read_dcd : c++ class + main file example for reading a CHARMM dcd file
 *  Copyright (C) 2013  Florent Hedin
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//***************** Fully contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#include <cstdlib>
#include <iostream>
#include <vector>
#include <math.h>
#include <cstring>
#include <sstream>
#include "atom.hpp"


using namespace std;


int ATOM::string_comp(vector<string> str) {
    int flag = 0;
    
    if (str[0] == "atom_index") {
        flag = string_comp_atom_index(str);
    } else if (str[0] == "segid") {
        flag = string_comp_segid(str);
    } else if (str[0] == "segname") {
        flag = string_comp_segname(str);
    } else if (str[0] == "resid") {
        flag = string_comp_resid(str);
    } else if (str[0] == "resname") {
        flag = string_comp_resname(str);
    } else if (str[0] == "atomname") {
        flag = string_comp_atomname(str);
    } else if (str[0] == "atomtype") {
        flag = string_comp_atomtype(str);
    } else if (str[0] == "charge") {
        flag = string_comp_charge(str);
    } else if (str[0] == "mass") {
        flag = string_comp_mass(str);
    } else if (str[0] == "beta") {
        flag = string_comp_beta(str);
    } else {
        cout << "Please check your selection argument." << endl;
    }
    return flag;
}

int ATOM::string_comp_atom_index(vector<string> str) {

    int flag = 0;

    if (str[0] == "atom_index") {
        if (str[1] == "<") {
	    if (atom_index < stoi(str[2])) flag = 1;
        } else if (str[1] == "<>") {
	    if (atom_index > stoi(str[2]) && atom_index < stoi(str[3])) flag = 1;
        } else if (str[1] == ">") {
	    if (atom_index > stoi(str[2])) flag = 1;
  	} else { 
            for (unsigned i = 1; i < str.size(); i++) {
                if (atom_index == stoi(str[i])) flag = 1;
     	    }
        }
    } else {
        error_exit( "Error calling string_comp_atom_index()!! ");
    }
    return flag;
}


int ATOM::string_comp_segid(vector<string> str) {

    int flag = 0;

    if (str[0] == "segid") {
        if (str[1] == "<") {
	    if (segid < stoi(str[2])) flag = 1;
        } else if (str[1] == "<>") {
	    if (segid > stoi(str[2]) && segid < stoi(str[3])) flag = 1;
        } else if (str[1] == ">") {
	    if (segid > stoi(str[2])) flag = 1;
  	} else { 
            for (unsigned i = 1; i < str.size(); i++) {
                if (segid== stoi(str[i])) flag = 1;
     	    }
        }
    } else {
        error_exit( "Error calling string_comp_segid()!! ");
    }
    return flag;
}


int ATOM::string_comp_resid(vector<string> str) {

    int flag = 0;

    if (str[0] == "resid") {
        if (str[1] == "<") {
	    if (resid < stoi(str[2])) flag = 1;
        } else if (str[1] == "<>") {
	    if (resid > stoi(str[2]) && resid < stoi(str[3])) flag = 1;
        } else if (str[1] == ">") {
	    if (resid > stoi(str[2])) flag = 1;
  	} else { 
            for (unsigned i = 1; i < str.size(); i++) {
                if (resid== stoi(str[i])) flag = 1;
     	    }
        }
    } else {
        error_exit( "Error calling string_comp_resid()!! ");
    }
    return flag;
}


int ATOM::string_comp_charge(vector<string> str) {

    int flag = 0;

    if (str[0] == "charge") {
        if (str[1] == "<") {
	    if (charge < stof(str[2])) flag = 1;
        } else if (str[1] == "<>") {
	    if (charge > stof(str[2]) && charge < stof(str[3])) flag = 1;
        } else if (str[1] == ">") {
	    if (charge > stof(str[2])) flag = 1;
  	} else { 
            for (unsigned i = 1; i < str.size(); i++) {
                if (charge == stof(str[i])) flag = 1;
     	    }
        }
    } else {
        error_exit("Error calling string_comp_charge()!! ");
    }
    return flag;
}


int ATOM::string_comp_mass(vector<string> str) {

    int flag = 0;

    if (str[0] == "mass") {
        if (str[1] == "<") {
	    if (mass < stof(str[2])) flag = 1;
        } else if (str[1] == "<>") {
	    if (mass > stof(str[2]) && mass < stof(str[3])) flag = 1;
        } else if (str[1] == ">") {
	    if (mass > stof(str[2])) flag = 1;
  	} else { 
            for (unsigned i = 1; i < str.size(); i++) {
                if (mass == stof(str[i])) flag = 1;
     	    }
        }
    } else {
        error_exit( "Error calling string_comp_mass()!! ");
    }
    return flag;
}



int ATOM::string_comp_resname(vector<string> str) {

    int flag = 0;

    //cout << "MSe: the atomselect name is: " << str[1] << "; the actual name is: " << resname << endl;

    if (str[0] == "resname") {
        for (unsigned i = 1; i < str.size(); i++) {
            if (resname == str[i]) flag = 1;
     	}
    } else {
        error_exit ("Error calling string_comp_resname()!! ");
    }
    return flag;
}


int ATOM::string_comp_segname(vector<string> str) {

    int flag = 0;


    if (str[0] == "segname") {
        for (unsigned i = 1; i < str.size(); i++) {
            if (segname == str[i]) flag = 1;
      }
    } else {
        error_exit ("Error calling string_comp_segname()!! ");
    }
    return flag;
}

int ATOM::string_comp_atomname(vector<string> str) {

    int flag = 0;

    if (str[0] == "atomname") {
        for (unsigned i = 1; i < str.size(); i++) {
            if (atomname == str[i]) flag = 1;
     	}
    } else {
        error_exit( "Error calling string_comp_atomname()!! ");
    }
    return flag;
}


int ATOM::string_comp_atomtype(vector<string> str) {

    int flag = 0;

    if (str[0] == "atomtype") {
        for (unsigned i = 1; i < str.size(); i++) {
            if (atomtype == str[i]) flag = 1;
     	}
    } else {
        error_exit( "Error calling string_comp_atomtype()!! ");
    }
    return flag;
}


int ATOM::string_comp_beta(vector<string> str) {

    int flag = 0;

    if (str[0] == "beta") {
        if (str[1] == "<") {
	    if (beta < stoi(str[2])) flag = 1;
        } else if (str[1] == "<>") {
	    if (beta > stoi(str[2]) && beta < stoi(str[3])) flag = 1;
        } else if (str[1] == ">") {
	    if (beta > stoi(str[2])) flag = 1;
  	} else { 
            for (unsigned i = 1; i < str.size(); i++) {
                if (beta == stoi(str[i])) flag = 1;
     	    }
        }
    } else {
        error_exit( "Error calling string_comp_beta()!! ");
    }
    return flag;
}


bool ATOM::is_selection(string str) {
    bool isop = false;
        if (str == "atom_index" || str == "segid" || str == "resid" || str == "resname" || str == "atomname" || str == "atomtype" || str == "charge" || str == "mass" || str == "beta") isop = true;
    return isop;
}

bool ATOM::is_operator(string str) {
    bool isop = false;
        if (str == "n" || str == "*" || str == "+" || str == "(" ) isop = true;
    return isop;
}



int ATOM::select_atoms(vector<string> str) {
    int flag = 0;
    int parenthesis_count = 0;
//    istringstream iss(str);
    string temp1;
    vector<string> arg1;
    vector<string> arithmetic_analog;
        //while (iss >> temp) {
    

    
        for (auto &temp : str) {
// Convert the atom selection to an equation. For ex., "not resid a aa aaa and segname sss ss or atomname bb" is (1 - x1) * x2 + x3
            if (temp == "all") {
		flag = 1;
	        return flag;
		break;
	    }
            if (temp == "(") {
	        arithmetic_analog.push_back(temp);
		parenthesis_count++;
            } else if (temp == "not") {
	        arithmetic_analog.push_back("n");
            } else if (temp == "and") {
              //   if (arithmetic_analog.size() == 0) error_exit("Error: 0 array size for and!!!!!!!!!!"); //debug
                 if (arithmetic_analog.size() == 0 || (arithmetic_analog.size() > 0 && arithmetic_analog[arithmetic_analog.size() - 1] != ")")) {//bug
                     if (arg1.size() == 0) {
    		         error_exit( "Error: No selection before and");
			
                     } else {
  			 if (!is_selection(arg1[0])) error_exit("Error: not a selection option!! Check: need a space between ( and other arguments.");
                         int val = string_comp(arg1);
                         arithmetic_analog.push_back(to_string(val));
 	     	         arg1.clear();
		     }
                 }
                 arithmetic_analog.push_back("*");
            } else if (temp == "or") {
                 //if (arithmetic_analog.size() == 0) error_exit("Error: 0 array size for or!!!!!!!!!!"); //debug
                 if (arithmetic_analog.size() == 0 || (arithmetic_analog.size() > 0 && arithmetic_analog[arithmetic_analog.size() - 1] != ")")) {//bug
                     if (arg1.size() == 0) {
    		         error_exit( "Error: No selection before or");
                     } else {
  			 if (!is_selection(arg1[0])) error_exit("Error: not a selection option!! Check: need a space between ( and other arguments.");
                         int val = string_comp(arg1);
                         arithmetic_analog.push_back(to_string(val));
 	     	         arg1.clear();
		     }
                 }
                 arithmetic_analog.push_back("+");
            } else if (temp == ")") {
                 if (arithmetic_analog.size() == 0 || (arithmetic_analog.size() > 0 && arithmetic_analog[arithmetic_analog.size() - 1] != ")")) {//bug
  	     	     if (!is_selection(arg1[0])) error_exit("Error: not a selection option!! Check: need a space between ( and other arguments. ");
                     int val = string_comp(arg1);
                     arithmetic_analog.push_back(to_string(val));
 	     	     arg1.clear();
                }
	        arithmetic_analog.push_back(temp);
		parenthesis_count--;
            }
              else {
              arg1.push_back(temp);
            }
            temp1 = temp;
    }



    if (temp1 != ")") {
        if (!is_selection(arg1[0])) error_exit("Error: not a selection option!! Check: need a space between ( and other arguments.");
        int val = string_comp(arg1);
        arithmetic_analog.push_back(to_string(val));
 	arg1.clear();
    }

    if (parenthesis_count != 0) error_exit( "Error: parenthesis unbalanced!");

   
    if ( debug_flag > 0 ) {
        cout << "Selection string: ";
	for (auto &string_individual : str) {
            cout << string_individual << " ";
        }
        cout << endl; 
        cout << "The selection formula is: " << endl;
    }


    flag = arithmetic_calc(arithmetic_analog);

    if ( debug_flag > 0 ) {
        for (unsigned i = 0; i < arithmetic_analog.size(); i++) {
            cout << arithmetic_analog[i] << " ";
        }
    
        cout << "= " << flag << endl;
    }
 
    return flag;

}

int ATOM::arithmetic_calc(vector<string> str) {
    int flag = 0;
    stack<string> ops;
    stack<int> vals;
    int val = 0;



    for (unsigned i = 0; i < str.size(); i++) {
        if ( str[i] == "(" ) {
            ops.push(str[i]);
        } else if ( ( ops.size() > 0 ) && ( str[i] == "n" || str[i] == "*" || str[i] == "+" )) {
            while ( ops.size() != 0 && has_precedence( str[i], ops.top() ) ) {
                if ( ops.top() == "n" ) {
		    val = 1 - vals.top();   
		    vals.pop();
		    vals.push(val);
		} else if ( ops.top() == "*" || ops.top() == "+" ) {
		    int val1 = vals.top();
		    vals.pop();
		    int val2 = vals.top();
		    vals.pop();
	            val = op_apply(ops.top(), val1, val2);
		    vals.push(val);
	        }
	        ops.pop();	
            }
            ops.push(str[i]);
        } else if ( str[i] == ")" ) {
   	    while ( ops.top()!= "(") {
                if ( ops.top() == "n" ) {
		    val = 1 - vals.top();   
		    vals.pop();
		    vals.push(val);
		} else if ( ops.top() == "*" || ops.top() == "+" ) {
		    int val1 = vals.top();
		    vals.pop();
		    int val2 = vals.top();
		    vals.pop();
	            val = op_apply(ops.top(), val1, val2);
		    vals.push(val);
	        }	
	     ops.pop();
	  } 
	  ops.pop();
 	    //this is the break point. Note in the end the parenthesis are all popped out.
        } else if (isdigit(str[i][0])) {
            int val_temp = stoi( str[i] );
  	    vals.push(val_temp);
        } 
    }
    //in the end the last operator need to be taken care of
    while ( ops.size() != 0 ) { 
        if ( ops.top() == "n" ) {
	    val = 1 - vals.top();   
	    vals.pop();
	    vals.push(val);
	} else if ( ops.top() == "*" || ops.top() == "+" ) {
	    int val1 = vals.top();
	    vals.pop();
	    int val2 = vals.top();
	    vals.pop();
	    val = op_apply(ops.top(), val1, val2);
	    vals.push(val);
	 }
	 ops.pop();	
    }
    flag = vals.top();
    return flag;
}

bool ATOM::has_precedence(string op1, string op2) { 
    bool has_prec = true;
    unordered_map<string,int> operator_order;
    operator_order["+"] = 0;
    operator_order["*"] = 1;
    operator_order["n"] = 2;
    operator_order["("] = -1;
    if ( operator_order[op1] > operator_order[op2]) has_prec = false;
    return has_prec;
}

int ATOM::op_apply(string op1, int val1, int val2) {
    int op_val = 0;
    if ( op1 == "*" ) {
        op_val = val1 * val2;
    } else if ( op1 == "+" ) {
        op_val = val1 + val2;
    }
    return op_val;
}


void ATOM::error_exit(string str) {
    cout << str << endl;
    exit(1);
}

ATOM::~ATOM()
{
}
