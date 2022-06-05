#include "checkpath.h"

bool checkPath(vector<int> path){

	bool check = false;
	int size = path.size();

	for(int i=0; i<size; i++){
		for(int j=0; j<i; j++){
			if(path[j] == path[i]){
				check = true;
			}
		}
	}
	return check;
}

