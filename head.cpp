
#include <cstdlib>
#include <iostream>
#include <bits/stdc++.h>

using namespace std;



class matrix {
private:
  int nx,ny;
  std::string id;
  double** data;
public:
  matrix (int nx, int ny, std::string id);
  void print ();
};

matrix :: matrix (int nx, int ny, std::string id = "unknown") {
    std::srand(std::time(nullptr));
    
    this -> nx = nx;
    this -> ny = ny;
    this -> id = id;
    //cout << "init matrix " << id << " with size " << nx << "x" << ny << std::endl;
    
    data = new double* [ nx ];
    
    for (int i=0; i < nx; i++) {
      data[i] = new double [ny];
      
      //memset(data[i], 0.01, ny * sizeof(double));
      
      for (int j = 0; j < ny; j++) {
        data[i][j] = double(1 + std::rand()/((RAND_MAX + 1u)/9));
      }
    }
    
};

void matrix :: print () {
    cout << "matrix " << id << " with size " << nx << "x" << ny << std::endl;
    
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        cout << " " << data[i][j];
      }
      cout << std::endl;
    }
};

class vect {
private:
  int size;
  std::string id;
  double* data;
public:
  vect (int size, std::string id);
  void print ();
};

vect :: vect (int size, std::string id = "unknown") {
    this -> size = size;
    this -> id = id;
    cout << "init vect " << id << " with size " << size << std::endl;
    
    data = new double [ size ];
    
    for (int i = 0; i < size; i++) {
      data[i] = 0.0001;
      cout << i << std::endl;
    }
};

void vect :: print () {
    cout << id << " " << size << std::endl;
    
    for (int i = 0; i < size; i++) {
      cout << i << " " << data[i] << std::endl;
    }
};

int main()
{
  //vect v1(10, "test1");
  //vect v2(2 , "test2");
  
  matrix m1(5,8,"testm");
  
  //v1.print();
  //v2.print();
  
  m1.print();
  
  cout << "done" << endl;
}
