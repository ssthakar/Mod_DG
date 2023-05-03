// class to store 3d array in a 1d array with simple acess and config methods
#ifndef mMatrix3_H
#define mMatrix3_H
#include "vmatrix2.h"
#include <iostream>
#include <cassert>
#include <fstream>
#include <vector>
//declaration
template<class T>
class mMatrix3
{
	public:
		mMatrix3(); //default
		mMatrix3(int x,int y,int z); // construct a 3d array of x by y by z dimensions and init to 0
		mMatrix3(const mMatrix3<T> &matrix); //copy constructor
		// config/init methods
		void init(int x,int y,int z);
		//getter methods 
		int x(); //get sizes of the 3d matrix 
		int y();
		int z();
    int size(); // get size of the matrix, mainly used to see if size>1 or resize has happened somewhere in the code
		//access methods for contigency

		// operator overload for () to access elements
		T &operator()(const int i,const int j,const int k);
    void test_1(mMatrix3<T> &mat); //function to test if the indexing is correct or not
	  void reset(); //resets the entire vector to zero
  public:
  	std::vector<T> m_vec; //actual array that stores in data 
		int m_i ,m_j,m_k,m_n; // dimensions of 3d matrix
		int index(int i, int j, int k); // get index for 1d array given the location of a point in the 3d array
};

template<class T>
void mMatrix3<T>::test_1(mMatrix3<T> &mat)
{
  //mat.init(3, 4, 3);
  for(int i=0;i<mat.size();i++)
  {
    //std::cout<<"i"<<std::endl;
    mat.m_vec[i] = i;
  }
  //std::cout<<mat(1,1,1);
}

//definations
template<class T> //default constructor
mMatrix3<T>::mMatrix3()
{
	m_i = 1;
	m_j = 1;
	m_k = 1;
	m_n = m_i*m_j*m_k; 
	std::vector<T> m_vec(1);
}

template<class T> //empty init constructor
mMatrix3<T>::mMatrix3(int x, int y, int z)
{
	m_i = x;
	m_j = y;
	m_k = z;
	m_n = m_i*m_j*m_k; 
	m_vec.resize(m_n,0.0);
}

template<class T>
mMatrix3<T>::mMatrix3(const mMatrix3<T> &matrix)
{
	m_vec = matrix.m_vec;
}


//resize method to init most arrays for dimensions not known at compile time
template<class T>
void mMatrix3<T>::init(int x, int y, int z)
{
	m_i = x;
	m_j = y; //rows of the 2d slice 
	m_k = z; //cols of the 2d slice 
	m_n= x*y*z; //total elements
  m_vec.resize(m_n,0.0);
}


// get linear index from 3d index 
template<class T>
int mMatrix3<T>::index(int i, int j, int k)
{
	assert(i<m_i && j<m_j && k<m_k and i>=0 and j>=0 and k>=0); //make sure index query is within bounds
	return i*m_j*m_k + j*m_k + k;
}


//getter functions for dimensions of array 
template<class T>
int mMatrix3<T>::x()
{
	return m_i;
}
template<class T>
int mMatrix3<T>::y()
{
	return m_j;
}

template<class T>
int mMatrix3<T>::z()
{
	return m_k;
}

template<class T>
int mMatrix3<T>::size()
{
  return m_n;
}

//overload () operator to return reference to value stored at requested location 
template<class T>
T &mMatrix3<T>::operator()(const int i, const int j, const int k)
{
	assert(i<m_i && j<m_j && k<m_k and i>=0 and j>=0 and k>=0);	
  int lindex = index(i,j,k);
	return m_vec[lindex];
}

template<class T>
void mMatrix3<T>::reset()
{
  std::fill(m_vec.begin(),m_vec.end(),0.0);	
}


#endif 
