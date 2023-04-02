// class to store 3d array in a 1d array with simple acess and config methods
#ifndef mMatrix3_H
#define mMatrix3_H

#include <iostream>
#include <cassert>
#include <fstream>

//declaration
template<class T>
class mMatrix3
{
	public:
		mMatrix3(); //default
		mMatrix3(int x,int y,int z); // construct a 3d array of x by y by z dimensions and init to 0
		mMatrix3(const mMatrix3<T> &arr); //copy constructor
		~mMatrix3(); //destructor to free up memory
		// config/init methods
		void resize(int x,int y,int z);
		//getter methods 
		int getx(); //get sizes of the 3d matrix 
		int gety();
		int getz();
		//access methods for contigency
		T getval(int i, int j, int k);
		void setval(int i, int j, int k, T val);

		// operator overload for () to access elements
		T &operator()(const int i,const int j,const int k);
	private:
		T *m_matrix; //actual array that stores in data 
		int m_i ,m_j,m_k,m_n; // dimensions of 3d matrix
		int getindex(int i, int j, int k); // get index for 1d array given the location of a point in the 3d array
};

//definations
template<class T> //default constructor
mMatrix3<T>::mMatrix3()
{
	m_i = 1;
	m_j = 1;
	m_k = 1;
	m_n = m_i*m_j*m_k; 
	m_matrix = new T[m_n];
	m_matrix[0] = 0.0;
}

template<class T> //empty init constructor
mMatrix3<T>::mMatrix3(int x, int y, int z)
{
	m_i = x;
	m_j = y;
	m_k = z;
	m_n = m_i*m_j*m_k; 
	m_matrix = new T[m_n];
	for(int r=0;r<m_n;r++)
		m_matrix[r] = 0.0;
}

template<class T>
mMatrix3<T>::mMatrix3(const mMatrix3<T> &mtx) //copy constructor
{
	m_i = mtx.m_i;
	m_j = mtx.m_j;
	m_k = mtx.m_k;
	m_n = m_i*m_j*m_k;
	m_matrix = new T[m_n];
	for(int r=0;r<m_n;r++)
		m_matrix[r] = mtx.m_matrix[r];
}

//destructor to free up memory occupied by instance
template<class T>
mMatrix3<T>::~mMatrix3()
{
	if(m_matrix != nullptr)
		delete [] m_matrix;
}

//resize method to init most arrays
template<class T>
void mMatrix3<T>::resize(int x, int y, int z)
{
	m_i = x;
	m_j = y;
	m_k = z;
	m_n= x*y*z;
	delete [] m_matrix; //delete old instance of memory
	m_matrix  = new T[m_n];
	if(m_matrix != nullptr)
	{
		for(int r=0;r<m_n;r++)
		{
			m_matrix[r] = 0.0;
		}
	}
}


// get linear index from 3d index 
template<class T>
int mMatrix3<T>::getindex(int i, int j, int k)
{
	assert(i<m_i && j<m_j && k<m_k and i>=0 and j>=0 and k>=0); //make sure index query is within bounds
	return i*m_i + j*m_j + k;
}


//getter functions for dimensions of array 
template<class T>
int mMatrix3<T>::getx()
{
	return m_i;
}
template<class T>
int mMatrix3<T>::gety()
{
	return m_j;
}

template<class T>
int mMatrix3<T>::getz()
{
	return m_k;
}

template<class T>
T mMatrix3<T>::getval(int i,int j, int k)
{
	assert(i<m_i && j<m_j && k<m_k and i>=0 and j>=0 and k>=0);//make sure index query is within bounds
	int linear_index = getindex(i,j,k);
	return m_matrix[linear_index];
}

template<class T>
void mMatrix3<T>::setval(int i,int j, int k, T val)
{
	assert(i<m_i && j<m_j && k<m_k and i>=0 and j>=0 and k>=0);	int linear_index = getindex(i,j,k);
	m_matrix[linear_index] = val;
}

//overload () operator to return reference to value stored at requested location 
template<class T>
T &mMatrix3<T>::operator()(const int i, const int j, const int k)
{
	assert(i<m_i && j<m_j && k<m_k and i>=0 and j>=0 and k>=0);	int linear_index = getindex(i,j,k);
	return m_matrix[linear_index];
}

#endif 
