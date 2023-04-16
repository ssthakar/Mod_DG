/* A simple template class to store in 2d arrays in flattened format in 1d vectors*/
#ifndef VMATRIX_H
#define VMATRIX_H
#include <initializer_list>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

//TODO implement iterator functionality for more control over access
template <class T>
class vmatrix2
{
public:
  //constructors
  vmatrix2();                                               //++defualt constructor
  vmatrix2(const int nrows, const int ncols);// constrcutor for 2d matrix having nrows and ncols all values inited to zero
  vmatrix2(std::initializer_list<T> input); //list initialization
  vmatrix2(const vmatrix2<T> &vector); //copy constructor
	// methods
  void init(int nrows, int ncols);                       // initialize a vector and set all values to 0 using resize
  void printMatrix(vmatrix2<T> &matrix, std::string s1); // method to print matrix to file
  int rows();
  int cols();
  int size();
  T &operator()(const int row, const int col); // operator overload for acess and assignment
private:
  int index(int row, int col);
  std::vector<T> m_vec;    // main container
  int m_rows, m_cols, m_n; // rows, cols, size of matrix
};

// defualt constrcutor
template <class T>
vmatrix2<T>::vmatrix2()
{
  m_rows = 1;
  m_cols = 1;
  m_n = 1;
  std::vector<T> m_vec(1);
}

// static initialization using constructor
template <class T>
vmatrix2<T>::vmatrix2(const int nrows, const int ncols)
{
  m_rows = nrows;
  m_cols = ncols;
  m_n = m_rows * m_cols;
  m_vec.resize(m_n, 0.0);
}

template<class T>
vmatrix2<T>::vmatrix2(const vmatrix2<T> &vector)
{
	m_vec = vector.m_vec;
}

//list initializer constructor to construct a matrix
template<class T>
vmatrix2<T>::vmatrix2(std::initializer_list<T> input)
{
	m_rows = input.size();
	m_cols = 1;
	m_n = m_rows;
	m_vec.insert(m_vec.begin(),input.begin(),input.end());
}


// function to get striding index for 1d vector from 2d matrix rows and cols
template <class T>
int vmatrix2<T>::index(int row, int col)
{
  return row * m_cols + col;
}

// need this because m_vec is private, size of the 2d Matrix;
template <class T>
int vmatrix2<T>::size()
{
  return m_n;
}

// methods to return rows and cols of a matrix stored in a 1d vector
template <class T>
int vmatrix2<T>::cols()
{
  return m_cols;
}

template <class T>
int vmatrix2<T>::rows()
{
  return m_rows;
}

//operator overload for () for assignment and access of the data
template <class T>
T &vmatrix2<T>::operator()(int row, int col)
{
  assert(row < m_rows && col < m_cols && row >= 0 && col >= 0);
  int lindex = index(row, col);
  return m_vec[lindex];
}

//dynamic initialization
template <class T>
void vmatrix2<T>::init(int nrows, int ncols)
{
  m_rows = nrows;
  m_cols = ncols;
  m_n = nrows * ncols;
  m_vec.resize(m_n, 0.0);
}

//simple function to print out arrays in .dat format for debugging and plotting
template <class T>
void printMatrix(vmatrix2<T> &matrix, std::string s1)
{
  std::ofstream file(s1);
  int nrows = matrix.rows();
  int ncols = matrix.cols();
  // assert(nrows*ncols <64);
  for (int row = 0; row < nrows; row++)
  {
    for (int col = 0; col < ncols; col++)
    {
      file << matrix(row, col) << " ";
    }
    file << std::endl;
  }
}
 
#endif
