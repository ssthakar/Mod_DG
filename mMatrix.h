#ifndef mMatrix2_H 
#define mMatrix2_H
#include <cassert>
#include <fstream>
#include <iostream>

/*----------------Matrix class to store a 2 Dimensional array in a 1d array and use striding to access 2d elements---------------------------------------*/


// might  overhaul this to use std::array instead of standary  C arrays
//NOTE: To do to consecutive binary operations of the sort of A*B*C, use brackets 
// for eg: (A*B)*C , otherwise assertions will fail leading to an error

//TODO


// * check some assertions, they need to be in tune with c++ 0 index notation
//
// * assertions need to be placed in some methods and operator overloads to make code more rigourous and less prone to user error 
// * config methods
//    **   upper diagonal storage reduced by half
//    **   lower diagonal storage reduced by half

// make matrix matrix operator overloads members instead of friend functions
// two template classes to binary operate different types // really need this but complete code overhaul??
// vectors style iterator functionality
// vector functionality like push_back maybe???

/*------------------------------------------------------------------------------------------------------------------------------------------------------*/


// declaration
template <class T>
class mMatrix2
{
  public:
    // constructors
    mMatrix2();//default constructor
	  mMatrix2(int n); //constructor for standard 1d array;
    mMatrix2(int nrows, int ncols); // constructor for a matrix of user defined size
    mMatrix2(int nrows, int ncols, const T val); //construtor for a diagonal matrix with all diagonal elems =val
    mMatrix2(int nrows, int ncols, const T *arr); // constructor with user input as matrix
    mMatrix2(const mMatrix2<T> &arr); // copy constructor 
    ~mMatrix2(); //destructor to delete allocated memory

    //config methods
	  void resize(int nrows,int ncols);

    //acess methods not really need these with the () operator overload but just in case 
    T getval(int row, int col);
	  void setval(int row, int col, T val);
    int getnrows();
    int getncols();

    //operator overloads
	  T &operator()(const int nrows, const int ncols); //() return by reference overload to access and assign

    bool operator== (const mMatrix2<T> &arr);
    template <class U> friend mMatrix2 <U> operator + (const mMatrix2 <U> &arr1,const mMatrix2 <U> &arr2);
    template <class U> friend mMatrix2 <U> operator + (const U &scalar, const mMatrix2 <U> &arr2);
    template <class U> friend mMatrix2 <U> operator + (const mMatrix2 <U> &arr1,const U &scalar);
   
    template <class U> friend mMatrix2 <U> operator - (const mMatrix2 <U> &arr1,const mMatrix2 <U> &arr2);
    template <class U> friend mMatrix2 <U> operator - (const U &scalar, const mMatrix2 <U> &arr2);
    template <class U> friend mMatrix2 <U> operator - (const mMatrix2 <U> &arr1,const U &scalar);
    
    template <class U> friend mMatrix2 <U> operator * (const mMatrix2 <U> &arr1,const mMatrix2 <U> &arr2);
    template <class U> friend mMatrix2 <U> operator * (const U &scalar, const mMatrix2 <U> &arr2);
    template <class U> friend mMatrix2 <U> operator * (const mMatrix2 <U> &arr1,const U &scalar);
 
  // methods   
		private:
    int getindex(int row, int col); // private method to get index of flat array from matrix format 

  // members 
		private:
    T *m_matrix;
    int m_rows,m_cols,m_n; //rows, columns and number of elements 
};

//definintions

//default constructor, constructs a 1x1 matrix with the value of zero
template <class T>
mMatrix2 <T>::mMatrix2()
{
  m_rows = 1;
  m_cols = 1;
  m_n = m_rows*m_cols;
  m_matrix = new T[m_n];
  m_matrix[0] = 0.0;
}
template<class T>
mMatrix2<T>::mMatrix2(int n)
{
	m_rows = n;
	m_cols = 1;
	m_n = n;
	m_matrix = new T[n];//simple 1d array
	for(int i=0;i<n;i++)
		m_matrix[i] = 0.0;
}

// constructor for a diagonal matrix with all diag elems = val used mainly for identity matrix 
template <class T>
mMatrix2 <T>::mMatrix2(int nrows, int ncols,const T val)
{
  m_rows = nrows;
  m_cols = ncols;
  m_n = nrows*ncols;
  m_matrix = new T[m_n];
  for(int i=0;i<m_n;i++)
  {
    m_matrix[i] = 0.0;
  }
  for(int i =0;i<nrows;i++)
  {
    for(int j=0;j<ncols;j++)
    {
      if(i==j)
      {
        int linearindex = i*m_cols + j;
        m_matrix[linearindex] = val;
      }
    }
  }
}



// constructs a nrows x ncols matrix and inits all vals ot lim_zero or zero 
template <class T>
mMatrix2 <T>::mMatrix2(int nrows, int ncols)
{
  m_rows = nrows;
  m_cols = ncols;
  m_n = m_rows*m_cols;
  m_matrix = new T[m_n];
  for(int i=0;i<m_n;i++)
  {
    m_matrix[i] = 0.0;
  }
}


//construct matrix from user input of flat array of type T
template <class T>
mMatrix2 <T>::mMatrix2(int nrows, int ncols, const T *arr1) // pointer to const type t array 
{
  m_rows = nrows;
  m_cols = ncols;
  m_n = m_rows*m_cols;
  m_matrix = new T[m_n];
  for(int i=0;i<m_n;i++)
  {
    m_matrix[i] = arr1[i];
  }
}

//copy constructor, creates copy of a mMatrix2 class instance
template <class T>
mMatrix2 <T>::mMatrix2(const mMatrix2<T> &mtx1)
{
  m_rows = mtx1.m_rows;
  m_cols = mtx1.m_cols;
  m_n = mtx1.m_n;
  m_matrix = new T[m_n];
  for(int i=0;i<m_n;i++)
  {
    m_matrix[i] = mtx1.m_matrix[i];
  }
}

//destructor to free up memory
template<class T>
mMatrix2<T>::~mMatrix2()
{
  if(m_matrix != nullptr)
    delete[] m_matrix;
}

//resize config method
template<class T>
void mMatrix2<T>::resize(int nrows, int ncols)
{
	m_rows = nrows;
	m_cols = ncols;
	m_n = nrows*ncols;
	delete [] m_matrix; //delete old matrix
	m_matrix = new T[m_n];
	if(m_matrix != nullptr)
	{
		for(int i=0;i<m_n;i++)
		{
			m_matrix[i] = 0.0;
		}
	}
}

//get linear array index from 2d row, col matrix representation
template<class T>
T mMatrix2<T>::getval(int row, int col)
{
  	
  assert(row < m_rows && col < m_cols && row>=0 && col>=0); // note there is no equal to row in the row assertion in order with c++ zero indexing
  int index = getindex(row, col);
  return m_matrix[index];
}

// set value of a particular element given its row and col number 
template<class T>
void mMatrix2<T>::setval(int row, int col, T val)
{
  int index = getindex(row, col);
  assert(index >= 0);
  m_matrix[index] = val;
}

//get number of rows of a matrix
template<class T>
int mMatrix2<T>::getnrows()
{
  return m_rows;
}

//get number of cols of a matrix 
template<class T>
int mMatrix2<T>::getncols()
{
  return m_cols;
}

/*______________________________________________________________________________________________________*/


//matrix + matrix  + operator overload
//

template<class T>
T &mMatrix2<T>::operator ()(int row, int col)
{
	assert(row < m_rows && col < m_cols && row>=0 && col>=0); 
	int linearindex = getindex(row, col);
	return m_matrix[linearindex];
}

template <class T>
mMatrix2<T> operator + (const mMatrix2<T> &arr1 , const mMatrix2<T> &arr2)
{
  int rows1 = arr1.m_rows;
  int rows2 = arr2.m_rows;
  int cols1 = arr1.m_cols;
  int cols2 = arr2.m_cols;
  assert(rows1 == rows2 && cols1 == cols2); //make sure matrices are addable;
  int n = rows1*cols1;
  T *temp = new T[n]; //temp array to store in data for result, use input matrix constructor to generate mMatrix2 later
  for(int i=0;i<n;i++)
    temp[i] = arr1.m_matrix[i] + arr2.m_matrix[i];
  mMatrix2<T> result(rows1, cols1,temp);
  delete[] temp;
  return result;
}


//matrix - matrix  - operator overload
template <class T>
mMatrix2<T> operator - (const mMatrix2<T> &arr1 , const mMatrix2<T> &arr2)
{
  int rows1 = arr1.m_rows;
  int rows2 = arr2.m_rows;
  int cols1 = arr1.m_cols;
  int cols2 = arr2.m_cols;
  assert(rows1 == rows2 && cols1 == cols2); //make sure matrices are addable;
  int n = rows1*cols1;
  T *temp = new T[n]; //temp array to store in data for result, use input matrix constructor to generate mMatrix2 later
  for(int i=0;i<n;i++)
    temp[i] = arr1.m_matrix[i] - arr2.m_matrix[i];
  mMatrix2<T> result(rows1, cols1,temp);
  delete[] temp;
  return result;
}

//scalar plus matrix
template<class T>
mMatrix2<T> operator + (const T &scalar,const mMatrix2<T> &arr1)
{
  int rows = arr1.m_rows;
  int cols = arr1.m_cols;
  int n = rows*cols;
  T *temp = new T[n];
  for(int i=0;i<n;i++)
    temp[i] = scalar +arr1.m_matrix[i];
  mMatrix2<T> result(rows, cols, temp); //input matrix construtor 
  delete[] temp;
  return result;

}


//scalar - matrix
template<class T>
mMatrix2<T> operator - (const T &scalar,const mMatrix2<T> &arr1)
{
  int rows = arr1.m_rows;
  int cols = arr1.m_cols;
  int n = rows*cols;
  T *temp = new T[n];
  for(int i=0;i<n;i++)
    temp[i] = scalar - arr1.m_matrix[i];
  mMatrix2<T> result(rows, cols, temp); //input matrix construtor 
  delete[] temp;
  return result;

}

//matrix plus scalar operator overload 
template<class T>
mMatrix2<T> operator + (const mMatrix2<T> &arr1, const T &scalar)
{
  int rows = arr1.m_rows;
  int cols = arr1.m_cols;
  int n = rows*cols;
  T *temp = new T[n];
  for(int i=0;i<n;i++)
    temp[i] = scalar +arr1.m_matrix[i];
  mMatrix2<T> result(rows, cols, temp); //input matrix construtor 
  delete[] temp;
  return result;
}

//matrix minus scalar operator overload
template<class T>
mMatrix2<T> operator - (const mMatrix2<T> &arr1, const T &scalar)
{
  int rows = arr1.m_rows;
  int cols = arr1.m_cols;
  int n = rows*cols;
  T *temp = new T[n];
  for(int i=0;i<n;i++)
    temp[i] = arr1.m_matrix[i] - scalar;
  mMatrix2<T> result(rows, cols, temp); //input matrix construtor 
  delete[] temp;
  return result;
}

// multiplication overloads

//matrix*matrix operator overload
template<class T>
mMatrix2<T> operator * (const mMatrix2<T> &lhs, const mMatrix2<T> &rhs)
{
  int l_numrows = lhs.m_rows;
  int r_numrows = rhs.m_rows;
  int l_numcols = lhs.m_cols;
  int r_numcols = rhs.m_cols;
  assert(l_numcols == r_numrows); // condition for matrix multiplication to be possible
  T *tempresult = new T[l_numrows*r_numcols]; // array to store in results
  for(int lhsrow=0;lhsrow<l_numrows;lhsrow++)
  {
    for(int rhscol=0;rhscol<r_numcols;rhscol++)
    {
      T elemres = 0.0;
      for(int lhscol=0;lhscol<l_numcols;lhscol++)
      {
        int lhsindex = (lhsrow*l_numcols) + lhscol;
        int rhsindex = (lhscol*r_numcols) + rhscol;
        elemres = lhs.m_matrix[lhsindex]*rhs.m_matrix[rhsindex] +elemres;
      }
      int resultindex = lhsrow*r_numcols+rhscol;
      tempresult[resultindex] = elemres;
    }
  }
  mMatrix2<T> result(l_numrows,r_numcols,tempresult);
  delete[] tempresult;
  return result;
}

//matrix times scalar
template<class T>
mMatrix2<T> operator * (const mMatrix2<T> &arr1, const T &scalar)
{
  int rows = arr1.m_rows;
  int cols = arr1.m_cols;
  int n = rows*cols;
  T *temp = new T[n];
  for(int i=0;i<n;i++)
    temp[i] = arr1.m_matrix[i]*scalar;
  mMatrix2<T> result(rows, cols, temp); //input matrix construtor 
  delete[] temp;
  return result;
}


//scalar times matrix
template<class T>
mMatrix2<T> operator * (const T &scalar,const mMatrix2<T> &arr1)
{
  int rows = arr1.m_rows;
  int cols = arr1.m_cols;
  int n = rows*cols;
  T *temp = new T[n];
  for(int i=0;i<n;i++)
    temp[i] = arr1.m_matrix[i]* scalar;
  mMatrix2<T> result(rows, cols, temp); //input matrix construtor 
  delete[] temp;
  return result;
}



//get linear index given number of rows and cols from 2d matrix 
template<class T>
int mMatrix2<T>::getindex(int row, int col)
{
  //assert(row<m_rows && col<m_cols && row>=0 && col >=0);
  return  row*m_cols+col;
}

// == operator overload 
template<class T> //note this is not a friend function, because its lhs argument is a matrix class itself
bool mMatrix2<T>::operator == (const mMatrix2<T> &arr1)
{
  int rhs_rows = arr1.m_rows;
  int rhs_cols = arr1.m_cols;
  int lhs_rows = this->m_rows;
  int lhs_cols = this->m_cols;
  assert(rhs_rows == lhs_rows && rhs_cols == lhs_cols); //check if the matrix has the same size
  bool flag = true;
  for(int i=0;i<this->m_n;i++)
  {
    if(this->m_matrix[i] != arr1.m_matrix[i])
    {
      flag = false;
    }
  }
  return flag;
}

//small function to print out matrices to file 
//USE FOR UNIT TESTS ONLY!!! Dont use for large matrices!!!!
template <class T>
void printMatrix(mMatrix2<T> matrix,std::string s1)
{
  std::ofstream file(s1);
  int nrows = matrix.getnrows();
  int ncols = matrix.getncols();
  //assert(nrows*ncols <64);
  for(int row=0;row<nrows;row++)
  {
    for(int col=0;col<ncols;col++)
    {
      file<<matrix.getval(row, col)<<" ";

    }
    file<<std::endl;
  }
}


#endif
