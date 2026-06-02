/** @file Cascade-square-matrix.hh
 *  @brief Definition of an Cascade square matrix
 *
 *  Description:
 *
 *  An Cascade square matrix is a n x n square matrix (where n <= capacity)
 *  which can be extended to the (n+1) x (n+1) matrix (if (n+1) <= capacity)
 *  whenever it is necessary.
 *  
 *  Such a matrix is able to represent a binary relation over a set, a matrix of
 *  counters etc. The data type of the matrix cells is templated.
 *
 *  This file contains an abstract class CascadeSquareMatrix which does
 *  not contain the exact inner representation of the matrix.
 *
 *  The file also contains several concrete subclasses of the 
 *  CascadeSquareMatrix class, namely:
 *  - CascadeSquareMatrix
 *  - DynamicSquareMatrix
 *  - HashedSquareMatrix
 *  which implement the inner representation of the matrix on its own. The
 *  data cells will be accessible exclusively through the methods which are
 *  virtual in context of the abstract class CascadeSquareMatrix and which
 *  are implemented within the subclasses.
 *
 *  Using this class, it is possible to:
 *  - get a value of the matrix cell using two indices in O(1)
 *  - set a value of the matrix cell using two indices in O(1)
 *  - extend a n x n matrix to a (n+1) x (n+1) matrix in O(n)
 *  - implement and use custom inner representation of the matrix
 *  - choose a data type of the matrix cells
 *
 *  @author Tomáš Kocourek
 */

#ifndef _Cascade_SQUARE_MATRIX_HH_
#define _Cascade_SQUARE_MATRIX_HH_

#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <cassert>
#include "mata/utils/sparse-set.hh"
#include <type_traits>

namespace mata::utils {

/************************************************************************
*
*
*                        Cascade SQUARE MATRIX
*
*
*                        (RELATIONS AND COUNTERS)                   
*
*
*************************************************************************/

/**
 * CascadeSquareMatrix
 * 
 * @brief interface for Cascade square matrix implementations
 *
 * Square matrix "n x n" which can be extended to "(n+1) x (n+1)" matrix
 * if n is less than the maximal capacity. Such a class allows us to
 * represent binary relations over carrier set with n elements and adjust
 * it to n+1 elements whenever a new element of the carrier set is created
 * (for example when we represent relation over partition and a block of
 * partition is split in two) or matrices of counters etc.
 * 
 * This abstract class declares methods for accessing elements of the 
 * matrix, assigning to the cells of matrix and extending the matrix by one row
 * and one column (changing the size). 
 * It defines attributes for maximal capacity (which cannot be changed) and
 * current size.
 * It does not define the data structure for storing data. Each subclass
 * which inherits from this abstract class should:
 * - contain the storage for data of datatype T which represents n x n matrix 
 * - implement methods set, get and extend
 * - implement a method clone which creates a deep copy of a matrix
 * Then, the CascadeSquareMatrix can be used independently of the inner
 * representation of the matrix. Therefore, one can dynamically choose from
 * various of implementations depending on the situation. If any new 
 * subclass is implemented, one should also modify the 'create' 
 * function and extend 'MatrixType' enumerator.
 *
 * Note that in context of an n x n matrix, this implementation uses the word
 * 'size' to refer to the number n (number of rows or columns). The word 
 * 'capacity' corresponds to the maximal allowed size (maximal number
 * of rows or columns). 
 *
**/

template <typename T>
class CascadeSquareMatrix {
    protected:
        
        // number of rows (or columns) of the current square matrix
        size_t size_{0};
        
        // maximal allowed number of rows (or columns) of the square matrix        
        size_t capacity_{0};
        
        std::vector<T> data_{};

    public:

        CascadeSquareMatrix(size_t max_rows, size_t init_rows) {
            assert(init_rows <= max_rows && 
                   "Initial size of the matrix cannot be"
                   "bigger than the capacity");
            
            this->capacity_ = max_rows;
            //data_.reserve(this->capacity_ * this->capacity_);
            
            // creating the initial size and filling the data cells with
            // default values
            for(size_t i = 0; i < init_rows; ++i) {extend();}
        }
        
        /** This method provides a way to create a copy of a given
        * CascadeSquareMatrix and preserves the reserved capacity of the vector
        * 'data_'. This goal is achieved using the custom assignment operator.
        * @brief copy constructor of a CascadeSquareMatrix
        * @param other matrix which should be copied
        */
        CascadeSquareMatrix(const CascadeSquareMatrix<T>& other) {
            *this = other;
        }

        
        //
        //  GETTERS
        //
        
        /** Returns a size of the matrix. In this context,
        * the size of an n x n matrix corresponds to the value n.
        * @brief returns the size of the matrix
        * @return size of the matrix
        */
        size_t size(void) const { return size_; }
        
        /** Returns a capacity of the matrix. In this context,
        * the capacity of an n x n matrix corresponds to the value n_max such
        * that if the matrix is extended to the n_max x n_max matrix, it cannot
        * be extended anymore.
        * @brief returns the capacity of the matrix
        * @return capacity of the matrix
        */
        size_t capacity(void) const { return capacity_; }
        
        //
        // VIRTUAL METHODS 
        // 
        // These virtual methods will be implemented in the subclasses according
        // to allow to the concrete representation of the matrix. These methods
        // provide a way to access the contents of the matrix
        //
        
        /**
        * @brief Assigns a value to the cell of the matrix
        * @param[in] i row of the matrix 
        * @param[in] j column of the matrix
        * @param[in] value a value which will be assigned to the memory cell
        */        
        void set(const size_t i, const size_t j, const T value) {
            assert(i < this->size_ && "Nonexisting row cannot be accessed");
            assert(j < this->size_ && "Nonexisting column cannot be accessed");
            
            // accessing the matrix in the cascading way
            data_[i >= j ? i * i + j : j * j + 2 * j - i] = value;
        }

        /**
        * @brief Finds a value of the matrix memory cell
        * @param[in] i row of the matrix 
        * @param[in] j column of the matrix
        * @return a found value in the matrix
        */
        const T get(const size_t i, const size_t j) const {
            assert(i < this->size_ && "Nonexisting row cannot be accessed");
            assert(j < this->size_ && "Nonexisting column cannot be accessed");

            // accessing the matrix in the cascading way
            return data_[i >= j ? i * i + j : j * j + 2 * j - i];
        }
        
        /**
        * @brief changes the n x n matrix to the (n+1) x (n+1) matrix 
        * @param[in] placeholde value which will be assigned to the 
        * newly allocated memory cells 
        */
        void extend(T placeholder = T()) {
            assert(this->size_ < this->capacity_ 
                   && "The matrix cannot be extended anymore");

            // allocation of 2 * size + 1 new data cells
            data_.insert(data_.end(), 2 * this->size_ + 1, placeholder);
            // the size increases
            ++this->size_;
        }

        /**
        * Changes the n x n matrix to the (n+1) x (n+1) matrix by duplicating
        * existing row and column. The row parameter is an index of the row
        * which should be duplicated and added as a (n+1)th row, while the
        * col parameter is an index of the column which should be duplicated
        * and added as an (n+1)th row. If the row parameter equals n, then
        * it will be initialized using default values of the type T. If the
        * col parameter equals n, the new column will be also initialized
        * with the defauls values of the type T. Using this approach, one is
        * able to copy only a row or column and initialize the other one with
        * default values. Calling extend_and_copy(n, n) has the same effect as
        * calling extend(). 
        * The element at the position [n, n] will be always initialized using
        * the default value of the type T. 
        * @brief changes the n x n matrix to the (n+1) x (n+1) matrix with
        * copying of the existing data 
        * @param[in] placeholde value which will be assigned to the 
        * newly allocated memory cells 
        */        
        void extend_and_copy(size_t row, size_t col) {
            assert(this->size_ < this->capacity_ 
                   && "The matrix cannot be extended anymore");
            assert(row <= this->size_
                   && "Index of the copied row cannot be bigger than the size");
            assert(col <= this->size_
                   && "Index of the copied col cannot be bigger than the size");

            // if the row index corresponds to the index of the newly created
            // row, it will be initialized using default values of the type T
            if(row == this->size_) {
                for(size_t i = 0; i < this->size_; ++i) {
                    data_.push_back(T());
                }
            // otherwise, the row with the index 'row' will be duplicated
            // to create a new row
            } else {
                for(size_t i = 0; i < this->size_; ++i) {
                    data_.push_back(get(row, i));
                }
            }
            
            // element at the position [n, n] will be always initialized using
            // the default value of the type T
            data_.push_back(T());
            
            // if the column index corresponds to the index of the newly created
            // column, it will be initialized using default values of the type T
            if(col == this->size_) {
                for(size_t i = 0; i < this->size_; ++i) {
                    data_.push_back(T());
                }            
            // otherwise, the column with the index 'col' will be duplicated
            // to create a new column
            } else {
                for(size_t i = 0; i < this->size_; ++i) {
                    data_.push_back(get(this->size_ - 1 - i, col));
                }
            }
            
            // the size of the matrix increases after extending
            ++this->size_;
        }
        
        CascadeSquareMatrix<T>& operator=(const CascadeSquareMatrix<T>& other) {
            // initialization of the matrix
            this->capacity_ = other.capacity();
            this->size_ = 0;
            this->data_ = std::vector<T>();
            //this->data_.reserve(this->capacity_ * this->capacity_);
            size_t other_size = other.size();
            for(size_t i = 0; i < other_size; ++i) {this->extend();}
            // copying memory cells
            for(size_t i = 0; i < this->size_; ++i) {
                for(size_t j = 0; j < this->size_; ++j) {
                    this->set(i, j, other.get(i, j));
                }
            }
            return *this;
        }
        
        //
        //  MATRIX PROPERTIES
        //
        
        /** This function checks whether the matrix is reflexive. In this
        * context, the matrix is reflexive iff none of the elements on the main
        * diagonal is the zero element of the type T
        * @brief checks whether the Cascade square matrix is reflexive
        * @return true iff the matrix is reflexive
        */
        bool is_reflexive(void) const {
            size_t size = this->size();
            for(size_t i = 0; i < size; ++i) {
                if(!get(i, i)) { return false; }
            }
            return true;
        }

        /** This function checks whether the matrix is antisymetric. In this
        * context, the matrix is antisymetric iff there are no indices i, j
        * where i != j and both matrix[i][j], matrix[j][i] contain nonzero 
        * elementes of the type T
        * @brief checks whether the Cascade square matrix is antisymetric
        * @return true iff the matrix is antisymetric
        */
        bool is_antisymetric(void) const {
            size_t size = this->size();
            for(size_t i = 0; i < size; ++i) {
                for(size_t j = 0; j < size; ++j) {
                    if(i == j) [[unlikely]] { continue; }
                    if(get(i, j) && get(j, i)) { return false; }
                }
            }
            return true;
        }

        /** This function checks whether the matrix is transitive. In this
        * context, the matrix is transitive iff it holds that the input matrix
        * casted to the matrix of booleans (false for zero values of type T, 
        * otherwise true) remains the same if it is multiplied by itself.
        * @brief checks whether the Cascade square matrix is transitive
        * @return true iff the matrix is transitive
        */
        bool is_transitive(void) const {
            size_t size = this->size();
            for(size_t i = 0; i < size; ++i) {
                for(size_t j = 0; j < size; ++j) {
                    bool found = false;
                    for(size_t k = 0; k < size; ++k) {
                        if(get(i, k) && get(k, j)) { 
                            found = true;
                            break; 
                        }
                    }
                    if(!found == static_cast<bool>(get(i, j))) { return false; }
                }
            }
            return true;
        }

        /** Equality operator for the Cascade square matrices. Two matrices
        * are considered to be equal iff they have the same size (both are
        * n x n matrices) and they contain the same values in all the cells.
        * In context of equality, the type and the capacity of matrices
        * do not matter.   
        * @brief checks whether two Cascade square matrices are equal
        * @param[in] other Cascade square matrix
        * @return true iff both matrices are equal
        */        
        bool operator==(const CascadeSquareMatrix<T>& other) const {
            if(size_ != other.size_) { return false; }
            for(size_t row = 0; row < size_; ++row) {
                for(size_t col = 0; col < size_; ++col) {
                    if(get(row, col) != other.get(row, col)) { return false; }
                }
            }
            return true;
        }
        
        CascadeSquareMatrix<T>& operator=(CascadeSquareMatrix<T>&& other) noexcept {
            if (this != &other) {
                this->capacity_ = other.capacity_;
                this->size_ = other.size_;
                this->data_ = std::move(other.data_);
            }
            return *this;
        }
        
                
}; // CascadeSquareMatrix

/**
* @brief debugging function which allows us to print text representation of
* the Cascade square matrix
* @param[out] output stream
* @param[in] maxtrix which will be printed
* @return output stream
*/
template <typename T>
inline std::ostream& operator<<(std::ostream& os, 
    const CascadeSquareMatrix<T>& matrix) {
    
    size_t size = matrix.size();
    size_t capacity = matrix.capacity();    
    std::string result = "\nSIZE: " + std::to_string(size) + "\n";
    result += "CAPACITY: " + std::to_string(capacity) + "\n";
    result += "MATRIX:\n";
    for(size_t i = 0; i < size; ++i) {
        for(size_t j = 0; j < size; ++j) {
            result += std::to_string(
                static_cast<unsigned long>(matrix.get(i, j))) + " ";
        }
        result += "\n";
    }
    return os << result;
}

}

#endif
