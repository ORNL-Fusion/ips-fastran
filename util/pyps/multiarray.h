/**
 * $Id: multiarray.h 868 2011-04-10 15:46:01Z cary $
 *
 */

#ifndef MULTIARRAY_H
#define MULTIARRAY_H

// std includes
#include <cstdlib>  // To get size_t
#include <vector>

/**
 * A multi-dimensional array class such that thye data are
 * stored contiguously in memory. 
 */

template <typename TYPE>
class MultiArray {

  public:
/**
 * Constructor 
 */
    MultiArray();

/**
 * Destructor
 */
    virtual ~MultiArray();

/** 
 * Set dimensions
 * @param dms sizes along each dimension
 */
    void setDims(const std::vector<size_t>& dms);

/** 
 * Get dimensions
 * @return sizes along each dimension
 */
    std::vector<size_t> getDims() const;

/** 
 * Get size
 * @return total number of elements
 */
    size_t getSize() const;

/** 
 * Set data
 * @param dta flat array of data
 */
    void setData(const std::vector<TYPE>& dta);

/** 
 * Set data
 * @param dta flat array of data
 */
    const TYPE* getDataPointer() const;

/** 
 * Left hand side global indexer
 * @param global index
 * @return the value of the MultiArray at index
 */
    TYPE& operator[] (const size_t index);

/** 
 * Right hand side global indexer
 * @param global index
 * @return the value of the MultiArray at index
 */
    const TYPE& operator[] (const size_t index) const;

/** 
 * Indexer
 * @param indices an array of indices
 * @return the value of the MultiArray at indices
 */
    TYPE& operator[] (const std::vector<size_t>& indices);

  private:

    size_t getGlobalIndex(const std::vector<size_t>& indices) const;

    std::vector<TYPE> _data;
    std::vector<size_t> _sizes;
    size_t _ntot;
};

#ifdef SWIG
   %template(MultiArrayInt)   MultiArray<int>;
   %template(MultiArrayDouble)   MultiArray<double>;
   %template(MultiArrayChar)   MultiArray<char>;
#endif
#endif // MULTIARRAY_H
