#ifndef DATA_HDF_BUFFERED_HDF_ARRAY_H_
#define DATA_HDF_BUFFERED_HDF_ARRAY_H_

// HDF5 library includes
#include "hdf5.h"
#include "H5Cpp.h"

// PBHDF5 library includes
#include "HDFConfig.h"
#include "HDFData.h"
#include "HDFGroup.h"
#include "HDFWriteBuffer.h"
#include "HDFFile.h"
#include <assert.h>
#include <iostream>
#include <DNASequence.h>
#include <FASTQSequence.h>
using namespace std;
using namespace H5;

/*
 *
 * Implementation of a 1-D array for IO from an HDF array.
 * This is templated, but specialized for a few data types, so that 
 * the HDF data types do not need to be specified by anybody.
 *
 *  Two examples of the usage of this class follow:
 *
 *	HDFArray<int> nElemArray;
 * 	nElemArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/NumEvent");
 *  nElemArray.Read(i, i+1, &nElem);
 *
 * 	HDFArray<unsigned char> qualArray;
 *	qualArray.Initialize(hdfFile, "PulseData/BaseCalls/QualityValue");
 *  qualArray.Read(cur, cur + nElem, qual);
 *
 */

template<typename T>
class BufferedHDFArray : public HDFData, public HDFWriteBuffer<T> {
 protected:
 public:
	hsize_t   nDims;
	hsize_t   *dimSize;
	int       maxDims;
	UInt      arrayLength;

	/*
	 * Constructor meant to be used for data that will be written.  
	 * This allocates the write buffer.
	 */
 BufferedHDFArray(int pBufferSize=1024) : HDFData() {
		nDims = 0;
    maxDims = 0;
    arrayLength = 0;
		dimSize = NULL;
		this->bufferIndex = 0;
		this->InitializeBuffer(pBufferSize);
	}
	

	BufferedHDFArray(CommonFG* _container, string _datasetName) : HDFData(_container, _datasetName) {
    // no-op
	}

	~BufferedHDFArray() {
		//
		// Clean up the write buffer.
		//
		if (dimSize != NULL) {
			delete[] dimSize;
      dimSize = NULL;
		}
		this->HDFWriteBuffer<T>::~HDFWriteBuffer<T>();
	}
  
  void SetBufferSize(int _bufferSize) {
		this->InitializeBuffer(_bufferSize);
  }

  /*
	void Initialize(HDFGroup parentGroup, string _datasetName, int _bufferSize = 0) {
		container   = &parentGroup.group;
		datasetName = _datasetName;
    //
    // Try and initialize the dataset if it exists.
    //
		this->InitializeBuffer(_bufferSize);
	}
  */


    /*	
 	void Write(const T *data, int dataLength) {
		// Fill the buffer with data. When there is overflow, write
		// that out to disk.
		//
		int dataIndex = 0;
		int bufferCapacity;
		int bufferFillSize = 0;
		bool flushBuffer;
		while(dataIndex < dataLength) {
			bufferCapacity = this->bufferSize - this->bufferIndex;
			flushBuffer = false;
			if (bufferCapacity  > dataLength - dataIndex) {
				bufferFillSize = dataLength - dataIndex;
			}
			else {
				bufferFillSize = bufferCapacity;
				flushBuffer = true;
			}
			memcpy((void*) &this->writeBuffer[this->bufferIndex], (void*) &data[dataIndex], sizeof(T)*bufferFillSize);
			dataIndex   += bufferFillSize;
			this->bufferIndex += bufferFillSize;
			if (flushBuffer) {
				Flush();
			}
		}
	}
    */
 	
    void Write(const T *data, UInt dataLength, bool append=true, UInt writePos = 0) {
		// Fill the buffer with data. When there is overflow, write
		// that out to disk.
		//
		UInt dataIndex = 0;
		int bufferCapacity;
		int bufferFillSize = 0;
		bool flushBuffer;
		while(dataIndex < dataLength) {
			bufferCapacity = this->bufferSize - this->bufferIndex;
			flushBuffer = false;
			if (bufferCapacity  > dataLength - dataIndex) {
				bufferFillSize = dataLength - dataIndex;
			}
			else {
				bufferFillSize = bufferCapacity;
				flushBuffer = true;
			}
			memcpy((void*) &this->writeBuffer[this->bufferIndex], (void*) &data[dataIndex], sizeof(T)*bufferFillSize);
			dataIndex   += bufferFillSize;
			this->bufferIndex += bufferFillSize;
			if (flushBuffer) {
				Flush(append, writePos);
			}
		}
	}

	void Flush(bool append=true, UInt writePos = 0) {
		//
		// Flush contents of current buffer to the file.
		//
		if (this->WriteBufferEmpty()) {
			// 
			// There is no data in the buffer, so nothing can be written.
			// HDF does not support empty arrays (as far as I can tell), so
			// don't even bother trying to create the dataspace.
			// 
			return;
		}
    
		// fetch the current size of the dataspace
		if (fileDataSpaceInitialized == false) {
      cout << "ERROR, trying to flush a dataset that has not been created or initialized" << endl;
      exit(1);
			fileDataSpaceInitialized = true;
		}

    DataSpace fileSpace;
    fileSpace = dataset.getSpace();

    //
    // Load the current size of the array on disk.
    //
    hsize_t fileArraySize[1], blockStart;
    fileArraySize[0] = fileSpace.getSimpleExtentNpoints();
    if (append) {
      blockStart = fileSpace.getSimpleExtentNpoints();
      fileArraySize[0] += this->bufferIndex;
      //
      // Make room in the file for the array.
      //
      dataset.extend(fileArraySize);
    }
    else {
      blockStart = writePos;
      if (blockStart + this->bufferIndex > fileArraySize[0]) {
        fileArraySize[0] = blockStart + this->bufferIndex;
        dataset.extend(fileArraySize);
      }
    }
		
    DataSpace extendedSpace = dataset.getSpace();
    int extendedSize        = extendedSpace.getSimpleExtentNpoints();
    //
    // Configure the proper addressing to append to the array.
    //
    hsize_t dataSize[1];
    hsize_t offset[1];
    dataSize[0] = this->bufferIndex;
    offset[0]   = blockStart;
    extendedSpace.selectHyperslab(H5S_SELECT_SET, dataSize, offset);
    DataSpace memorySpace(1, dataSize);

    //
    // Finally, write out the data.  
    // This uses a generic function which is specialized with
    // templates later on to t
    // memorySpace addresses the entire array in linear format
    // fileSpace addresses the last dataLength blocks of dataset.
    //
    try {
      TypedWrite(this->writeBuffer, memorySpace, extendedSpace);
    }
    catch(DataSetIException e) {
      cout <<"ERROR! Could not write HDF5 data." << endl;
      e.printError();
      exit(1);
    }
    memorySpace.close();
    extendedSpace.close();
    fileSpace.close();

		// Clear the buffer.
		this->ResetWriteBuffer();
	}

	void TypedWrite(const char **data, const DataSpace &memorySpace, const DataSpace &extendedSpace) {
		StrType varStrType(0,H5T_VARIABLE);
		dataset.write(data, varStrType, memorySpace, extendedSpace);
	}

	void TypedWrite(const T*data, const DataSpace &memorySpace, const DataSpace &extendedSpace) {
		assert("Calling TypedWrite on an unsupported type" == 0);
	}

	void TypedCreate(DataSpace &fileSpace, DSetCreatPropList &cparms) {
		cout << "DEFAULT typed create " << endl;
	}
  

	void Create(HDFGroup &parentGroup, string _datasetName) {
    return Create(&parentGroup.group, _datasetName);
  }


  void Create(CommonFG* _container, string _datasetName) {
    //
    // Initialize where the dataset will go.
    container   = _container;
    datasetName = _datasetName;

		hsize_t dataSize[]    = {0};
		hsize_t maxDataSize[] = {H5S_UNLIMITED};
		DataSpace fileSpace(1, dataSize, maxDataSize);
		DSetCreatPropList cparms;

		/*
		 * For some reason, chunking must be enabled when creating a dataset
		 * that  has an unlimited dimension.  Of course, this is not
		 * mentioned in the hdf5 c++ documentation, because that
		 * docuemntation was written for people who enjoy learning how to
		 * use an API by reading comments in source code.
		 */
		hsize_t      chunk_dims[1] = { 16384 };
		cparms.setChunk( 1, chunk_dims );
		TypedCreate(fileSpace, cparms);
		
		//
		// Since TypedCreate created an assigned a dataset, this array is
		// now initialized.  Do the bookkeeping here.
		//

		isInitialized = true;
    fileDataSpaceInitialized = true;
		fileSpace.close();
	}
  


	/*
	 * Initialize for reading.
	 *
	 * Open a dataset in an hdf file. Only call this on datasets that
	 * exist, since this currently handles errors with opening datasets
	 * by ungracefully exiting the program. 
	 */

  int InitializeForReading(HDFGroup &parentGroup, const string datasetName) {
    return Initialize(parentGroup, datasetName, false);
  }
  
  int Initialize(HDFGroup &parentGroup, const string &datasetName) {
    return Initialize(parentGroup, datasetName, true);
  }

  int Initialize(HDFGroup &parentGroup, const string &datasetName, 
            bool createIfMissing, UInt newArrayLength=0) {
    //
    // For writing to this dataset, start at the first position in the
    // write buffer.
    //
		this->bufferIndex = 0;
		//
		// It's possible that the group may be asked to initialize this
		// dataset when the dataset does not exist.  Check that here.
		//
    bool parentHasObject = parentGroup.ContainsObject(datasetName);
		if ( parentHasObject and InitializeDataset(parentGroup, datasetName) == 0) {
      // 
      // The parent group already contains this dataset.  Try to
      // initialize this dataset and if it does not exist, flag fail.
      //
      return 0;
		}

    //
    // This is a hack to create in read/write mode.  If the parent
    // does not have the object, try and create it.  The problem with
    // trying to open a dataset in append mode is it will fail if the
    // dataset does not exist.
    //
    if (parentHasObject == false) {
      if (createIfMissing) {
        Create(parentGroup, datasetName);
      }
      else {
        //
        // Trying to open a dataset to read only, but it does not
        // exist.  Bail.
        //
        return 0;
      }
    }
    int ret = UpdateH5Dataspace();
    if (newArrayLength > 0) {
        ret *= Resize(newArrayLength);
    }
    return ret;
  }

  int UpdateH5Dataspace() {
      try {
          dataspace = dataset.getSpace();
      }
      catch(H5::DataSetIException &e) { 
          e.printError();
          return 0;
      }
      maxDims = MAX_DIMS;
      try {
          nDims = dataspace.getSimpleExtentNdims();
          /*
           * Prevent abuse of this class for multidimensional IO.
           */
          if (nDims != 1) {
              cout << "ERROR in HDF format: dataset: " << datasetName << " should be 1-D, but it is not." << endl;
              exit(1);
          }

          /*
           * Load in the size of this dataset, and make a map to the whole thing.
           */
          if (dimSize != NULL) {
              delete[] dimSize;
              dimSize = NULL;
          }
          dimSize = new hsize_t[nDims];
          dataspace.getSimpleExtentDims(dimSize);
          arrayLength = dimSize[0];
          if (dimSize[0] == 0) {
              // DONT create a real dataspace if the size is 0
              // cout << "WARNING, trying to open a zero sized dataspace." << endl;
              dataspace.close();
              return 1;
          }
          fullSourceSpace = DataSpace(1, dimSize);
          dataspace.close();
      }
      catch(Exception &e) {
          e.printError();
          return 0;
      }
      return 1;
  }

  int Resize(UInt newArrayLength) {
      //
      // Resize this dataset. May or may not allocate space in file.
      // May or may not write fill value.
      //
      try{
          DataSpace fileSpace;
          fileSpace = dataset.getSpace();

          hsize_t fileArraySize[1];
          fileArraySize[0] = newArrayLength;
          arrayLength = newArrayLength;
          dataset.extend(fileArraySize);
          fileSpace.close();
      } catch(H5::DataSetIException &e) { 
          e.printError();
          return 0;
      }
      return 1;
  }


  void Close() {
    if (dimSize != NULL) {
      delete[] dimSize;
      dimSize = NULL;
      HDFData::Close();
    }
  }

	UInt size() {
		dataspace = dataset.getSpace();
		hsize_t dimSizeArray[1];
		dataspace.getSimpleExtentDims(dimSizeArray);
		dataspace.close();
		return dimSizeArray[0];
	}

	/*
	 * Unspecialized form of read.
	 * Read cannot be called on a type T* that does not have a
	 * specialized template definition.  This is all determined at
	 * compile time.  To ensure this, the following
	 * default definition is provided that gives a nasty warning and
	 * exits the code.
	 */

  virtual
	void Read(UInt start, UInt end, T* dest) {
		assert("ERROR, calling Read with an unsupported type. Use Read(start,end,datatype, dest) instead." == 0);
		exit(1); // this is in case the assert statement is removed.
	}
  /*  
	virtual void ReadPoly(int start, int end, T* dest) {
    assert("ERROR, calling Read with an unsupported type." == 0);
    exit(1);
  }
  */
  /*
	 * Read in type T from the opened dataset from the interval (start,
	 * end].
	 */
	
	void ReadDataset(vector<T> &dest) {
		assert("ERROR, calling ReadDataset with an unsupported type.");
		exit(1); // this is in case the assert statement is removed.
	}

  
	void Read(UInt start, UInt end, DataType typeID, T *dest) {
		if (end - start == 0) {
			return;
		}
		hsize_t memSpaceSize[] = {0};
		memSpaceSize[0] = end - start;
		hsize_t sourceSpaceOffset[] = {0};
		sourceSpaceOffset[0] = start;
		DataSpace destSpace(1, memSpaceSize);
		fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
		dataset.read(dest, typeID, destSpace, fullSourceSpace);
		destSpace.close();
	}

	void ReadCharArray(UInt start, UInt end, string*dest) {
		hsize_t memSpaceSize[] = {0};
		memSpaceSize[0] = end - start;
		hsize_t sourceSpaceOffset[] = {0};
		sourceSpaceOffset[0] = start;
		DataSpace destSpace(1, memSpaceSize);
		StrType strType(0, H5T_VARIABLE);
		fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
		vector<char*> tmpStringArray;
		tmpStringArray.resize(end-start);
		dataset.read(&tmpStringArray[0], strType, destSpace, fullSourceSpace);
		UInt i;
		for (i = 0; i < tmpStringArray.size(); i++) {
			dest[i] = tmpStringArray[i];
		}
		destSpace.close();
	}
		
};






/*
 * Type specializations for some standard types. Use the macro for
 * vanilla specializations (that only require the HDF type ID to be
 * specified). 
 */
#define DEFINE_TYPED_READ_ARRAY(T, Pred) template<>  \
   	void BufferedHDFArray<T>::Read(UInt start, UInt end, T* dest) { \
   	Read(start,end, Pred, dest); \
	}


DEFINE_TYPED_READ_ARRAY(int, PredType::NATIVE_INT);
DEFINE_TYPED_READ_ARRAY(char, PredType::NATIVE_INT8);
DEFINE_TYPED_READ_ARRAY(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_READ_ARRAY(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_READ_ARRAY(uint16_t, PredType::NATIVE_UINT16);
DEFINE_TYPED_READ_ARRAY(float, PredType::NATIVE_FLOAT);
DEFINE_TYPED_READ_ARRAY(char*, PredType::C_S1);


#define DEFINE_TYPED_READ_DATASET(T, Pred) template<>  \
	void BufferedHDFArray<T>::ReadDataset(vector<T>  &dest) {	 \
	dest.resize(arrayLength); \
  Read(0,arrayLength, Pred, &dest[0]);										\
}

DEFINE_TYPED_READ_DATASET(int, PredType::NATIVE_INT);
DEFINE_TYPED_READ_DATASET(char, PredType::NATIVE_INT8);
DEFINE_TYPED_READ_DATASET(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_READ_DATASET(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_READ_DATASET(uint16_t, PredType::NATIVE_UINT16);
DEFINE_TYPED_READ_DATASET(float, PredType::NATIVE_FLOAT);
DEFINE_TYPED_READ_DATASET(char*, PredType::C_S1);

template<>
void BufferedHDFArray<string>::Read(UInt start, UInt end, string*dest) {
	vector<char*> tmpDestCharPtrs;
	if (end == start) return;
	assert(end > start);
	tmpDestCharPtrs.resize(end-start);
	//	ReadCharArray(start, end, (char**) &tmpDestCharPtrs[0]);
	ReadCharArray(start, end, dest);
	//	unsigned int ptrIndex;
	//	for (ptrIndex = 0; ptrIndex < end - start; ptrIndex++) {
	//		dest[ptrIndex] = tmpDestCharPtrs[ptrIndex];
	//	}
	// Unset the values of the tmp so that they are not destructed when
	// removed from the stack.
	//	tmpDestCharPtrs.resize(0);
}

#define DEFINE_TYPED_CREATE_ARRAY(T,Pred) template<> \
	void BufferedHDFArray<T>::TypedCreate(DataSpace &fileSpace,  DSetCreatPropList &cparms) { \
	T zero; zero = 0;\
	cparms.setFillValue(Pred,&zero);\
	dataset = container->createDataSet(datasetName.c_str(), Pred, fileSpace, cparms); \
}
	//	isInitialized = true;																							\

DEFINE_TYPED_CREATE_ARRAY(int, PredType::NATIVE_INT);
DEFINE_TYPED_CREATE_ARRAY(char, PredType::NATIVE_INT8);
DEFINE_TYPED_CREATE_ARRAY(char*, StrType(0,H5T_VARIABLE));
DEFINE_TYPED_CREATE_ARRAY(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_CREATE_ARRAY(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_CREATE_ARRAY(float, PredType::NATIVE_FLOAT);
DEFINE_TYPED_CREATE_ARRAY(uint16_t, PredType::NATIVE_UINT16);

template<>
void BufferedHDFArray<string>::TypedCreate(DataSpace &space, DSetCreatPropList &cparms) {
	StrType varStrType(0,H5T_VARIABLE);
	dataset = container->createDataSet(datasetName.c_str(), varStrType, space, cparms);
}

#define DEFINE_TYPED_WRITE_ARRAY(T, Pred) template<>													\
	void BufferedHDFArray<T>::TypedWrite(const T *data, const DataSpace &memorySpace, const DataSpace &fileSpace) {	\
		dataset.write(data, Pred, memorySpace, fileSpace);									\
	}


DEFINE_TYPED_WRITE_ARRAY(int, PredType::NATIVE_INT)
DEFINE_TYPED_WRITE_ARRAY(unsigned int, PredType::NATIVE_UINT)
DEFINE_TYPED_WRITE_ARRAY(unsigned char, PredType::NATIVE_UINT8)
DEFINE_TYPED_WRITE_ARRAY(char, PredType::NATIVE_INT8)
DEFINE_TYPED_WRITE_ARRAY(float, PredType::NATIVE_FLOAT)
DEFINE_TYPED_WRITE_ARRAY(uint16_t, PredType::NATIVE_UINT16)
DEFINE_TYPED_WRITE_ARRAY(char*, StrType(0,H5T_VARIABLE))
DEFINE_TYPED_WRITE_ARRAY(string, StrType(0,H5T_VARIABLE))

/*
 * This is a nonstandard definition because it requires the creation
 * of a special datatype for variable length string type.
 */

//
// Use this as the base class for other lists.
//
typedef BufferedHDFArray<void> BaseHDFArray;


#define DEFINE_TYPED_CLASS(CLASS_NAME, TEMPLATE_TYPE) \
class CLASS_NAME : public BufferedHDFArray< TEMPLATE_TYPE  > { \
 public: \
  void Read(UInt start, UInt end, TEMPLATE_TYPE *dest) { \
    BufferedHDFArray<TEMPLATE_TYPE>::Read(start, end, dest); \
  } \
}; 

DEFINE_TYPED_CLASS(HDFIntArray, int)
DEFINE_TYPED_CLASS(HDFUIntArray, unsigned int)
DEFINE_TYPED_CLASS(HDFUCharArray, unsigned char)
DEFINE_TYPED_CLASS(HDFCharArray, char)
DEFINE_TYPED_CLASS(HDFUShortArray, uint16_t )
DEFINE_TYPED_CLASS(HDFStringArray, string)
DEFINE_TYPED_CLASS(HDFFloatArray, float)


#endif
