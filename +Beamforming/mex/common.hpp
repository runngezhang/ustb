#define NOMINMAX
#include "mex.h"

inline void* MallocPersistent(mwSize pSz)
{
    void* iPtr = mxMalloc(pSz);
    if(iPtr == NULL)
        mexErrMsgTxt("Error allocating memory");

    mexMakeMemoryPersistent(iPtr);

    return iPtr;
}

inline void FreePersistent(void* pPtr)
{
    if(pPtr != NULL)
        mxFree(pPtr);
}

template<class T>
inline bool GetScalar(const mxArray * pArray,T& pValue)
{
    if(pArray == NULL)
        return false;

    pValue = static_cast<T>(mxGetScalar(pArray));
    return true;
}

inline void Init()
{
    //UxLogger::InitLogger("surfprolab.log",UxLogDEBUG);
}

template<class T,class S>
void CopyData(S* pDataDst,T* pDataSrc,const int pLength)
{
    for(int k=0;k<pLength;k++)
	{
        pDataDst[k] = static_cast<S>(pDataSrc[k]);
	}
}