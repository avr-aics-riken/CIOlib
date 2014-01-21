#ifndef _CIO_DFI_INLINE_H_
#define _CIO_DFI_INLINE_H_

/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI.h
 * @brief  cio_DFI Class Header
 * @author kero    
 */

#ifdef CIO_INLINE
 #undef CIO_INLINE
#endif

#ifndef CIO_NO_INLINE
 #define CIO_INLINE inline
#else
 #define CIO_INLINE
#endif

// #################################################################
// フィールドデータの読込み(読み込んだデータのポインタを戻り値として
// 返す）

//template<class T, class TimeT, class TimeAvrT> 
//CIO_INLINE T*
template<class TimeT, class TimeAvrT> 
CIO_INLINE void*
cio_DFI::ReadData(CIO::E_CIO_ERRORCODE &ret,
                  const unsigned step, 
                  const int gc, 
                  const int Gvoxel[3], 
                  const int Gdivision[3], 
                  const int head[3], 
                  const int tail[3],
                  TimeT &time,
                  const bool mode, 
                  unsigned &step_avr, 
                  TimeAvrT &time_avr)
{

   int sz[3];
   for(int i=0; i<3; i++) sz[i]=tail[i]-head[i]+1;
   cio_Array *data = cio_Array::instanceArray
                     ( DFI_Finfo.DataType
                     , DFI_Finfo.ArrayShape
                     , sz
                     , gc
                     , DFI_Finfo.Component);

   double d_time = (double)time;
   double d_time_avr = (double)time_avr;

//   int ret = ReadData(data, step, gc, Gvoxel, Gdivision, head, tail,
   ret = ReadData(data, step, gc, Gvoxel, Gdivision, head, tail,
                       d_time, mode, step_avr, d_time_avr);

   if( ret != CIO::E_CIO_SUCCESS ) {
     delete data;
     return NULL;
   }

//   T* ptr = (T*)data->getData(true);
   void* ptr = data->getData(true);
   delete data;
   time = d_time;
   time_avr = d_time_avr;

   return ptr;
}

// #################################################################
// フィールドデータの読込み(引数で渡された配列にデータを読込む）
template<class T, class TimeT, class TimeAvrT> 
CIO_INLINE
CIO::E_CIO_ERRORCODE cio_DFI::ReadData(T *val,
                                       const unsigned step,
                                       const int gc,
                                       const int Gvoxel[3],
                                       const int Gdivision[3],
                                       const int head[3],
                                       const int tail[3],
                                       TimeT &time,
                                       const bool mode,
                                       unsigned &step_avr,
                                       TimeAvrT &time_avr)
{

   int sz[3];
   for(int i=0; i<3; i++) sz[i]=tail[i]-head[i]+1;

   cio_Array *data = cio_Array::instanceArray
                     ( val
                     , DFI_Finfo.ArrayShape
                     , sz
                     , gc
                     , DFI_Finfo.Component);

   double d_time = (double)time;
   double d_time_avr = (double)time_avr;

   CIO::E_CIO_ERRORCODE ret;
   ret = ReadData(data, step, gc, Gvoxel, Gdivision, head, tail,
                  d_time, mode, step_avr, d_time_avr);

   if( ret == CIO::E_CIO_SUCCESS ) {
     time = d_time;
     time_avr = d_time_avr;
   }

   //data->getData(true);
   delete data;

   return ret; 
}

// #################################################################
// フィールドデータの出力
template<class T, class TimeT, class TimeAvrT> 
CIO_INLINE
CIO::E_CIO_ERRORCODE
cio_DFI::WriteData(const unsigned step, 
                   TimeT time, 
                   const int sz[3],
                   const int nComp,
                   const int gc, 
                   T* val, 
                   T* minmax,
                   const bool avr_mode, 
                   const unsigned step_avr, 
                   TimeAvrT time_avr)
{

  cio_Array *data = cio_Array::instanceArray
                    ( val
                    , DFI_Finfo.ArrayShape
                    , DFI_Process.RankList[m_RankID].VoxelSize[0]
                    , DFI_Process.RankList[m_RankID].VoxelSize[1]
                    , DFI_Process.RankList[m_RankID].VoxelSize[2]
                    , gc
                    , DFI_Finfo.Component);

  double d_time = (double)time;
  double d_time_avr = (double)time_avr;
  double *d_minmax=NULL;
  if( minmax ) {
    if( DFI_Finfo.Component>1 ) {
      d_minmax = new double[DFI_Finfo.Component*2+2];
      for(int i=0; i<DFI_Finfo.Component*2+2; i++) {
        d_minmax[i] = minmax[i];
      }
    } else { 
      d_minmax = new double[2];
      d_minmax[0] = minmax[0];
      d_minmax[1] = minmax[1];
    }
  }

  CIO::E_CIO_ERRORCODE ret;
  ret = WriteData(step, gc, d_time, data, d_minmax, avr_mode, step_avr, d_time_avr);

  //val = (T*)data->getData(true);
  //data->getData(true);

  if( d_minmax ) delete [] d_minmax;

  delete data;
  return ret;
                                          
}

#endif // _CIO_DFI_INLINE_H_
