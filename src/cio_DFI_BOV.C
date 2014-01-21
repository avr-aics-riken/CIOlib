/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_BOV.C
 * @brief  cio_DFI_BOV Class
 * @author kero    
 */

#include "cio_DFI.h"
#include "cio_DFI_BOV.h"

// #################################################################
// コンストラクタ
cio_DFI_BOV::cio_DFI_BOV()
{

}


// #################################################################
// デストラクタ
cio_DFI_BOV::~cio_DFI_BOV()
{

}

// #################################################################
// ファイルのヘッダーレコード読込み
CIO::E_CIO_ERRORCODE
cio_DFI_BOV::read_HeaderRecord(FILE* fp,
                               bool matchEndian,
                               unsigned step,
                               const int head[3],
                               const int tail[3],
                               int gc,
                               int voxsize[3],
                               double &time)
{

  time=0.0;
  for(int i=0; i<DFI_TimeSlice.SliceList.size(); i++) {
     if( DFI_TimeSlice.SliceList[i].step == step ) {
       time=(double)DFI_TimeSlice.SliceList[i].time;
     }
  }

  for(int i=0; i<3; i++) voxsize[i]=tail[i]-head[i]+1+(2*gc);

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// ファイルのデータレコード読込み
CIO::E_CIO_ERRORCODE
cio_DFI_BOV::read_Datarecord(FILE* fp,
                             bool matchEndian,
                             cio_Array* buf,
                             int head[3],
                             int nz,
                             cio_Array* &src)
{

  //１層ずつ読み込み
  int hzB = head[2];

  for( int k=0; k<nz; k++ ) {
    //headインデクスをずらす
    head[2]=hzB+k;
    buf->setHeadIndex(head);

    //１層読み込
    size_t ndata = buf->getArrayLength();
    if( buf->readBinary(fp,matchEndian) != ndata ) return CIO::E_CIO_ERROR_READ_FIELD_DATA_RECORD;

    // コピー
    buf->copyArray(src);
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// Averaged レコードの読込み
CIO::E_CIO_ERRORCODE
cio_DFI_BOV::read_averaged(FILE* fp,
                           bool matchEndian,
                           unsigned step,
                           unsigned &step_avr,
                           double &time_avr)
{

  step_avr=0;
  time_avr=0.0;

  for(int i=0; i<DFI_TimeSlice.SliceList.size(); i++) {
     if( DFI_TimeSlice.SliceList[i].step == step ) {
       step_avr=(int)DFI_TimeSlice.SliceList[i].AveragedStep;
       time_avr=(double)DFI_TimeSlice.SliceList[i].AveragedTime;
     }
  }

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// ヘッダーレコード出力BOVは何も出力しない
CIO::E_CIO_ERRORCODE
cio_DFI_BOV::write_HeaderRecord(FILE* fp,
                                const unsigned step,
                                const double time,
                                const int n)
{
  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// BOVデータレコード出力
CIO::E_CIO_ERRORCODE
cio_DFI_BOV::write_DataRecord(FILE* fp, 
                              cio_Array* val, 
                              const int gc, 
                              const int n)
{

  CIO::E_CIO_DTYPE Dtype = (CIO::E_CIO_DTYPE)DFI_Finfo.DataType;
  int Real_size = get_cio_Datasize(Dtype);

  int size[3];
  for(int i=0; i<3; i++ ) size[i] = (int)DFI_Process.RankList[n].VoxelSize[i]+(int)(2*gc);

  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.Component > 1 ) dLen *= 3;

  unsigned int dmy = dLen * Real_size;

  if( val->writeBinary(fp) != dLen ) return CIO::E_CIO_ERROR_WRITE_FIELD_HEADER_RECORD;

  return CIO::E_CIO_SUCCESS;
}

// #################################################################
// 平均の出力BOVは何も出力しない
CIO::E_CIO_ERRORCODE
cio_DFI_BOV::write_averaged(FILE* fp,
                            const unsigned step_avr,
                            const double time_avr)
{
  return CIO::E_CIO_SUCCESS;
}
