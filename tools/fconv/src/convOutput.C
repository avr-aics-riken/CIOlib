/*
 * fconv (File Converter)
 *
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   convOutput.C
 * @brief  convOutput Class
 * @author aics
 */

#include "convOutput.h"
#include "convOutput_SPH.h"
#include "convOutput_BOV.h"
#include "convOutput_AVS.h"
#include "convOutput_VTK.h"
#include "convOutput_PLOT3D.h"

// #################################################################
// コンストラクタ
convOutput::convOutput()
{

}

// #################################################################
// デストラクタ
convOutput::~convOutput()
{


}

// #################################################################
// InputParamのポインタをコピー
bool convOutput::importInputParam(InputParam* InputCntl)
{
  if( !InputCntl ) return false;
  m_InputCntl = InputCntl;
  return true;
}

//#################################################################
// convOutput インスタンス
convOutput*
convOutput::OutputInit(const CIO::E_CIO_FORMAT out_format)
{

  convOutput *OutConv = NULL;

  if     ( out_format == CIO::E_CIO_FMT_SPH    ) OutConv = new convOutput_SPH();
  else if( out_format == CIO::E_CIO_FMT_BOV    ) OutConv = new convOutput_BOV();
  else if( out_format == CIO::E_CIO_FMT_AVS    ) OutConv = new convOutput_AVS();
  else if( out_format == CIO::E_CIO_FMT_VTK    ) OutConv = new convOutput_VTK();
  else if( out_format == CIO::E_CIO_FMT_PLOT3D ) OutConv = new convOutput_PLOT3D();

  return OutConv;

}

//#################################################################
//
bool convOutput::WriteFieldData(FILE* fp, cio_Array* src, size_t dLen)
{
  if( src->writeBinary(fp) != dLen ) Exit(0);
  //if( src->writeAscii(fp) != dLen ) Exit(0);
  return true;
}

