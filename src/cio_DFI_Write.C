/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_Write.C
 * @brief  cio_DFI Class
 * @author kero     
 */

#include "cio_DFI.h"

// #################################################################
// Index DFIファイルの出力
CIO::E_CIO_ERRORCODE
cio_DFI::WriteIndexDfiFile(const std::string dfi_name)
{

  if ( dfi_name.empty() ) return CIO::E_CIO_ERROR_WRITE_INDEXFILENAME_EMPTY;
  if ( DFI_Finfo.Prefix.empty() ) return CIO::E_CIO_ERROR_WRITE_PREFIX_EMPTY;

  FILE* fp = NULL;

  // File exist ?
  bool flag = false;
  if ( fp = fopen(dfi_name.c_str(), "r") )
  {
    flag = true;
    fclose(fp);
  }

  if( !(fp = fopen(dfi_name.c_str(), "w")) )
  {
    fprintf(stderr, "Can't open file.(%s)\n", dfi_name.c_str());
    return CIO::E_CIO_ERROR_WRITE_INDEXFILE_OPENERROR;
  }

  //FileInfo {} の出力
  if( DFI_Finfo.Write(fp, 0) != CIO::E_CIO_SUCCESS )
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FILEINFO;
  }

  //FilePath {} の出力
  if( DFI_Fpath.Write(fp, 1) != CIO::E_CIO_SUCCESS )
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FILEPATH;
  }


  //Unit {} の出力
  if( DFI_Unit.Write(fp, 0) != CIO::E_CIO_SUCCESS ) 
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_UNIT;
  }

  //TimeSlice {} の出力
  if ( DFI_TimeSlice.Write(fp, 1) != CIO::E_CIO_SUCCESS ) 
  {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_TIMESLICE;
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// proc DFIファイルの出力コントロール (float 版)
CIO::E_CIO_ERRORCODE
cio_DFI::WriteProcDfiFile(const MPI_Comm comm,
                          bool out_host,
                          float* org)
{

  //orign の再設定
  double d_org[3];
  if( org != NULL ) {
    for(int i=0; i<3; i++) {
      d_org[i]=(double)org[i];
    }
  } else {
    for(int i=0; i<3; i++) {
      d_org[i]=DFI_Domain.GlobalOrigin[i];
    }
  }

  return WriteProcDfiFile(comm, out_host, d_org);

}
// #################################################################
// proc DFIファイルの出力コントロール (double 版)
CIO::E_CIO_ERRORCODE
cio_DFI::WriteProcDfiFile(const MPI_Comm comm,
                          bool out_host,
                          double* org)
{

  //procファイル名の生成
  std::string procFileName = CIO::cioPath_DirName(m_indexDfiName)+"/"+CIO::cioPath_FileName(DFI_Fpath.ProcDFIFile,".dfi");

  if( procFileName.empty() ) return CIO::E_CIO_ERROR_WRITE_PROCFILENAME_EMPTY;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  cio_MPI out_mpi;
  out_mpi.NumberOfRank = nrank;
  out_mpi.NumberOfGroup = 1;

  cio_Domain out_domain;
  cio_Process out_Process;

  //出力するProcess情報の生成
  cio_Create_dfiProcessInfo(comm, out_Process);

  //orign の設定
  if( org!=NULL ) {
    for(int i=0; i<3; i++) {
      out_domain.GlobalOrigin[i] = org[i];
    }
  } else {
    for(int i=0; i<3; i++) {
      out_domain.GlobalOrigin[i] = DFI_Domain.GlobalOrigin[i];
    }
  }

  //Domain の設定
  for(int i=0; i<3; i++) {
    out_domain.GlobalVoxel[i]    = DFI_Domain.GlobalVoxel[i];
    out_domain.GlobalDivision[i] = DFI_Domain.GlobalDivision[i];
    out_domain.GlobalRegion[i]   = DFI_Domain.GlobalRegion[i];
  }

  //ホスト名出力指示ありの時、各ランクのホスト名を集める
  if( out_host ) {
    const int LEN=256;
    char *recbuf = new char[out_Process.RankList.size()*LEN];
    char  sedbuf[LEN];
    //sprintf(sedbuf,"%s",hostname.c_str());
    sprintf(sedbuf,"%s",DFI_Process.RankList[RankID].HostName.c_str());
    MPI_Gather(sedbuf,LEN,MPI_CHAR,recbuf,LEN,MPI_CHAR,0,MPI_COMM_WORLD);

    for( int i=0; i<out_Process.RankList.size(); i++ ) {
     char* hn =&(recbuf[i*LEN]);
     out_Process.RankList[i].HostName=(std::string(hn));
    }

    if( recbuf ) delete [] recbuf;
  }

  //proc.dfの出力
  if( RankID == 0 ) {

    FILE* fp = NULL;
    if( !(fp = fopen(procFileName.c_str(), "w")) )
    {
      fprintf(stderr, "Can't open file.(%s)\n", procFileName.c_str());
      return CIO::E_CIO_ERROR_WRITE_PROCFILE_OPENERROR;
    }

    //Domain {} の出力
    if( out_domain.Write(fp, 0) != CIO::E_CIO_SUCCESS )
    {
      if (fp) fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_DOMAIN;
    }

    //MPI {} の出力
    if( out_mpi.Write(fp, 0) != CIO::E_CIO_SUCCESS )
    {
      fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_MPI;
    }

    //Process {} の出力
    if( out_Process.Write(fp, 0) != CIO::E_CIO_SUCCESS )
    {
      fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_PROCESS;
    }

    fclose(fp);
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// fileld data 出力
CIO::E_CIO_ERRORCODE
cio_DFI::WriteData(const unsigned step,
                   const int gc,
                   double time,
                   cio_Array* val,
                   double* minmax,
                   const bool avr_mode,
                   const unsigned step_avr,
                   double time_avr,
                   bool force)
{

  //インターバルチェック
  if( !m_intervalMngr.isTriggered(step, time, force ) ) return CIO::E_CIO_SUCCESS;

  bool mio=false;
  if( DFI_MPI.NumberOfRank > 1 ) mio=true;
  std::string outFile;
  if( CIO::cioPath_isAbsolute(DFI_Finfo.DirectoryPath) ){
    outFile = Generate_FieldFileName(m_RankID,step,mio);
  } else {
    outFile = m_directoryPath + "/"+ Generate_FieldFileName(m_RankID,step,mio);
  }

  std::string dir = CIO::cioPath_DirName(outFile);

  if( MakeDirectory(dir) != 1 ) return CIO::E_CIO_ERROR_MAKEDIRECTORY;

  cio_Array *outArray = val;
  if( gc != DFI_Finfo.GuideCell ) {
    //出力用バッファのインスタンス
    outArray = cio_Array::instanceArray
               ( DFI_Finfo.DataType
               , DFI_Finfo.ArrayShape
               , DFI_Process.RankList[m_RankID].VoxelSize
               , DFI_Finfo.GuideCell
               , DFI_Finfo.Component); 
    //配列のコピー val -> outArray
    int ret = val->copyArray(outArray);
  }

  //フィールデータの出力
  CIO::E_CIO_ERRORCODE err = CIO::E_CIO_SUCCESS;
  err = WriteFieldData(outFile, step, time, outArray, avr_mode, step_avr, time_avr);

  //出力バッファのメモリ解放
  if( val != outArray ) {
    delete outArray;
  }

  if( err != CIO::E_CIO_SUCCESS ) return err;

  //index dfi ファイルのディレクトリ作成
  cio_DFI::MakeDirectory(m_directoryPath);
  std::string dfiname = CIO::cioPath_FileName(m_indexDfiName,".dfi");
  std::string fname = CIO::cioPath_ConnectPath( m_directoryPath, dfiname );

  //Slice へのセット
  DFI_TimeSlice.AddSlice(step, time, minmax, DFI_Finfo.Component, avr_mode,
                         step_avr, time_avr);

  //index dfi のファイル出力
  if( m_RankID == 0 ) {
    err = WriteIndexDfiFile(fname);
  }

  return err;
}

// #################################################################
// フィールドデータ出力
CIO::E_CIO_ERRORCODE
cio_DFI::WriteFieldData(std::string fname,
                        const unsigned step,
                        double time,
                        cio_Array *val,
                        const bool avr_mode,
                        const unsigned step_avr,
                        const double time_avr)
{

  FILE* fp;
  if( (fp = fopen(fname.c_str(),"wb")) == NULL ) {
    fprintf(stderr,"Can't open file.(%s)\n",fname.c_str());
    return CIO::E_CIO_ERROR_OPEN_FIELDDATA;
  }

  //ヘッダー出力
  if( write_HeaderRecord(fp, step, time, m_RankID) != CIO::E_CIO_SUCCESS ) {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FIELD_HEADER_RECORD;
  }

  //データ出力
  if( write_DataRecord(fp, val, DFI_Finfo.GuideCell, m_RankID) != CIO::E_CIO_SUCCESS) {
    fclose(fp);
    return CIO::E_CIO_ERROR_WRITE_FIELD_DATA_RECORD;
  }

  //average 出力
  if( !avr_mode ) {
    if( write_averaged(fp, step_avr, time_avr) != CIO::E_CIO_SUCCESS ) {
      fclose(fp);
      return CIO::E_CIO_ERROR_WRITE_FIELD_AVERAGED_RECORD;
    }
  }

  fclose(fp);

  return CIO::E_CIO_SUCCESS;

}

