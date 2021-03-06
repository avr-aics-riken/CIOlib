#ifndef _CIO_DFI_VTK_H_
#define _CIO_DFI_VTK_H_

/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013-2014 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_VTK.h
 * @brief  cio_DFI_VTK Class Header
 * @author aics    
 */

#include "cio_DFI.h"

class cio_DFI_VTK : public cio_DFI {

protected:

public:

  /** コンストラクタ */
  cio_DFI_VTK();

  /** 
   * @brief コンストラクタ 
   * @param [in] F_Info  FileInfo
   * @param [in] F_Path  FilePath
   * @param [in] unit    Unit
   * @param [in] domain  Domain
   * @param [in] mpi     MPI
   * @param [in] TSlice  TimeSlice
   * @param [in] process Process
   */
  cio_DFI_VTK(const cio_FileInfo F_Info, 
              const cio_FilePath F_Path, 
              const cio_Unit unit, 
              const cio_Domain domain, 
              const cio_MPI mpi,
              const cio_TimeSlice TSlice, 
              const cio_Process process)
  {
    DFI_Finfo      = F_Info; 
    DFI_Fpath      = F_Path;
    DFI_Unit       = unit;
    DFI_Domain     = domain;
    DFI_MPI        = mpi;
    DFI_TimeSlice  = TSlice;
    DFI_Process    = process;
    m_bgrid_interp_flag = true;
  };
  
  /**　デストラクタ */
  ~cio_DFI_VTK();

public:

protected:

  /**
   * @brief sphファイルのヘッダーレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian エンディアンチェックフラグ true:合致
   * @param[in]  step        ステップ番号
   * @param[in]  head        dfiのHeadIndex
   * @param[in]  tail        dfiのTailIndex
   * @param[in]  gc          dfiのガイドセル数
   * @param[out] voxsize     voxsize
   * @param[out] time        時刻
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  read_HeaderRecord(FILE* fp, 
                    bool matchEndian,
                    unsigned step,
                    const int head[3],
                    const int tail[3],
                    int gc, 
                    int voxsize[3],
                    double &time)
  { return CIO::E_CIO_SUCCESS; };

  /**
   * @brief フィールドデータファイルのデータレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  buf         読込み用バッファ
   * @param[in]  head        読込みバッファHeadIndex
   * @param[in]  nz          z方向のボクセルサイズ（実セル＋ガイドセル＊２）
   * @param[out] src         読み込んだデータを格納した配列のポインタ
   * @return error code
   */ 
  CIO::E_CIO_ERRORCODE
  read_Datarecord(FILE* fp,
                  bool matchEndian,
                  cio_Array* buf,
                  int head[3],
                  int nz,
                  cio_Array* &src)
  { return CIO::E_CIO_SUCCESS; };

  /**
   * @brief sphファイルのAverageデータレコードの読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        読込みstep番号
   * @param[out] avr_step    平均ステップ   
   * @param[out] avr_time    平均タイム
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  read_averaged(FILE* fp,
                bool matchEndian,
                unsigned step, 
                unsigned &avr_step,
                double &avr_time) 
  { return CIO::E_CIO_SUCCESS; };

  /**
   * @brief VTKヘッダファイルの出力
   * @param[in] fp     ファイルポインタ
   * @param[in] step   ステップ番号
   * @param[in] time   時刻
   * @param[in] RankID ランク番号
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  write_HeaderRecord(FILE* fp, 
                     const unsigned step, 
                     const double time, 
                     const int RankID); 

  /**
   * @brief VTKデータレコードの出力
   * @param[in]  fp ファイルポインタ
   * @param[in]  val データポインタ
   * @param[in]  gc ガイドセル
   * @param[in]  RankID ランク番号
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  write_DataRecord(FILE* fp, 
                   cio_Array* val, 
                   const int gc, 
                   const int RankID); 

  /**
   * @brief Averageレコードの出力
   * @param[in] fp       ファイルポインタ
   * @param[in] step_avr 平均ステップ番号
   * @param[in] time_avr 平均時刻
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  write_averaged(FILE* fp,
                 const unsigned step_avr,
                 const double time_avr) 
  { return CIO::E_CIO_SUCCESS; };
};

#endif // _cio_DFI_VTK_H_
