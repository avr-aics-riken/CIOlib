#ifndef _CIO_DFI_BOV_H_
#define _CIO_DFI_BOV_H_

/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_BOV.h
 * @brief  cio_DFI_BOV Class Header
 * @author kero    
 */

#include "cio_DFI.h"

using namespace std;


class cio_DFI_BOV : public cio_DFI {

public:
  /** コンストラクタ */
  cio_DFI_BOV();

  cio_DFI_BOV(cio_FileInfo F_Info, cio_FilePath F_Path, cio_Unit unit, cio_Domain domain, cio_MPI mpi,
              vector<cio_Slice> TSlice, vector<cio_Rank> RInfo)
  {
    DFI_Finfo  = F_Info;
    DFI_Fpath  = F_Path;
    DFI_Unit   = unit;
    DFI_Domain = domain;
    DFI_MPI    = mpi;
    TimeSlice  = TSlice;
    RankInfo   = RInfo;
  };

  /**　デストラクタ */
  ~cio_DFI_BOV();

  /**
   * @brief read sph data
   * @param[in]  step         読込むstep番号
   * @param[in]  gc           仮想セル数
   * @param[in]  Gvoxel[3]    グローバルボクセルサイズ
   * @param[in]  Gdivision[3] 領域分割数
   * @param[in]  head[3]      計算領域の開始位置
   * @param[in]  tail[3]      計算領域の終了位置
   * @param[out] val          フィールドデータポインタ
   */
  void ReadData(int step, int gc,
                int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                void *val, double &time,
                const bool mode, unsigned &step_avr, double &time_avr);

  void *ReadData(int step, int gc,
                int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                double &time,
                const bool mode, unsigned &step_avr, double &time_avr);

  /**
   * @brief bov S3D ファイルの読込み 
   * @param[in]  fname bovファイル名
   * @param[in]  step  読込むstep番号
   * @param[in]  gc    仮想セル数
   * @param[out] val   フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */
  bool read_S3D(const char* fname,
                int         step,
                int         gc,
                void*       val,
                int         voxsize[3],
                double&     time
                );

  /**
   * @brief 密データの読込みコピー（同一粒子）
   * @param[in]  fname  sphファイル名
   * @param[in]  step   読込むstep番号
   * @param[in]  gc     仮想セル数
   * @param[in]  head   自領域の起点
   * @param[in]  tail   自領域の終点
   * @param[in]  R_head 読み込む領域の起点
   * @param[in]  R_tail 読み込む領域の終点
   * @param[in]  sta[3] start
   * @param[in]  end[3] end
   * @param[out] val    フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */
  bool read_MxN(const char* fname,
                int         step,
                int         gc,
                int         head[3],
                int         tail[3],
                int         R_head[3],
                int         R_tail[3],
                int         sta[3],
                int         end[3],
                void*       val,
                double&     time,
                int         n,
                int         voxsize[3]);

  template<class T>
  void setval(T* val, int p1, T* val2, int p2)
  {
     val[p1] = val2[p2];
  };

  /**
   * @brief 粗い粒子、１対１の読込みコピー&補間
   * @param[in]  fname     bovファイル名 
   * @param[in]  step      読込むstep番号
   * @param[in]  gc        仮想セル数（自）
   * @param[in]  head      開始インデックス（自）
   * @param[in]  tail      終了インデックス（自）
   * @param[in]  R_head    開始インデックス（dfi）
   * @param[in]  R_tail    終了インデックス（dif）
   * @param[in]  Voxelsize ボクセルサイズ（dif）
   * @param[in]  dfi_gc    仮想セル数（dfi）
   * @param[in]  ncomp     コンポーネント数
   * @param[in]  sta[3]    start
   * @param[in]  end[3]    end
   * @param[out] val       フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */
  bool read_Coarse(const char* fname,
                int         step,
                int         gc,
                int         head[3],
                int         tail[3],
                int         G_Voxelsize[3],
                int         R_head[3],
                int         R_tail[3],
                int         Voxelsize[3],
                int         dfi_gc,
                int         ncomp,
                int         sta[3],
                int         end[3],
                void*       val,
                double&     time);

  /**
   * @brief 粗い粒子、１対多の読込みコピー&補間
   * @param[in]  fname     bovファイル名
   * @param[in]  step      読込むstep番号
   * @param[in]  gc        仮想セル数（自）
   * @param[in]  head      開始インデックス（自）
   * @param[in]  tail      終了インデックス（自）
   * @param[in]  R_head    開始インデックス（dfi）
   * @param[in]  R_tail    終了インデックス（dif）
   * @param[in]  Voxelsize ボクセルサイズ（dif）
   * @param[in]  dfi_gc    仮想セル数（dfi）
   * @param[in]  ncomp     コンポーネント数
   * @param[in]  sta[3]    start
   * @param[in]  end[3]    end
   * @param[out] val       フィールドデータポインタ
   * @return true:出力成功 false:出力失敗
   */
  bool read_Coarse_MxN(const char* fname,
                int         step,
                int         gc,
                int         head[3],
                int         tail[3],
                int         G_Voxelsize[3],
                int         R_head[3],
                int         R_tail[3],
                int         Voxelsize[3],
                int         dfi_gc,
                int         ncomp,
                int         sta[3],
                int         end[3],
                void*       val,
                double&     time,
                int         n);

  /**
   * @brief write sph data
   * @param[in]  step 出力step番号
   * @param[in]  gc      仮想セル数　　　
   * @param[in]  time    出力時刻　　　　
   * @param[in]  val     フィールドデータポインタ
   * @param[in]  minmax
   * @param[in]  interval 出力間隔
   * @param[in]  force    強制出力指示
   */
  void WriteData(int step, int gc, void* time,
                 void* val, void* minmax, int interval, 
                 const bool mode, const unsigned step_avr, const double time_avr,
                 bool force);

  /**
   * @brief BOVデータ出力 
   * @param[in]  fp ファイルポインタ
   * @param[in]  val データポインタ
   * @param[in]  RankID ランク番号
   * @return true:出力成功 false:出力失敗
   */
  bool write_data(FILE* fp, void* val, int gc, int RankID);
  
};

#endif // _cio_DFI_BOV_H_
