#ifndef _CIO_DFI_H_
#define _CIO_DFI_H_

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
 * @author aics
 */

#include "cio_Define.h"
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <typeinfo>
#include <set>
#include <map>
#include <string>

#include "cio_Version.h"


#include "cio_PathUtil.h"
#include "cio_TextParser.h"
#include "cio_ActiveSubDomain.h"
#include "cio_endianUtil.h"
#include "cio_TypeArray.h"

#include "cio_FileInfo.h"
#include "cio_FilePath.h"
#include "cio_Unit.h"
#include "cio_TimeSlice.h"
#include "cio_Domain.h"
#include "cio_MPI.h"
#include "cio_Process.h"

/** CIO main class */
class cio_DFI {

public:

protected :
  MPI_Comm          m_comm;          ///< MPI コミュニケータ
  std::string       m_directoryPath; ///< index dfi ファイルのディレクトリパス
  std::string       m_indexDfiName;  ///< index dfi ファイル名
  CIO::E_CIO_READTYPE m_read_type;   ///< 読込みタイプ

  int               m_RankID;        ///< ランク番号

  cio_FileInfo      DFI_Finfo;       ///< FileInfo class
  cio_FilePath      DFI_Fpath;       ///< FilePath class
  cio_Unit          DFI_Unit;        ///< Unit class
  cio_Domain        DFI_Domain;      ///< Domain class
  cio_MPI           DFI_MPI;         ///< MPI class
  cio_TimeSlice     DFI_TimeSlice;   ///< TimeSlice class
  cio_Process       DFI_Process;     ///< Process class

  vector<int>m_readRankList;         ///< 読込みランクリスト

  bool m_bgrid_interp_flag;               ///< 節点への補間フラグ
  CIO::E_CIO_OUTPUT_TYPE  m_output_type;  ///< 出力形式(ascii,binary,FortarnBinary)
  CIO::E_CIO_OUTPUT_FNAME m_output_fname; ///< 出力ファイル命名規約(step_rank,rank_step)   

public:
  /** コンストラクタ */
  cio_DFI();
  
  /**　デストラクタ */
  ~cio_DFI();

  /**
   * @brief read インスタンス
   * @param [in]  comm    MPIコミュニケータ
   * @param [in]  dfifile DFIファイル名
   * @param [in]  G_Voxel 計算空間全体のボクセルサイズ
   * @param [in]  G_Div   計算空間の領域分割数
   * @param [out] ret     終了コード
   * @return インスタンスされたクラスのポインタ
   */
  static cio_DFI*
  ReadInit(const MPI_Comm comm, 
           const std::string dfifile,
           const int G_Voxel[3],
           const int G_Div[3],
           CIO::E_CIO_ERRORCODE &ret); 

  /**
   * @brief cioFileInfoクラスのポインタを取得
   * @return cio_FileInfoクラスポインタ
   */
  const cio_FileInfo* GetcioFileInfo();

  /**
   * @brief cio_FilePathクラスのポインタを取得
   * @return cio_FilePathクラスポインタ
   */
  const cio_FilePath* GetcioFilePath();  

  /**
   * @brief cio_FilePathクラスのセット
   */
  void SetcioFilePath(cio_FilePath FPath); 

  /**
   * @brief cio_Unitクラスのポインタを取得
   * @return cio_Unitクラスポインタ
   */
  const cio_Unit* GetcioUnit(); 

  /**
   * @brief cio_Unitクラスのセット
   */
  void SetcioUnit(cio_Unit unit); 


  /**
   * @brief cio_Domainクラスのポインタ取得
   * @return cio_Domainクラスポインタ
   */
  const cio_Domain* GetcioDomain(); 

  /**
   * @brief cio_Domainクラスのセット
   */
  void SetcioDomain(cio_Domain domain); 

  /**
   * @brief cio_MPIクラスのポインタ取得
   * @return cio_MPIクラスポインタ 
   */
  const cio_MPI* GetcioMPI();

  /**
   * @brief cio_MPIクラスセット
   */
  void SetcioMPI(cio_MPI mpi);

  /**
   * @brief cio_TimeSliceクラスのポインタ取得
   * @return cio_TimeSliceクラスポインタ
   */
  const cio_TimeSlice* GetcioTimeSlice(); 

  /**
   * @brief cio_TimeSlice クラスセット
   */
  void SetcioTimeSlice(cio_TimeSlice TSlice);
 
  /**
   * @brief cio_Processクラスのポインタ取得
   * @return cio_Processクラスポインタ
   */
  const cio_Process* GetcioProcess(); 

  /**
   * @brief cio_Processクラスセット
   */
  void SetcioProcess(cio_Process Process);

  /**
   * @brief 出力DFIファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @return DFIファイル名
   */ 
  static std::string
  Generate_DFI_Name(const std::string prefix);

  /**
   * @brief フィールドデータ（SPH,BOV)ファイル名の作成(ディレクトリパスが付加されている）
   * @param [in] RankID ランク番号
   * @param [in] step   読込みステップ番号
   * @param [in] mio    並列判定フラグ（逐次or並列の判定用）
   * @return 生成されたファイル名　　　　　　　
   */
  std::string Generate_FieldFileName(int RankID,
                                int step, 
                                const bool mio);

  /**
   * @brief ファイル名生成
   * @param [in] prefix ベースファイル名
   * @param [in] RankID ランク番号  
   * @param [in] step   出力ステップ番号（負のとき、ステップ番号が付加されない）
   * @param [in] ext    拡張子
   * @param [in] output_fname step_rank,rank_step指示
   * @param [in] mio    並列判定フラグ
   * @param [in] TimeSliceDirFlag Time Slice 毎の出力指示
   * @return 生成されたファイル名
   */
  static
  std::string Generate_FileName(std::string prefix,
                                int RankID,
                                int step,
                                std::string ext,
                                CIO::E_CIO_OUTPUT_FNAME output_fname,
                                bool mio,
                                CIO::E_CIO_ONOFF TimeSliceDirFlag);

  /**
   * @brief write インスタンス float型
   * @param [in] comm        MPIコミュニケータ
   * @param [in] DfiName     DFIファイル名
   * @param [in] Path        フィールドデータのディレクトリ
   * @param [in] prefix      ベースファイル名
   * @param [in] format      ファイルフォーマット
   * @param [in] GCell       出力仮想セル数　　　
   * @param [in] DataType    データタイプ　　　　
   * @param [in] ArrayShape  配列形状　　　　　　
   * @param [in] nComp       成分数　　　　　　　
   * @param [in] proc_fname  proc.dfiファイル名
   * @param [in] G_size      グローバルボクセルサイズ　
   * @param [in] pitch       ピッチ　　　　　　　　　　
   * @param [in] G_origin    原点座標値　　　　　　　　
   * @param [in] division    領域分割数　　　　　　　　
   * @param [in] head        計算領域の開始位置　　　　
   * @param [in] tail        計算領域の終了位置　　　　
   * @param [in] hostname    ホスト名
   * @param [in] TSliceOnOff TimeSliceフラグ
   * @return インスタンスされたクラスのポインタ
   */
  static cio_DFI*
  WriteInit(const MPI_Comm comm,
            const std::string DfiName,
            const std::string Path,
            const std::string prefix,
            const CIO::E_CIO_FORMAT format,
            const int GCell,
            const CIO::E_CIO_DTYPE DataType,
            const CIO::E_CIO_ARRAYSHAPE ArrayShape,
            const int nComp,
            const std::string proc_fname,
            const int G_size[3],
            const float pitch[3],
            const float G_origin[3],
            const int division[3],
            const int head[3],
            const int tail[3],
            const std::string hostname,
            const CIO::E_CIO_ONOFF TSliceOnOff);

  /**
   * @brief write インスタンス double型
   * @param [in] comm        MPIコミュニケータ
   * @param [in] DfiName     DFIファイル名
   * @param [in] Path        フィールドデータのディレクトリ
   * @param [in] prefix      ベースファイル名
   * @param [in] format      ファイルフォーマット
   * @param [in] GCell       出力仮想セル数　　　
   * @param [in] DataType    データタイプ　　　　
   * @param [in] ArrayShape  配列形状　　　　　　
   * @param [in] nComp       成分数　　　　　　　
   * @param [in] proc_fname  proc.dfiファイル名
   * @param [in] G_size      グローバルボクセルサイズ　
   * @param [in] pitch       ピッチ　　　　　　　　　　
   * @param [in] G_origin    原点座標値　　　　　　　　
   * @param [in] division    領域分割数　　　　　　　　
   * @param [in] head        計算領域の開始位置　　　　
   * @param [in] tail        計算領域の終了位置　　　　
   * @param [in] hostname    ホスト名　　　　　　　　　
   * @param [in] TSliceOnOff TimeSliceフラグ
   * @return インスタンスされたクラスのポインタ
   */
  static cio_DFI* 
  WriteInit(const MPI_Comm comm,
            const std::string DfiName,
            const std::string Path,
            const std::string prefix,
            const CIO::E_CIO_FORMAT format,
            const int GCell,
            const CIO::E_CIO_DTYPE DataType,
            const CIO::E_CIO_ARRAYSHAPE ArrayShape,
            const int nComp,
            const std::string proc_fname,
            const int G_size[3],
            const double pitch[3],
            const double G_origin[3],
            const int division[3],
            const int head[3],
            const int tail[3],
            const std::string hostname,
            const CIO::E_CIO_ONOFF TSliceOnOff);

  /**
   * @brief RankIDをセットする
   * @param[in] rankID RankID
   */
  void set_RankID(const int rankID)
  { m_RankID = rankID; };
  /**
   * @brief 出力形式(ascii,binary,FortranBinary)をセット
   * @param [in] output_type 出力形式
   */
  void set_output_type(CIO::E_CIO_OUTPUT_TYPE output_type)
  {  m_output_type = output_type; };

  /**
   * @brief 出力ファイル命名規約(step_rank,rank_step)をセット
   * @param [in] output_fname 出力ファイル命名規約
   */
  void set_output_fname(CIO::E_CIO_OUTPUT_FNAME output_fname)
  { m_output_fname = output_fname; };

  /**
   * @brief DFIファイル名の取り出し
   * @return dfiファイル名
   */ 
  std::string get_dfi_fname()
  { return m_indexDfiName; };


  /**
   * @brief read field data record (template function)
   * @details 読み込んだデータのポインタを戻り値として返す
   * @param [out] ret       終了コード 1:正常、1以外：エラー  
   * @param [in] step       入力ステップ番号
   * @param [in] gc         仮想セル数　　　
   * @param [in] Gvoxel     グローバルボクセルサイズ　
   * @param [in] Gdivision  領域分割数　　　　　　　　
   * @param [in] head       計算領域の開始位置　　　　
   * @param [in] tail       計算領域の終了位置　　　　
   * @param [out] time      読み込んだ時間
   * @param [in]  mode      平均ステップ＆時間読込みフラグ　false : 読込み
   *                                                           true  : 読み込まない
   * @param [out] step_avr  平均ステップ
   * @param [out] time_avr  平均時間
   * @return 読みんだフィールドデータのポンタ
   */
//  template<class T, class TimeT, class TimeAvrT> T*
  template<class TimeT, class TimeAvrT> void*
  ReadData(CIO::E_CIO_ERRORCODE &ret,
           const unsigned step, 
           const int gc, 
           const int Gvoxel[3], 
           const int Gdivision[3], 
           const int head[3], 
           const int tail[3],
           TimeT &time,
           const bool mode, 
           unsigned &step_avr, 
           TimeAvrT &time_avr);
  /**
   * @brief read field data record (template function)
   * @details 引数で渡された配列ポインタにデータを読込む
   * @param [out] val        読み込んだデータポインタ　
   * @param [in]  step       入力ステップ番号
   * @param [in]  gc         仮想セル数　　　
   * @param [in]  Gvoxel     グローバルボクセルサイズ　
   * @param [in]  Gdivision  領域分割数　　　　　　　　
   * @param [in]  head       計算領域の開始位置　　　　
   * @param [in]  tail       計算領域の終了位置　　　　
   * @param [out] time       読み込んだ時間
   * @param [in]  mode       平均ステップ＆時間読込みフラグ　false : 読込み
   *                                                           true  : 読み込まない
   * @param [out] step_avr      平均ステップ
   * @param [out] time_avr      平均時間
   * @return 終了コード 1:正常 1以外:エラー
   */
  template<class T, class TimeT, class TimeAvrT>
  CIO::E_CIO_ERRORCODE 
  ReadData(T *val,
           const unsigned step,
           const int gc,
           const int Gvoxel[3],
           const int Gdivision[3],
           const int head[3],
           const int tail[3],
           TimeT &time,
           const bool mode,
           unsigned &step_avr,
           TimeAvrT &time_avr);

  /**
   * @brief read field data record 
   * @details template ReadData関数で型に応じた配列を確保した後、呼び出される
   * @param [out] val        読み込み先の配列をポインタで渡す　
   * @param [in]  step       読み込むステップ番号
   * @param [in]  gc         仮想セル数　　　
   * @param [in]  Gvoxel     グローバルボクセルサイズ　
   * @param [in]  Gdivision  領域分割数　　　　　　　　
   * @param [in]  head       計算領域の開始位置　　　　
   * @param [in]  tail       計算領域の終了位置　　　　
   * @param [out] time       読み込んだ時間
   * @param [in]  mode       平均ステップ＆時間読込みフラグ　false : 読込み
   *                                                           true  : 読み込まない
   * @param [out] step_avr      平均ステップ
   * @param [out] time_avr      平均時間
   * @return 終了コード 1:正常 1以外:エラー
   */
  CIO::E_CIO_ERRORCODE 
  ReadData(cio_Array *val,
           const unsigned step,
           const int gc,
           const int Gvoxel[3],
           const int Gdivision[3],
           const int head[3],
           const int tail[3],
           double &time,
           const bool mode,
           unsigned &step_avr,
           double &time_avr);

  /**
   * @brief write field data record (template function)
   * @details スカラーのとき、minmax[0]=min 
   *                          minmax[1]=max
   *          ベクトルのとき、minmax[0]   =成分1のminX
   *                          minmax[1]   =成分1のmaxX
   *                               ...
   *                          minmax[2n-2]=成分nのminX
   *                          minmax[2n-1]=成分nのmaxX
   *                          minmax[2n  ]=合成値のmin
   *                          minmax[2n+1]=合成値のmax
   * @param [in] step     出力ステップ番号
   * @param [in] time     出力時刻　　　　
   * @param [in] sz       valの実ボクセルサイズ
   * @param [in] nComp    valの成分数（1or3)
   * @param [in] gc       valの仮想セル数　　　
   * @param [in] val      出力データポインタ
   * @param [in] minmax   フィールデータのMinMax
   * @param [in] avr_mode 平均ステップ＆時間出力　false : 出力 true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   */ 
  template<class T, class TimeT, class TimeAvrT>
  CIO::E_CIO_ERRORCODE
  WriteData(const unsigned step, 
            TimeT time,
            const int sz[3], 
            const int nComp,
            const int gc, 
            T* val, 
            T* minmax=NULL, 
            bool avr_mode=true, 
            unsigned step_avr=0, 
            TimeAvrT time_avr=0.0);

  /**
   * @brief write field data record
   * @details template WriteData関数で方に応じた配列を確保した後、呼び出される 
   * @param [in] step     出力ステップ番号
   * @param [in] gc       仮想セル数　　　
   * @param [in] time     出力時刻　　　　
   * @param [in] val      出力データポインタ
   * @param [in] minmax   フィールデータのMinMax
   * @param [in] avr_mode 平均ステップ＆時間出力　false : 出力
   *                                              true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   */ 
  CIO::E_CIO_ERRORCODE
  WriteData(const unsigned step, 
            const int gc, 
            double time, 
            cio_Array* val, 
            double* minmax, 
            const bool avr_mode, 
            const unsigned step_avr, 
            double time_avr);

  /**
   * @brief proc DFIファイル出力コントロール (float)
   * @param [in] comm      MPIコミュニケータ
   * @param [in] out_host  ホスト名出力フラグ　　　　
   * @return true:出力成功 false:出力失敗
   */
/*
  CIO::E_CIO_ERRORCODE
  WriteProcDfiFile(const MPI_Comm comm, 
                   bool out_host=false);
                   float* org=NULL);
*/
  /**
   * @brief proc DFIファイル出力コントロール
   * @param [in] comm          MPIコミュニケータ
   * @param [in] out_host      ホスト名出力フラグ　　　　
   * @return true:出力成功 false:出力失敗
   */
  CIO::E_CIO_ERRORCODE
  WriteProcDfiFile(const MPI_Comm comm, 
                   bool out_host=false);
                   //double* org=NULL);

  /**
   * @brief 配列形状を文字列で返す
   * @return 配列形状（文字列)
   */
  std::string 
  GetArrayShapeString();

  /**
   * @brief 配列形状を返す
   * @return 配列形状（e_num番号)
   */
  CIO::E_CIO_ARRAYSHAPE 
  GetArrayShape();

  /**
   * @brief get DataType （データタイプの取り出し関数）
   * @return データタイプ（文字列)
   */
  std::string 
  GetDataTypeString();

  /**
   * @brief get DataType （データタイプの取り出し関数）
   * @return データタイプ(e_num番号)
   */
  CIO::E_CIO_DTYPE 
  GetDataType();

  /** 
   * @brief get FileFormat （FileFormatの取り出し関数） 
   * @return FileFormat(文字列)
   */
  std::string
  GetFileFormatString();

  /** @brief get FileFormat (FileFormatの取り出し関数）
   *  @return FileFormat(e_num番号)
   */
  CIO::E_CIO_FORMAT
  GetFileFormat();

  /**
   * @brief get Number of Component （成分数の取り出し関数）
   * @return 成分数
   */
  int 
  GetNumComponent();

  /*
   * @brief get Number of GuideCell (仮想セル数の取り出し関数)
   * @return 仮想セル数
   */
  int
  GetNumGuideCell();

  /**
   * @brief データタイプを文字列からe_num番号に変換 
   * @param [in] datatype dfiから取得したデータタイプ
   * @return データタイプ(E_CIO_DTYPE)
   */
  static CIO::E_CIO_DTYPE 
  ConvDatatypeS2E(const std::string datatype); 

  /**
   * @brief データタイプをe_num番号から文字列に変換 
   * @param [in] Dtype データタイプ
   * @return データタイプ(string)
   */
  static std::string 
  ConvDatatypeE2S(const CIO::E_CIO_DTYPE Dtype); 

  /**
   * @brief DFI DomainのGlobalVoxelの取り出し
   * @return GlobalVoxelのポインタ
   */
  int* 
  GetDFIGlobalVoxel(); 

  /**
   * @brief DFI DomainのGlobalDivisionの取り出し
   * @return GlobalDivisionのポインタ
   */
  int* 
  GetDFIGlobalDivision(); 

  /**
   * @brief Uuitをセットする
   * @param [in] Name       追加する単位系("Length","Velocity",,,,)
   * @param [in] Unit       単位ラベル("M","CM","MM","M/S",,,)
   * @param [in] reference  規格化したスケール値
   * @param [in] difference 差の値
   * @param [in] BsetDiff   differenceの有無
   */
  void 
  AddUnit(const std::string Name,
          const std::string Unit,
          const double reference,
          const double difference= 0.0,
          const bool BsetDiff=false);

  /**
   * @brief UuitElemを取得する
   * @param[in]  Name 取得する単位系
   * @param[out] unit 取得したcio_UnitElem
   * @return error code
   */
  CIO::E_CIO_ERRORCODE GetUnitElem(const std::string Name,
                                   cio_UnitElem &unit); 

  /**
   * @brief UnitElemのメンバ変数毎に取得する
   * @param[in]  Name 取得する単位系
   * @param[out] unit 単位文字列
   * @param[out] ref  reference
   * @param[out] diff difference
   * @param[out] bSetDiff differenceの有無（true:あり false:なし）
   * @return error code
   */
  CIO::E_CIO_ERRORCODE GetUnit(const std::string Name,
                               std::string &unit,
                               double &ref,
                               double &diff,
                               bool &bSetDiff); 
 
  /**
   * @brief TimeSlice OnOff フラグをセットする
   * @param [in] ONOFF
   */
  void 
  SetTimeSliceFlag(const CIO::E_CIO_ONOFF ONOFF); 

  /**
   * @brief FileInfoの成分名を登録する
   * @param [in] pcomp    成分位置 0:u, 1:v, 2:w
   * @param [in] compName 成分名 "u","v","w",,,
   */
  void setComponentVariable(int pcomp, std::string compName); 

  /**
   * @brief FileInfoの成分名を取得する
   * @param [in] pcomp 成分位置 0:u, 1:v, 2:w
   * @return 成分名
   */
  std::string getComponentVariable(int pcomp);

  /**
   * @brief DFIに出力されているminmaxの合成値を取得
   * @param [in]  step 取得するステップ
   * @param [out] vec_min 取得したminmaxの合成値
   * @param [out] vec_max 取得したminmaxの合成値
   * @return error code 取得出来たときは E_CIO_SUCCESS
   */
  CIO::E_CIO_ERRORCODE getVectorMinMax(const unsigned step,
                                       double &vec_min,
                                       double &vec_max);

  /**
   *brief DFIに出力されているminmaxを取得
   * @param [in]  step 取得するステップ
   * @param [in]  compNo 成分No(0～n)
   * @param [out] min_value 取得したmin
   * @param [out] max_value 取得したmax
   * @return error code 取得出来たときは E_CIO_SUCCESS
   */
  CIO::E_CIO_ERRORCODE getMinMax(const unsigned step,
                                 const int compNo,
                                 double &min_value,
                                 double &max_value);

  /**
   * @brief 読込みランクリストの作成
   * @details RankListがあるかないか判定しないときは新規にRankListを生成し
   *          それをもとにランクマップの生成、読込みランクリストreadRankList
   *          を生成する
   * @param [in]  dfi_domain DFIのdomain情報
   * @param [in]  head       ソルバーのHeadIndex
   * @param [in]  tail       ソルバーのTailIndex
   * @param [in]  readflag   読込み方法
   * @param [out] readRankList 読込みランクリスト
   * @return error code
   */
  CIO::E_CIO_ERRORCODE
  CheckReadRank(cio_Domain dfi_domain,
                const int head[3],
                const int tail[3],
                CIO::E_CIO_READTYPE readflag,
                vector<int> &readRankList);

  /**
   * @brief 出力インターバルステップの登録
   * @details 登録しない（本メソッドがコールされない）場合はCIOでのインターバル
   *          制御は行わない
   * @param [in] interval_step インターバルステップ
   * @param [in] base_step     基準となるステップ（デフォルト0ステップ）
   * @param [in] start_step    セッション開始ステップ（デフォルト0ステップ）
   * @param [in] last_step     セッション最終ステップ（デフォルト、-1：最終ステップで出力しない）
   */
  void setIntervalStep(int interval_step,
                       int base_step =0, 
                       int start_step=0,
                       int last_step =-1); 

  /**
   * @brief インターバルタイムの登録
   * @param [in] interval_time 出力インターバルタイム
   * @param [in] dt            計算の時間間隔
   * @param [in] base_time     基準となるタイム（デフォルト0.0タイム）
   * @param [in] start_time    セッション開始タイム（デフォルト0.0タイム）
   * @param [in] last_time     せっしょん最終タイム（デフォルト、-1.0：最終タイムで出力しない）
   */
  void setIntervalTime(double interval_time,
                       double dt,
                       double base_time =0.0, 
                       double start_time=0.0,
                       double last_time =-1.0); 

  /**
   * @brief インターバルの計算に使われる全ての時間をスケールで無次元化する
   * @details (base_time, interval_time, start_time, last_time)
   * @param [in] scale スケール
   * return modeがStepのときはfalseを返す、無次元化しない
   */
  bool normalizeTime(const double scale);

  /**
   * @brief インターバルのbase_timeをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeBaseTime(const double scale);
 
  /**
   * @brief インターバルのintervalをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeIntervalTime(const double scale);
 
  /**
   * @brief インターバルのstart_timeをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeStartTime(const double scale);

  /**
   * @brief インターバルのlast_timeをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeLastTime(const double scale);

  /**
   * @brief インターバルのDetlaTをスケールで無次元化する
   * @param [in] scale スケール
   */
  void normalizeDelteT(const double scale);



  /**
   * @brief read field data record(sph or bov)
   * @param [in]  fname    FieldData ファイル名
   * @param [in]  step     読込みステップ番号
   * @param [out] time     読み込んだ時間
   * @param [in]  sta      読込みスタート位置
   * @param [in]  end      読込みエンド位置
   * @param [in]  DFI_head dfiのHeadIndex
   * @param [in]  DFI_tail dfiのTailIndex
   * @param [in]  avr_mode 平均ステップ＆時間読込みフラグ　false : 読込み
   * @details　                                             true  : 読み込まない
   * @param [out] avr_step 平均ステップ
   * @param [out] avr_time 平均時間
   * @param [out] ret      終了コード
   * @return 読み込んだ配列のポインタ
   */
  virtual
  cio_Array* 
  ReadFieldData(std::string fname,
                const unsigned step,
                double &time,
                const int sta[3],
                const int end[3],
                const int DFI_head[3],
                const int DFI_tail[3],
                bool avr_mode,
                unsigned &avr_step,
                double &avr_time,
                CIO::E_CIO_ERRORCODE &ret );

  /**
   * @brief フィールドデータファイルのヘッダーレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        ステップ番号
   * @param[in]  head        dfiのHeadIndex
   * @param[in]  tail        dfiのTailIndex
   * @param[in]  gc          dfiのガイドセル数
   * @param[out] voxsize     voxsize
   * @param[out] time        時刻
   * @return true:出力成功 false:出力失敗
   */
  virtual CIO::E_CIO_ERRORCODE 
  read_HeaderRecord(FILE* fp,
                    bool matchEndian,
                    unsigned step,
                    const int head[3],
                    const int tail[3],
                    int gc,
                    int voxsize[3],
                    double &time)=0;

  /**
   * @brief フィールドデータファイルのデータレコード読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  buf         読込み用バッファ
   * @param[in]  head        読込みバッファHeadIndex
   * @param[in]  nz          z方向のボクセルサイズ（実セル＋ガイドセル＊２）
   * @param[out] src         読み込んだデータを格納した配列のポインタ
   */
  virtual CIO::E_CIO_ERRORCODE 
  read_Datarecord(FILE* fp,
                  bool matchEndian,
                  cio_Array* buf,
                  int head[3],
                  int nz,
                  cio_Array* &src)=0;

  /**
   * @brief sphファイルのAverageデータレコードの読込み
   * @param[in]  fp          ファイルポインタ
   * @param[in]  matchEndian true:Endian一致
   * @param[in]  step        読込みstep番号
   * @param[out] avr_step    平均ステップ
   * @param[out] avr_time    平均タイム
   */
  virtual CIO::E_CIO_ERRORCODE
  read_averaged(FILE* fp,
                bool matchEndian,
                unsigned step,
                unsigned &avr_step,
                double &avr_time)=0;   

protected :
  /**
   * @brief write field data record (double)
   * @param [in] fname    出力フィールドファイル名
   * @param [in] step     出力ステップ番号
   * @param [in] time     出力時刻　　　　
   * @param [in] val      出力データポインタ
   * @param [in] mode     平均ステップ＆時間出力　false : 出力 true  : 出力しない
   * @param [in] step_avr 平均ステップ
   * @param [in] time_avr 平均時間
   * @return error code
   */
  virtual
  CIO::E_CIO_ERRORCODE 
  WriteFieldData(std::string fname,
                 const unsigned step, 
                 double time, 
                 cio_Array* val, 
                 const bool mode, 
                 const unsigned step_avr, 
                 const double time_avr);

  /**
   * @brief SPHヘッダファイルの出力
   * @param[in] fp     ファイルポインタ
   * @param[in] step   ステップ番号
   * @param[in] time   時刻
   * @param[in] RankID ランク番号
   * @return true:出力成功 false:出力失敗
   */
  virtual CIO::E_CIO_ERRORCODE
  write_HeaderRecord(FILE* fp,
                     const unsigned step,
                     const double time,
                     const int RankID)=0;

  /**
   * @brief SPHデータレコードの出力
   * @param[in]  fp ファイルポインタ
   * @param[in]  val データポインタ
   * @param[in]  gc ガイドセル
   * @param[in]  RankID ランク番号
   * @return true:出力成功 false:出力失敗
   */
  virtual CIO::E_CIO_ERRORCODE
  write_DataRecord(FILE* fp,
                   cio_Array* val,
                   const int gc,
                   const int RankID)=0;

  /**
   * @brief Averageレコードの出力
   * @param[in] fp       ファイルポインタ
   * @param[in] step_avr 平均ステップ番号
   * @param[in] time_avr 平均時刻
   * @return true:出力成功 false:出力失敗
   */
  virtual CIO::E_CIO_ERRORCODE
  write_averaged(FILE* fp,
                 const unsigned step_avr,
                 const double time_avr)=0;


//FEAST 20131125.s
  /**
   * @brief ascii ヘッダーレコード出力(bov,avs)
   * @param [in] step step番号
   * @param [in] time time
   */
  virtual
  bool 
  write_ascii_header(const unsigned step,
                     const double time)
  { return true; }; 

  /**
   * @brief データタイプ毎のサイズを取得
   * @param [in] Dtype データタイプ(Int8,Int16,,,,etc)
   * @return データサイズ
   * @return 0 エラー
   */
  static int 
  get_cio_Datasize(CIO::E_CIO_DTYPE Dtype); 

  /**
   * @brief Create Process 
   * @param [in] comm           MPIコミュニケータ
   * @param [out] G_Process     Process class　　　
   */
  void 
  cio_Create_dfiProcessInfo(const MPI_Comm comm,
                            cio_Process &G_Process);


  /**
   * @brief 読込み判定判定
   * @param [in] G_voxel            計算空間全体のボクセルサイズ（自）
   * @param [in] DFI_GlobalVoxel    計算空間全体のボクセルサイズ（DFI）
   * @param [in] G_Div              分割数（自）
   * @param [in] DFI_GlobalDivision 分割数（DFI）
   * @return 読込みタイプコード
   */
  //cio_EGlobalVoxel CheckGlobalVoxel(const int Gvoxel[3], 
  CIO::E_CIO_READTYPE 
  CheckReadType(const int G_voxel[3], 
                const int DFI_GlobalVoxel[3],
                const int G_Div[3],
                const int DFI_GlobalDivision[3]); 

  /**
   * @brief フィールドデータの読込み範囲を求める
   * @param [in] isSame   粗密フラグ true:密、false:粗
   * @param [in] head     計算領域の開始位置(自)　
   * @param [in] tail     計算領域の終了位置(自)　
   * @param [in] gc       仮想セル数(自)　
   * @param [in] DFI_head 計算領域の開始位置(DFI)　　
   * @param [in] DFI_tail 計算領域の終了位置(DFI)　　
   * @param [in] DFI_gc   仮想セル数(DFI)　
   * @param [in] readflag 読込み方法
   * @param [out] copy_sta コピー開始位置
   * @param [out] copy_end コピー終了位置　　
   * @param [out] read_sta 読込み開始位置
   * @param [out] read_end 読込み終了位置　　
   */
  void 
  CreateReadStartEnd(bool isSame,
                     const int head[3], 
                     const int tail[3], 
                     const int gc, 
                     const int DFI_head[3], 
                     const int DFI_tail[3],
                     const int DFI_gc, 
                     const CIO::E_CIO_READTYPE readflag, 
                     int copy_sta[3], 
                     int copy_end[3],
                     int read_sta[3],
                     int read_end[3]);

  /**
   * @brief index DFIファイル出力
   * @param [in] dfi_name  DFIファイル名
   * @return true:出力成功 false:出力失敗
   */
  CIO::E_CIO_ERRORCODE
  WriteIndexDfiFile(const std::string dfi_name);


public:

  /**
   * @brief セル中心データを格子点に値をセット
   * @param [out]  P 格子点データ
   * @param [in]   S セル中心data
   */
  template<class T1, class T2>
  bool setGridData(
                   cio_TypeArray<T1>* P,
                   cio_TypeArray<T2>* S);

  /**
   * @brief 内部の格子点のデータを重み付けでで割る
   * @param[out] P  格子点data
   */
  template<class T>
  void VolumeDataDivide(cio_TypeArray<T> *P);


public:

  /**
   * @brief ディレクトリパスの作成(MakeDirectorySubを呼出して作成)
   * @param [in] path パス
   * @return error code　　　　　　　
   */ 
  int MakeDirectory(const std::string path);

  /**
   * @brief ディレクトリパスの作成(MakeDirectory関数を呼出して作成)
   * @return error code　　　　　　　
   */ 
  int MakeDirectoryPath();

  /**
   * @brief ディレクトリパスの作成(system関数mkdirで作成)
   * @param [in] path パス
   * @return error code　　　　　　　
   */ 
  static int MakeDirectorySub( std::string path );

  /**
   * @brief dfiのパスとDirectoryPathを連結する関数
   * @return パス名
   */
  std::string Generate_Directory_Path(); 

  /** バージョンを出力する
   */
  static std::string getVersionInfo()
  {
    std::string str(CIO_VERSION_NO);
    return str;
  }
  
};

//inline 関数
#include "inline/cio_DFI_inline.h"


#endif // _cio_DFI_H_
