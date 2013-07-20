#ifndef _CIO_DFI_H_
#define _CIO_DFI_H_

/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI.h
 * @brief  cio_DFI Class Header
 * @author kero    
 */
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <typeinfo>
#include <set>
#include <map>

#include "cio_PathUtil.h"
#include "cio_TextParser.h"
#include "cio_Define.h"
#include "cio_ActiveSubDomain.h"

#include "cio_endianUtil.h"

using namespace std;

/** CIO main class */
class cio_DFI {


public:

  enum E_CIO_FORMAT
  {
    E_CIO_FMT_UNKNOWN = -1,
    E_CIO_FMT_SPH,
    E_CIO_FMT_BOV
  };

  enum E_CIO_ONOFF
  {
    E_CIO_OFF = 0,
    E_CIO_ON
  };

  enum E_CIO_DTYPE
  {
    E_CIO_DUMMY=0,
    E_CIO_INT8,
    E_CIO_INT16,
    E_CIO_INT32,
    E_CIO_INT64,
    E_CIO_UINT8,
    E_CIO_UINT16,
    E_CIO_UINT32,
    E_CIO_UINT64,
    E_CIO_FLOAT32,
    E_CIO_FLOAT64
  };

  enum E_CIO_ARRAYSHAPE
  {
    E_CIO_IJKN=0,
    E_CIO_NIJK
  };

  MPI_Comm m_comm;
  std::string m_directoryPath;
  std::string m_indexDfiName;
  std::string m_timeSliceDir;
  bool        m_outSlice;
  int         m_dfi_mng;
  int         m_start_type;
  int         m_RankID;

  /** index.dfi ファイルの FileInfo */
  struct cio_FileInfo
  {
    string DirectoryPath;                 ///<ディレクトリパス
    string TimeSliceDir;                  ///<TimeSlice on or off
    string Prefix;                        ///<ファイル接頭文字
    E_CIO_FORMAT FileFormat;              ///<ファイルフォーマット "bov","sph",,,
    int    GuideCell;                     ///<仮想セルの数
    string DataType;                      ///<配列のデータタイプ "float",,,,
    string Endian;                        ///<エンディアンタイプ "big","little"
    string ArrayShape;                    ///<配列形状
    int    Component;                     ///<成分数

    cio_FileInfo() 
    {
      DirectoryPath="";
      TimeSliceDir ="";
      Prefix       ="";
      FileFormat   = E_CIO_FMT_UNKNOWN;
      GuideCell    =0;
      DataType     ="";
      Endian       ="";
      ArrayShape   ="";
      Component    =0;
    }

    cio_FileInfo(string _DirectoryPath, string _TimeSliceDir, string _Prefix, 
             E_CIO_FORMAT _FileFormat,
             int _GuideCell, string _DataType, string _Endian, 
             string _ArrayShape, int _Component)
    {
      DirectoryPath=_DirectoryPath;
      Prefix       =_Prefix;
      TimeSliceDir =_TimeSliceDir;
      FileFormat   =_FileFormat;
      GuideCell    =_GuideCell;
      DataType     =_DataType;
      Endian       =_Endian;
      ArrayShape   =_ArrayShape;
      Component    =_Component;
    }   

  };


  /** index.dfi ファイルの FilePath */
  struct cio_FilePath
  {
    string Process;                       ///<proc.dfi ファイル名

    cio_FilePath()
    {
      Process="";
    }

    cio_FilePath(string _Process)
    {
      Process=_Process;
    }

  };


  /** index.dfi ファイルの Unit */
  struct cio_Unit  
  {
    bool out_Length;                   ///<Length,L0 出力フラグ
    string Length;                     ///<(NonDimensional, m, cm, mm)
    double L0;                         ///<規格化に用いた長さスケール

    bool out_Velocity;                 ///<Velocity,V0 出力フラグ
    string Velocity;                   ///<(NonDimensional, m/s)
    double V0;                         ///<代表速度(m/s)

    bool out_Pressure;                 ///<Presuure,P0,DiffPrs出力フラグ
    string Pressure;                   ///<(NonDimensional, Pa)
    double P0;                         ///<基準圧力(Pa)
    double DiffPrs;                    ///<圧力差(Pa)

    bool out_Temperature;              ///<Temperature,BaseTemp,DiffTemp出力フラグ
    string Temperature;                ///<(NonDimensional, C, K)
    double BaseTemp;                   ///<指定単力　　　　
    double DiffTemp;                   ///<指定単位　　　

    cio_Unit()
    {
      out_Length=false;
      Length    ="";
      L0        =0.0;

      out_Velocity=false;
      Velocity    ="";
      V0          =0.0;

      out_Pressure=false;
      Pressure    ="";
      P0          =0.0;
      DiffPrs     =0.0;

      out_Temperature=false;
      Temperature    ="";
      BaseTemp       =0.0;
      DiffTemp       =0.0;
    }

    cio_Unit(bool _out_Length, string _Length, double _L0,
         bool _out_Velocity, string _Velocity, double _V0,
         bool _out_Pressure, string _Pressure, double _P0, double _DiffPrs,
         bool _out_Temperature, string _Temperature, double _BaseTemp, double _DiffTemp)
    {

      out_Length=_out_Length;
      Length    =_Length;
      L0        =_L0;

      out_Velocity=_out_Velocity;
      Velocity    =_Velocity;
      V0          =_V0;

      out_Pressure=_out_Pressure;
      Pressure    =_Pressure;
      P0          =_P0;
      DiffPrs     =_DiffPrs;

      out_Temperature = _out_Temperature;
      Temperature     = _Temperature;
      BaseTemp        = _BaseTemp;
      DiffTemp        = _DiffTemp;
    }
  };

  /** index.dfi ファイルの Slice */
  struct cio_Slice
  {
    int    step;                          ///<ステップ番号
    double time;                          ///<時刻
    int    AveragedStep;
    double AveragedTime;
    vector<double> Min;                   ///<最小値
    vector<double> Max;                   ///<最大値
  };
  
  /** proc.dfi ファイルの Domain */
  struct cio_Domain
  {
    double GlobalOrigin[3];             ///<計算空間の起点座標
    double GlobalRegion[3];             ///<計算空間の各軸方向の長さ
    int GlobalVoxel[3];                    ///<計算領域全体のボクセル数
    int GlobalDivision[3];                 ///<計算領域の分割数
    string ActiveSubdomain;                ///<ActiveSubdomainファイル名

    cio_Domain()
    {
      for(int i=0; i<3; i++) GlobalOrigin[i]=0.0;
      for(int i=0; i<3; i++) GlobalRegion[i]=0.0;
      for(int i=0; i<3; i++) GlobalVoxel[i]=0;
      for(int i=0; i<3; i++) GlobalDivision[i]=0;
      ActiveSubdomain="";
    }

    cio_Domain(double* _GlobalOrigin, double* _GlobalRegion, int* _GlobalVoxel, 
           int* _GlobalDivision)
    {
      GlobalOrigin[0]=_GlobalOrigin[0];
      GlobalOrigin[1]=_GlobalOrigin[1];
      GlobalOrigin[2]=_GlobalOrigin[2];

      GlobalRegion[0]=_GlobalRegion[0];
      GlobalRegion[1]=_GlobalRegion[1];
      GlobalRegion[2]=_GlobalRegion[2];

      GlobalVoxel[0]=_GlobalVoxel[0];
      GlobalVoxel[1]=_GlobalVoxel[1];
      GlobalVoxel[2]=_GlobalVoxel[2];

      GlobalDivision[0]=_GlobalDivision[0];
      GlobalDivision[1]=_GlobalDivision[1];
      GlobalDivision[2]=_GlobalDivision[2];
    }

  };

  /** proc.dfi ファイルの MPI */
  struct cio_MPI
  {
    int NumberOfRank;                      ///<プロセス数
    int NumberOfGroup;                     ///<グループ数

    cio_MPI()
    {
       NumberOfRank=0;
       NumberOfGroup=1;
    }

    cio_MPI(int _NumberOfRank)
    {
       NumberOfRank=_NumberOfRank;
    }

  };

  /** proc.dfi ファイルの Process */
  struct cio_Rank
  {
    int RankID;                           ///<ランク番号
    string HostName;                      ///<ホスト名 
    int VoxelSize[3];                     ///<ボクセルサイズ
    int HeadIndex[3];                     ///<始点インデックス
    int TailIndex[3];                     ///<終点インデックス
  };

  cio_FileInfo DFI_Finfo;
  cio_FilePath DFI_Fpath;
  cio_Unit     DFI_Unit;
  cio_Domain   DFI_Domain;
  cio_MPI      DFI_MPI;
  vector<cio_Slice> TimeSlice;
  vector<cio_Rank> RankInfo;

  std::vector<cio_ActiveSubDomain> CIO_subDomainInfo; ///< CIO 活性サブドメイン情報
  int *CIO_rankMap;
  int (*CIO_HeadTail)[6];

  std::vector<cio_ActiveSubDomain> DFI_subDomainInfo; ///< DFI 活性サブドメイン情報
  int *DFI_rankMap;
  int (*DFI_HeadTail)[6];

  std::set<int>DFI_headX,DFI_headY,DFI_headZ;

  headT DFI_mapX,DFI_mapY,DFI_mapZ;


public:
  /** コンストラクタ */
  cio_DFI();
  
  /**　デストラクタ */
  ~cio_DFI();

  /**
   * @brief read インスタンス
   * @param [in] comm    MPIコミュニケータ
   * @param [in] dfifile DFIファイル名
   * @return インスタンスされたクラスのポインタ
   */
  static cio_DFI* ReadInit(MPI_Comm comm, string dfifile); 

  /**
   * ActiveSubdomainファイルのエンディアンをチェック
   * @param[in]  ident               ActiveSubdomainファイルのIdentifier
   * @retval     CIO_Match   一致
   * @retval     CIO_UnMatch 不一致
   * @retval     CIO_UnKnown フォーマットが異なる
   */
  static cio_EMatchType isMatchEndianSbdmMagick( int ident ); 

  /**
   * ActiveSubdomainファイルの読み込み(static関数)
   * @param[in]  subDomainFile ActiveSubdomainファイル名
   * @param[out] subDomainInfo 活性ドメイン情報
   * @param[out] div           ActiveSubdiomainファイル中の領域分割数
   * @return   終了コード(CIO_SUCCESS=正常終了)
   */
  static cio_ErrorCode ReadActiveSubdomainFile( std::string subDomainFile,
         std::vector<cio_ActiveSubDomain>& subDomainInfo, int div[3] );

  /**
   * 活性サブドメイン配列が空のとき、全領域が活性サブドメインになるため
   * このチェック関数内で活性サブドメイン情報を生成する.
   * @param[in] nRank 並列プロセス数
   * @param[in] div 領域分割数
   * @param[out] subDomainInfo 活性ドメイン情報
   * @return   終了コード(CIO_SUCCESS=正常終了)
   */ 
   static cio_ErrorCode CheckData( int nRank, int div[3],
             std::vector<cio_ActiveSubDomain>& subDomainInfo );

  /**
   * ランクマップを生成（非活性を含む）
   * @param[in] div 領域分割数
   * @param[in] subDomainInfo 活性ドメイン情報
   * @retval ランクマップ
   * @retval NULL
   */
   static int* CreateRankMap(int div[3],std::vector<cio_ActiveSubDomain> subDomainInfo); 

  /**
   *
   */
   static int* CreateActiveRankMap(vector<cio_Rank> &RankInfo, int (*HeadTail)[6], 
                                   headT mapX, headT mapY,headT mapZ);

  /**
   * head&tail情報の登録
   * @retval true  正常終了
   * @retval false エラー
   */
   static bool SetHeadTail( vector<cio_Rank> RankInfo, int (*HeadTail)[6] ,
                               std::set<int>&headx, std::set<int>&heady, std::set<int>&headz );
  
  /**
   * head&tail情報の作成
   * @retval true  正常終了
   * @retval false エラー
   */
   static bool CreateHeadTail( int div[3], int* rankMap, 
                               int gvox[3] ,int (*HeadTail)[6],
                               vector<cio_Rank> &RankInfo, cio_Rank rank);

  /**
   *
   */
   static bool CreateHeadMap( std::set<int>head, headT &map); 
 
  /**
   * @brief write インスタンス
   * @param [in] comm          MPIコミュニケータ
   * @param [in] DfiName       DFIファイル名
   * @param [in] DirectoryPath フィールドデータのディレクトリ
   * @param [in] Prefix        ベースファイル名
   * @param [in] FileFormat    ファイルフォーマット
   * @param [in] GuideCell     出力仮想セル数　　　
   * @param [in] DataType      データタイプ　　　　
   * @param [in] ArrayShape    配列形状　　　　　　
   * @param [in] Component     成分数　　　　　　　
   * @param [in] Process       proc.dfiファイル名
   * @param [in] G_size[3]     グローバルボクセルサイズ　
   * @param [in] pitch[3]      ピッチ　　　　　　　　　　
   * @param [in] G_origin[3]   原点座標値　　　　　　　　
   * @param [in] division[3]   領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [in] hostname      ホスト名　　　　　　　　　
   * @return インスタンスされたクラスのポインタ
   */
  template<class T>
  static cio_DFI* WriteInit(MPI_Comm comm, string DfiName, string Path, string prefix,
                            E_CIO_FORMAT format, int GCell, E_CIO_DTYPE DataType, 
                            E_CIO_ARRAYSHAPE ArrayShape,
                            int Comp, string process, int G_size[3], T pitch[3],
                            T G_origin[3], int division[3], int head[3], int tail[3],
                            string hostname, E_CIO_ONOFF TSliceOnOff)
  {

  cio_DFI *dfi = NULL;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  cio_FileInfo out_F_info;
  out_F_info.DirectoryPath = Path;
  if( TSliceOnOff == E_CIO_ON ) {
    out_F_info.TimeSliceDir  = D_CIO_ON;
  } else {
    out_F_info.TimeSliceDir  = D_CIO_OFF;
  }
  out_F_info.Prefix        = prefix;
  out_F_info.FileFormat    = format;
  out_F_info.GuideCell     = GCell;
  if     ( DataType == E_CIO_INT8    ) out_F_info.DataType = D_CIO_INT8;
  else if( DataType == E_CIO_INT16   ) out_F_info.DataType = D_CIO_INT16;
  else if( DataType == E_CIO_INT32   ) out_F_info.DataType = D_CIO_INT32;
  else if( DataType == E_CIO_INT64   ) out_F_info.DataType = D_CIO_INT64;
  else if( DataType == E_CIO_UINT8   ) out_F_info.DataType = D_CIO_UINT8;
  else if( DataType == E_CIO_UINT16  ) out_F_info.DataType = D_CIO_UINT16;
  else if( DataType == E_CIO_UINT32  ) out_F_info.DataType = D_CIO_UINT32;
  else if( DataType == E_CIO_UINT64  ) out_F_info.DataType = D_CIO_UINT64;
  else if( DataType == E_CIO_FLOAT32 ) out_F_info.DataType = D_CIO_FLOAT32;
  else if( DataType == E_CIO_FLOAT64 ) out_F_info.DataType = D_CIO_FLOAT64;
  if(      ArrayShape == E_CIO_IJKN ) out_F_info.ArrayShape = D_CIO_IJNK;
  else if( ArrayShape == E_CIO_NIJK ) out_F_info.ArrayShape = D_CIO_NIJK;
  out_F_info.Component     = Comp;

  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  if( cdumy[0] == 0x01 ) out_F_info.Endian = "little";
  if( cdumy[0] == 0x00 ) out_F_info.Endian = "big";

  cio_FilePath out_F_path;
  out_F_path.Process = process;

  cio_Unit out_unit;

  cio_MPI out_mpi;
  out_mpi.NumberOfRank = nrank;
  out_mpi.NumberOfGroup = 1;

  cio_Domain out_domain;
  vector<cio_Rank> out_RankInfo;
  cio_Rank out_Rank;

  for(int i=0; i<nrank; i++ ) {
     out_RankInfo.push_back(out_Rank);
  }      

  out_RankInfo[RankID].RankID=RankID;
  for(int i=0; i<3; i++) {
    out_RankInfo[RankID].HeadIndex[i]=head[i];
    out_RankInfo[RankID].TailIndex[i]=tail[i];
    out_RankInfo[RankID].VoxelSize[i]=tail[i]-head[i]+1;
  }

  //cio_Create_Domain(comm, G_size, division, head, tail, out_domain, out_RankInfo, out_Rank);

  for(int i=0; i<3; i++) {
    out_domain.GlobalVoxel[i]  = G_size[i];
    out_domain.GlobalDivision[i] = division[i];
    out_domain.GlobalOrigin[i] = (double)G_origin[i];
    out_domain.GlobalRegion[i] = pitch[i]*G_size[i];
  }

  vector<cio_Slice> out_TSlice;

  char tmpname[512];
  memset(tmpname,0x00,sizeof(char)*512);
  if( gethostname(tmpname, 512) != 0 ) printf("*** error gethostname() \n");

  dfi = get_dfi(out_F_info, out_F_path, out_unit, out_domain, out_mpi,
                          out_TSlice, out_RankInfo);

  if( dfi == NULL ) return NULL;

  dfi->m_indexDfiName = DfiName;
  dfi->m_directoryPath = CIO::cioPath_DirName(DfiName);
  dfi->m_comm = comm;
  dfi->m_RankID = RankID;

  return dfi;

  };


  /**
   * @brief field data  ファイル名の作成
   * @param [in] RankID ランク番号
   * @param [in] step   読込みステップ番号
   * @param [in] mio    並列判定フラグ（逐次or並列の判定用）
   * @return 生成されたファイル名　　　　　　　
   */
  std::string Generate_FileName(int RankID,int step, const bool mio);

  /**
   * @brief ディレクトリパスの作成
   * @param [in] path パス
   * @return error code　　　　　　　
   */ 
  static int MakeDirectory(std::string path);

  /**
   * @brief ディレクトリパスの作成
   * @return error code　　　　　　　
   */ 
  int MakeDirectoryPath();

  /**
   * @brief initialise dfi
   */ 
  void InitDFI();

  /**
   * @brief read FileInfo(inde.dfi)
   * @param [in]   dfifile index.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  finfo   読込んだcio_FileInfo 
   * @return error code
   */
  static int readFileInfo(string dfifile, cio_TextParser tpCntl, cio_FileInfo &finfo);

  /**
   * @brief read FilePath(inde.dfi)
   * @param [in]   dfifile index.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  fpath   読込んだFilePath 
   * @return error code
   */
  static int readFilePath(string dfifile, cio_TextParser tpCntl, cio_FilePath &fpath);

  /**
   * @brief read Unit(inde.dfi)
   * @param [in]   dfifile index.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  unit    読込んだUnit 
   * @return error code
   */
  static int readUnit(string dfifile, cio_TextParser tpCntl, cio_Unit &unit);

  /**
   * @brief read Slice(inde.dfi)
   * @param [in]      dfifile   index.dfiファイル名
   * @param [in]      tpCntl    cio_TextParserクラス 
   * @param [out]     TimeSlice 読込んだSliceを格納した領域 
   * @param [out,out] slice     TimeSlice読込み用領域 
   * @return error code
   */
  static int readSlice(string dfifile, cio_TextParser tpCntl, vector<cio_Slice> &TimeSlice, cio_Slice  slice);

  /**
   * @brief read Domain(proc.dfi)
   * @param [in]   dfifile proc.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  domain  読込んだDomain
   * @return error code
   */
  static int readDomain(string dfifile, cio_TextParser tpCntl, cio_Domain &domain);

  /**
   * @brief read MPI(proc.dfi)
   * @param [in]   dfifile proc.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  mpi     読込んだMPI
   * @return error code
   */
  static int readMPI(string dfifile, cio_TextParser tpCntl, cio_Domain domain, cio_MPI &mpi);

  /**
   * @brief read Rank(proc.dfi)
   * @param [in]   dfifile proc.dfiファイル名
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @param [out]  rank    読込んだProcess
   * @return error code
   */
  static int readRank(string dfifile, cio_TextParser tpCntl, vector<cio_Rank> &RankInfo, cio_Rank rank);

  /**
   * @brief get ArrayShape （配列形状の取り出し関数）
   * @return 配列形状
   */
  std::string getArrayShape();

  /**
   * @brief get DataType （データタイプの取り出し関数）
   * @return データタイプ
   */
  std::string getDataType();

  /**
   * @brief get Component （成分数の取り出し関数）
   * @return 成分数
   */
  int getComponent();

  /**
   * @brief get E_CIO_DTYPE （データタイプの取り出し関数） 
   * @param[in] datatype dfiから取得したデータタイプ
   * @return データタイプ番号
   */
  E_CIO_DTYPE get_cio_Datatype(string datatype); 

  /**
   *
   *
   */
  int get_cio_Datasize(E_CIO_DTYPE Dtype); 
  /**
   * @brief Create Domain & Process 
   * @param [in] comm          MPIコミュニケータ
   * @param [in] G_voxel[3]    グローバルボクセルサイズ　
   * @param [in] G_division[3] 領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [out]G_domain      Domain情報(構造体)　　　　
   * @param [out]G_RankInfo    Process情報(vector)　　　
   * @param [in] G_Rank        Process情報(構造体)　　　
   */
  static void cio_Create_Domain(MPI_Comm comm,
                                int G_voxel[3], int G_division[3],
                                int head[3], int tail[3],
                                cio_Domain &G_domain, vector<cio_Rank> &G_RankInfo, 
                                cio_Rank G_Rank);


  /**
   *
   *
   */
  static cio_DFI* get_dfi(cio_FileInfo F_Info, cio_FilePath F_Path, cio_Unit unit, 
                          cio_Domain domain,
                          cio_MPI mpi,vector<cio_Slice> TSlice, vector<cio_Rank> RInfo);  
  /**
   * @brief read field data record
   * @param [in] step 入力ステップ番号
   * @param [in] gc   仮想セル数　　　
   * @param [in] Gvoxel[3]    グローバルボクセルサイズ　
   * @param [in] Gdivision[3] 領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [out]val           読み込んだデータポインタ　
   *
   */
  virtual void ReadData(int step, int gc, 
                        int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                        void *val, double &time ,
                        const bool mode, unsigned &step_avr, double &time_avr) = 0; 

  virtual void* ReadData(int step, int gc, 
                        int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                        double &time,
                        const bool mode, unsigned &step_avr, double &time_avr ) = 0;

  /**
   * @brief 粗密データ判定
   * @param [in] Gvoxel     計算空間全体のボクセルサイズ（自）
   * @param [in] DFI_Gvoxel 計算空間全体のボクセルサイズ（DFI）
   * @return CIO_E_GV_SAME:密 CIO_E_GVX2_SAME:粗 CIO_E_OTHER:その他
   */
  cio_EGlobalVoxel CheckGlobalVoxel(int Gvoxel[3], int DFI_Gvoxel[3]); 

  /**
   * @brief 粗密ファイルのMxN判定
   * @param [in] rankList 読込みに必要なランク番号リスト
   * @param [in] head[3]  計算領域の開始位置　　　　
   * @param [in] tail[3]  計算領域の終了位置　　　　
   * @param [in] gc       仮想セル数(自)　　　
   * @param [in] dfi_gc   仮想セル数(DFI)　　　
   * @return  true:密 false:粗
   */
  bool CheckMxN(vector<int> &rankList, int head[3], int tail[3], int gc, int dfi_gc);
 
  /**
   * @brief 読込みランクファイルリストの作成
   * @param [in] head[3]  計算領域の開始位置　　　　
   * @param [in] tail[3]  計算領域の終了位置　　　　
   * @param [in] gc       仮想セル数(自)　　　
   * @param [in] readflag 粗密データ判定フラグ
   * @param [out]rankList 読込みに必要なランク番号リスト 
   */
  void CreateRankList(int head[3], int tail[3], int gc, cio_EGlobalVoxel readflag,
                      vector<int> &rankList);

  /**
   * @brief 読込み範囲を求める
   * @param [in] head[3]     計算領域の開始位置(自)　
   * @param [in] tail[3]     計算領域の終了位置(自)　
   * @param [in] gc          仮想セル数(自)　
   * @param [in] DEF_head[3] 計算領域の開始位置(DFI)　　
   * @param [in] DEF_tail[3] 計算領域の終了位置(DFI)　　
   * @param [in] DFI_gc      仮想セル数(DFI)　
   * @param [out]sta[3]      読込み開始位置
   * @param [out]end[3]      読込み終了位置　　
   * @return true:１対１
   */
  bool CheckReadArea(int head[3], int tail[3], int gc, int DEF_head[3], int DFI_tail[3],
                     int DFI_gc, cio_EGlobalVoxel readflag, int sta[3], int end[3]);


  /**
   * @brief write field data record
   * @param [in] step     出力ステップ番号
   * @param [in] gc       仮想セル数　　　
   * @param [in] tiem     出力時刻　　　　
   * @param [in] val      出力データポインタ
   * @param [in] interval 出力間隔
   * @param [in] force    強制出力指示
   */ 
  virtual void WriteData(int step, int gc, void* time, 
                        void *val, void *minmax, int interval,
                        const bool mode, const unsigned step_avr, const double time_avr,
                        bool force) = 0;

  /**
   * @brief index DFIファイル出力コントロール
   * @param [in] prefix  ファイル接頭文字
   * @param [in] step    ステップ
   * @param [in] time    時間　　
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
  bool WriteIndexDfiFile(string DfiName, int RabkID,const std::string prefix, const unsigned step, 
       void* time, void *minmax, const bool mio, const bool mode, const unsigned step_avr, 
       const double time_avr); 

  /**
   * @brief proc DFIファイル出力コントロール
   * @param [in] comm          MPIコミュニケータ
   * @param [in] procFileName  出力proc.dfiファイル名
   * @param [in] G_size[3]     グローバルボクセルサイズ　
   * @param [in] division[3]   領域分割数　　　　　　　　
   * @param [in] head[3]       計算領域の開始位置　　　　
   * @param [in] tail[3]       計算領域の終了位置　　　　
   * @param [in] hostname      ホスト名　　　　　　　　　
   * @param [in] out_host      ホスト名出力フラグ　　　　
   * @return true:出力成功 false:出力失敗
   */
  template<class T>
  static bool WriteProcDfiFile(MPI_Comm comm, string procFileName, 
                               int G_size[3],
                               int division[3], int head[3], int tail[3], T org[3],
                               T pch[3], 
                               string hostname, bool out_host)
  {

  if( procFileName.empty() ) return false;

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  cio_MPI out_mpi;
  out_mpi.NumberOfRank = nrank;
  out_mpi.NumberOfGroup = 1;

  cio_Domain out_domain;
  vector<cio_Rank> out_RankInfo;
  cio_Rank out_Rank;

  cio_Create_Domain(comm, G_size, division, head, tail, out_domain, out_RankInfo, out_Rank);
  for(int i=0; i<3; i++) {
    out_domain.GlobalOrigin[i] = org[i];
    out_domain.GlobalRegion[i] = pch[i]*G_size[i];
  }

  if( out_host ) {
    const int LEN=256;
    char *recbuf = new char[out_RankInfo.size()*LEN];
    char  sedbuf[LEN];
    sprintf(sedbuf,"%s",hostname.c_str());
    MPI_Gather(sedbuf,LEN,MPI_CHAR,recbuf,LEN,MPI_CHAR,0,MPI_COMM_WORLD);

    for( int i=0; i<out_RankInfo.size(); i++ ) {
     char* hn =&(recbuf[i*LEN]);
     out_RankInfo[i].HostName=(string(hn));
    }

    if( recbuf ) delete [] recbuf;
  }

  if(RankID != 0) return NULL;

  if( !Write_Proc_File(procFileName,out_domain,out_mpi,out_RankInfo) )
  {
    return false;
  }

  return true;

  };

  /**
   * @brief 出力DFIファイル名を作成する
   * @param [in] prefix ファイル接頭文字
   * @return DFIファイル名
   */ 
  static std::string Generate_DFI_Name(const std::string prefix);

  /**
   * @brief Directoryパスを生成する
   * @return パス名
   */
  std::string Generate_Directory_Path(); 


  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name  DFIファイル名
   * @param [in] prefix    ファイル接頭文字
   * @param [in] step      ステップ数
   * @param [in] time      時間　　
   * @param [in] dfi_mng   出力管理カウンタ
   * @param [in] mio       出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Index_File(const std::string dfi_name, const std::string prefix, const unsigned step, 
                        void* time, int& dfi_mng, void *minmax, const bool mio,
                        const bool avr_mode, const unsigned step_avr,const double time_avr);

  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name  DFIファイル名
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Proc_File(const std::string dfi_name); 

  /**
   * @brief DFIファイルを出力する
   * @param [in] dfi_name   DFIファイル名
   * @param [in] out_domain 出力Domain
   * @param [in] out_mpi    出力MPI
   * @param [in] out_RankInfo 出力Process
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Proc_File(const std::string dfi_name, cio_Domain out_domain, 
                              cio_MPI out_mpi, vector<cio_Rank> out_RankInfo); 

  /**
   * @brief DFIファイル:FileInfo要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   * @return true:出力成功 false:出力失敗
   */
   bool Write_FileInfo(FILE* fp, const unsigned tab, const std::string prefix);

  /**
   * @brief DFIファイル:Unit要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   * @return true:出力成功 false:出力失敗
   */
   bool Write_Unit(FILE* fp, const unsigned tab, const std::string prefix);

  /**
   * @brief DFIファイル:TimeSlice要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] step    ステップ数
   * @param [in] time    時間　　
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
   template<class T>
   bool Write_TimeSlice(FILE* fp, const unsigned tab, const unsigned step, T* time,
                        T* minmax, const bool avr_mode, const unsigned a_step, const double a_time)
   {

     string compname;

     Write_Tab(fp, tab);
     fprintf(fp, "Slice[@] {\n");

     Write_Step(fp,tab+1,step);

     Write_Time(fp,tab+1,time);

     if( !avr_mode ) {
       Write_Average(fp,tab+1,a_step,a_time);
     }

     if( DFI_Finfo.Component ) {
       Write_Tab(fp, tab+1);
       fprintf(fp, "MinMax[@] {\n");
       for(int i=0; i<1; i++){
         compname="Min";
         Write_Comp(fp,tab+2,compname,minmax[i*2]);
         compname="Max";
         Write_Comp(fp,tab+2,compname,minmax[i*2+1]);
       }
       Write_Tab(fp, tab+1);
       fprintf(fp, "}\n");
     }

     Write_Tab(fp, tab);
     fprintf(fp, "}\n");

     return true;
   };

  /**
   * @brief Tab(space２つ)を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント数
   */
  void Write_Tab(FILE* fp, const unsigned tab);


  /**
   * @brief DFIファイル:出力ファイル情報要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] prefix  ファイル接頭文字
   * @param [in] tab     インデント
   * @param [in] step    ステップ数
   * @param [in] mio     出力時の分割指定　 true = local / false = gather
   * @return true:出力成功 false:出力失敗
   */
  bool Write_OutFileInfo(FILE* fp, const unsigned tab, const std::string prefix, const unsigned step, const bool mio);


  /**
   * @brief DFIファイル:DirectoryPathを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] dirpath ディレクトリ　　
   */
 void Write_DirectoryPath(FILE* fp, const unsigned tab, const std::string dirpath); 

  /**
   * @brief DFIファイル:TimeSliceDirectoryを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] dirpath ディレクトリ　　
   */
 void Write_TimeSliceDir(FILE* fp, const unsigned tab, const std::string timeslicedir); 

  /**
   * @brief DFIファイル:BaseName要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] prefix  ファイル接頭文字
   */
 void Write_BaseName(FILE* fp, const unsigned tab, const std::string prefix); 

  /**
   * @brief DFIファイル:ファイルフォーマット要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_FileFormat(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:ガイドセル要素を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_GuideCell(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:データタイプを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_DataType(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Endianを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Endian(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:ArrayShapeを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_ArrayShape(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Componentを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Component(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Process(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Lengthを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Length(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:L0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_L0(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Velocityを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Velocity(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:V0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_V0(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Pressureを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Pressure(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:P0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_P0(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:P0を出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_DiffPrs(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Temperatureを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Temperature(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:BaseTempを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_BaseTemp(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:DiffTempを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_DiffTemp(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Stepを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] step   step番号
   */
  void Write_Step(FILE* fp, const unsigned tab, const unsigned step);

  /**
   * @brief DFIファイル:Timeを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] time   time 
   */
  template<class T>
  void Write_Time(FILE* fp, const unsigned tab, T time)
  {
    Write_Tab(fp, tab);
    fprintf(fp, "Time = %e\n",time[0]);
  };

  /**
   * @brief DFIファイル:Averageを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] a_step
   * @param [in] a_time
   */
  void Write_Average(FILE* fp, const unsigned tab, const unsigned a_step, const double a_time);

  /**
   * @brief DFIファイル:成分を出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @param [in] compname 成分名
   * @param [in] comp     出力成分 
   */
  template<class T>
  void Write_Comp(FILE* fp, const unsigned tab, const std::string compname,
                  T comp)
  {
    Write_Tab(fp, tab);
    fprintf(fp, "%s = %e\n",compname.c_str(),comp);
  }
  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Domain(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Domainを出力する
   * @param [in] fp         ファイルポインタ
   * @param [in] tab        インデント
   * @param [in] out_domain 出力Domain
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Domain(FILE* fp, const unsigned tab, cio_Domain out_domain);

  /**
   * @brief DFIファイル:MPIを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_MPI(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:MPIを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   * @param [in] out_mpi 出力MPI
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_MPI(FILE* fp, const unsigned tab, cio_MPI out_mpi);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Process_Rank(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] out_RankInfo 出力Process
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Process_Rank(FILE* fp, const unsigned tabi, vector<cio_Rank> out_RankInfo);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @return true:出力成功 false:出力失敗
   */
  bool Write_Rank(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:Processを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] rank   出力Process
   * @return true:出力成功 false:出力失敗
   */
  static bool Write_Rank(FILE* fp, const unsigned tab, cio_Rank rank);

  /**
   * @brief DFIファイル:Originを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Origin(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Originを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   * @param [in] org[3] 出力Origin
   */
  static void Write_Origin(FILE* fp, const unsigned tab, double org[3]);

  /**
   * @brief DFIファイル:Regionを出力する
   * @param [in] fp     ファイルポインタ
   * @param [in] tab    インデント
   */
  void Write_Region(FILE* fp, const unsigned tab);

  /**
   * @brief DFIファイル:Regionを出力する
   * @param [in] fp        ファイルポインタ
   * @param [in] tab       インデント
   * @param [in] Region[3] 出力Region
   */
  static void Write_Region(FILE* fp, const unsigned tab, double Region[3]);

  /**
   * @brief DFIファイル:ノード番号要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
  */
  void Write_MyID(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:ノード数要素を出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   */
  void Write_NodeNum(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:ノード数要素を出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] NumberOfRank ノード数
   */
  static void Write_NodeNum(FILE* fp, const unsigned tab, int NumberOfRank); 

  /**
   * @brief DFIファイル:全体ボクセルサイズ要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_WholeSize(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:全体ボクセルサイズ要素を出力する
   * @param [in] fp     　　    ファイルポインタ
   * @param [in] tab    　　    インデント
   * @param [in] GlobalVoxel[3] 全体ボクセルサイズ
   */
  static void Write_WholeSize(FILE* fp, const unsigned tab, int GlobalVoxel[3]); 

  /**
   * @brief DFIファイル:I,J,K分割数要素を出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_NumDivDomain(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:I,J,K分割数要素を出力する
   * @param [in] fp                ファイルポインタ
   * @param [in] tab               インデント
   * @param [in] GlobalDivision[3] 分割数
   */
  static void Write_NumDivDomain(FILE* fp, const unsigned tab, int GlobalDivision[3]); 

  /**
   * @brief DFIファイル:ActiveSubdomainを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  void Write_ActiveSubdomain_fname(FILE* fp, const unsigned tab); 

  /**
   * @brief DFIファイル:ActiveSubdomainを出力する
   * @param [in] fp      ファイルポインタ
   * @param [in] tab     インデント
   */
  static void Write_ActiveSubdomain_fname(FILE* fp, const unsigned tab, string ActiveSubdomain); 

  /**
   * @brief DFIファイル:IDを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号     
  */
  void Write_ID(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:IDを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号     
  */
  static void Write_RankID(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:Hostnameを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_Hostname(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:Hostnameを出力する
   * @param [in] fp       ファイルポインタ
   * @param [in] tab      インデント
   * @param [in] hostname ホスト名
  */
  static void Write_Hostname(FILE* fp, const unsigned tab, string hostname);

  /**
   * @brief DFIファイル:VoxelSizeを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_L_VoxelSize(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:VoxelSizeを出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] VoxelSize[3] VoxelSize
  */
  static void Write_L_VoxelSize(FILE* fp, const unsigned tab, int VoxelSize[3]);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_HeadIndex(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp　　　　   ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] HeadIndex[3] HeadIndex
  */
  static void Write_HeadIndex(FILE* fp, const unsigned tab, int HeadIndex[3]);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp  ファイルポインタ
   * @param [in] tab インデント
   * @param [in] n   ランク番号
  */
  void Write_TailIndex(FILE* fp, const unsigned tab, const int n);

  /**
   * @brief DFIファイル:HeadIndexを出力する
   * @param [in] fp           ファイルポインタ
   * @param [in] tab          インデント
   * @param [in] TailIndex[3] TailIndex
  */
  static void Write_TailIndex(FILE* fp, const unsigned tab, int TailIndex[3]);

  /**
   *
   */
  void SetUnitLength(bool out_length, string Length, double L0); 

  /**
   *
   */
  void SetUnitVelo(bool out_Velocity, string Velocity, double V0); 

  /**
   *
   */
  void SetUnitPres(bool out_Pressure, string Pressure, double P0, double DiffPrs); 

  /**
   *
   */
  void SetUnitTemp(bool out_Temp, string Temp, double Btemp, double DiffTemp); 

  /**
   * 
   *
   */
  bool dbwrite(int RankID);

  
  /** バージョンを出力する
   */
  static void VersionInfo(std::ostream &ofs)
  {
    ofs << std::endl
    << " CIOlib - Cartesian I/O Library \t\tVersion " << CIO_VERSION_NO << std::endl
    << std::endl;
  }

protected:
  static int MakeDirectorySub( std::string path );
  
};

#endif // _cio_DFI_H_
