#ifndef _CIO_DEFINE_H_
#define _CIO_DEFINE_H_

/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   cio_Define.h
 * @brief  CIOの定義マクロ記述ヘッダーファイル
 * @author kero
 */


#include "mpi.h"
#include "cio_Version.h"


#define D_CIO_EXT_SPH "sph"
#define D_CIO_EXT_BOV "dat"

#define D_CIO_ON  "on"
#define D_CIO_OFF "off"

#define D_CIO_INT8    "Int8"
#define D_CIO_INT16   "Int16"
#define D_CIO_INT32   "Int32"
#define D_CIO_INT64   "Int64"
#define D_CIO_UINT8   "UInt8"
#define D_CIO_UINT16  "UInt16"
#define D_CIO_UINT32  "UInt32"
#define D_CIO_UINT64  "UInt64"
#define D_CIO_FLOAT32 "Float32"
#define D_CIO_FLOAT64 "Float64"

#define D_CIO_IJNK "ijkn"
#define D_CIO_NIJK "nijk"


typedef std::map<int,int> headT;

/** CIOのエラーコード */
enum cio_ErrorCode
{
  CIO_SUCCESS                      = 1    ///<正常終了
, CIO_ERROR                        = -1   ///<エラー終了
, CIO_ERROR_MISMATCH_NP_SUBDOMAIN  = 3003 ///< 並列数とサブドメイン数が一致していない
, CIO_ERROR_OPEN_SBDM              = 3012 ///< ActiveSubdomainファイルのオープンに失敗
, CIO_ERROR_READ_SBDM_HEADER       = 3013 ///< ActiveSubdomainファイルのヘッダー読み込みに失敗
, CIO_ERROR_READ_SBDM_FORMAT       = 3014 ///< ActiveSubdomainファイルのフォーマットエラー
, CIO_ERROR_READ_SBDM_DIV          = 3015 ///< ActiveSubdomainファイルの領域分割数読み込みに失敗
, CIO_ERROR_READ_SBDM_CONTENTS     = 3016 ///< ActiveSubdomainファイルのContents読み込みに失敗
, CIO_ERROR_SBDM_NUMDOMAIN_ZERO    = 3017 ///< ActiveSubdomainファイルの活性ドメイン数が0
, CIO_ERROR_INVALID_DIVNUM         = 3011 ///< 領域分割数が不正
};

/** Endian check */
enum cio_EMatchType
{
  CIO_UnKnown = 0       ///<未定 
, CIO_Match   = 1       ///<一致
, CIO_UnMatch = 2       ///<不一致
};

/** 粗密判定コード */
enum cio_EGlobalVoxel
{
  CIO_E_GV_SAME
, CIO_E_GVX2_SAME
, CIO_E_OTHER
};

/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (long long)(_K+_VC) * (long long)(_NI+2*_VC) * (long long)(_NJ+2*_VC) \
+ (long long)(_J+_VC) * (long long)(_NI+2*_VC) \
+ (long long)(_I+_VC) \
)

/** 2次元インデクス(i,j) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_IJ(_I,_J,_NI,_NJ,_VC) \
( (long long)(_J+_VC) * (long long)(_NI+2*_VC) \
+ (long long)(_I+_VC) \
)

#define _IDX_NIJ(_N,_I,_J,_NI,_NJ,_NN,_VC) \
( (long long)(_NN)*_IDX_IJ(_I,_J,_NI,_NJ,_VC) \
 + (long long)(_N) \
)


/** 4次元インデクス(i,j,k,n) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  成分インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_IJKN(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) \
( (long long)(_N) * (long long)(_NI+2*_VC) * (long long)(_NJ+2*_VC) * (long long)(_NK+2*_VC) \
+ _IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
)

/** 4次元インデクス(n,i,j,k) -> 1次元インデクス変換マクロ
 *  @param[in] _N  成分インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NN 成分数
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_NIJK(_N,_I,_J,_K,_NN,_NI,_NJ,_NK,_VC) \
( (long long)(_NN) * _IDX_IJK(_I,_J,_K,_NI,_NJ,_NK,_VC) \
+ (long long)(_N) )



#endif /* _CIO_DEFINE_H_ */
