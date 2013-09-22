#ifndef _CIO_UNIT_H_
#define _CIO_UNIT_H_

/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_Unit.h
 * @brief  cio_UnitElem & cio_Unit Class Header
 * @author kero    
 */

class cio_UnitElem {

public:

  std::string Name;         ///<単位の種類名(Length,Velovity,,,)
  std::string Unit;         ///<単位のラベル(m,m/s,Pa,,,,)
  std::string BaseName;     ///<代表値のラベル(L0,P0,,,,)
  double      BaseValue;    ///<代表値
  std::string DiffName;     ///<差があるときの名前、空の時DiffValue未定義
  double      DiffValue;    ///<差の値

  /** コンストラクタ */
  cio_UnitElem();

  /** コンストラクタ */
  cio_UnitElem(const std::string _Name,
               const std::string _Unit,
               const std::string _BaseName,
               const double      _BaseValue,
               std::string       _DiffName = "",
               double            _DiffiValue = 0.0);

  /** デストラクタ */
  ~cio_UnitElem();

  /**
   * @brief Unit要素の読込み
   * @param [in]  tpCntl      cio_TextParserクラス 
   * @param [in]  _Name       単位種類("Length","Velocity","Pressure",,,,)
   * @param [in]  _BaseName   名前("L0","V0","P0",,,,)
   * @param [in]  _DiffName   差があるときの差の名前("DiffPrrs","DiffTemp",,,,) 
   * @return error code
   */
   CIO::E_CIO_ERRORCODE 
   Read(cio_TextParser tpCntl,
        const std::string _Name,
        const std::string _BaseName,
        const std::string _DiffName=""); 

  /**
   * @brief DFIファイル:Unit要素を出力する
   * @param [in]  fp      ファイルポインタ
   * @param [in]  tab     インデント
   * @return error code
   */
   CIO::E_CIO_ERRORCODE
   Write(FILE* fp, const unsigned tab);

};


/** index.dfi ファイルの Unit */
class cio_Unit { 

public:

  map<std::string,cio_UnitElem> UnitList;

  /** コンストラクタ **/
  cio_Unit();

  /** デストラクタ **/
  ~cio_Unit();

  /**
   * @brief read Unit(inde.dfi)
   * @param [in]   tpCntl  cio_TextParserクラス 
   * @return error code
   */
  CIO::E_CIO_ERRORCODE 
  Read(cio_TextParser tpCntl);

  /**
   * @brief 該当するUnitElemの取り出し
   * @param [in]  Name 取り出す単位の種類
   * @param [out] unit 取得したunitクラス
   * @return error code
   */ 
  CIO::E_CIO_ERRORCODE 
  GetUnitElem(const std::string Name,
              cio_UnitElem &unit);

  /**
   * @brief 単位の取り出し("m","cm",,,,,)
   * @param [in]  Name 取り出す単位の種類 
   * @param [out] ret  return code
   * @return 単位
   */
  std::string GetUnit(const std::string Name, int &ret);

  /**
   * @brief ベース名値の取り出し("L0","P0",,,,)
   * @param [in]  Name 取り出す単位の種類
   * @param [out] ret  return code
   * @return ベース名 
   */
  std::string GetBaseName(const std::string Name, int &ret);

  /**
   * @brief ベース値の取り出し
   * @param [in]  Name 取り出す単位の種類
   * @param [out] ret  return code
   * @return ベース値
   */
  double GetBaseValue(const std::string Name, int &ret);

  /**
   * @brief DiffNameの取り出し
   * @param [in]  Name 取り出す単位の種類
   * @param [out] ret  return code
   * @return true : DiffName
   */
  std::string GetDiffName(const std::string Name, int &ret);

  /**
   * @brief Diff Valueの取り出し
   * @param [in]  Name 取り出す単位の種類
   * @param [out] ret  return code
   * @return true : Diff Value
   */
  double GetDiffValue(const std::string Name, int &ret);
 
  /**
   * @brief DFIファイル:Unit要素を出力する
   * @param [in]  fp      ファイルポインタ
   * @param [in]  tab     インデント
   * @return error code
   */
   CIO::E_CIO_ERRORCODE
   Write(FILE* fp, const unsigned tab);

};

#endif // _CIO_UNIT_H_
