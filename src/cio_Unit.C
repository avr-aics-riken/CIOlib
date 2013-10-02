/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_Unit.C
 * @brief  cio_Unit Class
 * @author kero    
 */

#include "cio_DFI.h"
#include <unistd.h> // for gethostname() of FX10/K


/** cio_UnitElem class **/

// #################################################################
// コンストラクタ
cio_UnitElem::cio_UnitElem()
{

  Name="";
  Unit="";
  BaseName="";
  BaseValue=0.0;
  DiffName="";
  DiffValue=0.0;

}

// #################################################################
// コンストラクタ
cio_UnitElem::cio_UnitElem(const std::string _Name,
                           const std::string _Unit,
                           const std::string _BaseName,
                           const double      _BaseValue,
                           std::string       _DiffName,
                           double            _DiffValue)
{
  Name       = _Name;
  Unit       = _Unit;
  BaseName   = _BaseName;
  BaseValue  = _BaseValue;
  DiffName   = _DiffName;
  DiffValue  = _DiffValue;
}


// #################################################################
// デストラクタ
cio_UnitElem::~cio_UnitElem()
{


}

// #################################################################
// Unit要素の読込み
CIO::E_CIO_ERRORCODE 
cio_UnitElem::Read(cio_TextParser tpCntl,
                   const std::string _Name,
                   const std::string _BaseName,
                   const std::string _DiffName)
{

  std::string str,label;
  double dt;

  //単位系のの読込み
  Name = _Name;
  label="/Unit/"+Name;
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    return CIO::E_CIO_WARN_GETUNIT;
  }
  Unit=str;

  //値の読込み
  BaseName = _BaseName;
  label="/Unit/"+BaseName;
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    dt=0.0;
  }
  BaseValue=dt;

  //diffの読込み
  DiffName = _DiffName;
  if( !DiffName.empty() ) {
    label = "/Unit/"+DiffName;
    if ( !(tpCntl.GetValue(label, &dt )) )
    {
      dt=0.0;
    }
    DiffValue=dt;
  }

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// Unit要素の出力
CIO::E_CIO_ERRORCODE
cio_UnitElem::Write(FILE* fp, const unsigned tab)
{

  _CIO_WRITE_TAB(fp, tab);
  fprintf(fp, "%-11s = \"%s\"\n",Name.c_str(),Unit.c_str());
  _CIO_WRITE_TAB(fp, tab);
  fprintf(fp, "%-11s = %e\n",BaseName.c_str(),BaseValue);
  if( !DiffName.empty() ) {
    _CIO_WRITE_TAB(fp, tab);
    fprintf(fp, "%-11s = %e\n",DiffName.c_str(),DiffValue);
  }

  return CIO::E_CIO_SUCCESS;

}

/** cio_Uit class **/
// #################################################################
// コンストラクタ
cio_Unit::cio_Unit()
{

}

// #################################################################
// デストラクタ
cio_Unit::~cio_Unit()
{

  UnitList.clear();

}

// #################################################################
// Unitの読込み
CIO::E_CIO_ERRORCODE 
cio_Unit::Read(cio_TextParser tpCntl) 
{

  //Length
  cio_UnitElem unitLength;
  if( unitLength.Read(tpCntl, "Length","L0") == CIO::E_CIO_SUCCESS ) {
     UnitList.insert(map<std::string,cio_UnitElem>::value_type("Length",unitLength));
  }

  //Velocity
  cio_UnitElem unitVelocity;
  if( unitVelocity.Read(tpCntl,"Velocity","V0") == CIO::E_CIO_SUCCESS ) {
    UnitList.insert(map<std::string,cio_UnitElem>::value_type("Velocity",unitVelocity));
  }

  //Pressure
  cio_UnitElem unitPressure;
  if( unitPressure.Read(tpCntl,"Pressure","P0","DiffPrs") == CIO::E_CIO_SUCCESS ) {
    UnitList.insert(map<std::string,cio_UnitElem>::value_type("Pressure",unitPressure));
  }

  //Temperature
  cio_UnitElem unitTemperature;
  if( unitTemperature.Read(tpCntl,"Temperature","BaseTemp","DiffTemp") == CIO::E_CIO_SUCCESS ) {
    UnitList.insert(map<std::string,cio_UnitElem>::value_type("Temperature",unitTemperature));
  }

  return CIO::E_CIO_SUCCESS; 

}

// #################################################################
// 該当するUnitElemの取り出し
CIO::E_CIO_ERRORCODE 
cio_Unit::GetUnitElem(const std::string Name,
                      cio_UnitElem &unit)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合はNULLを返す
  if( it == UnitList.end() ) {
    return CIO::E_CIO_ERROR; 
  }

  //UnitElemを返す
  unit = (*it).second;

  return CIO::E_CIO_SUCCESS;

}

// #################################################################
// 単位の取り出し
std::string cio_Unit::GetUnit(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は空白を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return ""; 
  }

  //単位を返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.Unit;

}

// #################################################################
// ベース名の取り出し
std::string cio_Unit::GetBaseName(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は空白を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return ""; 
  }

  //ベース名を返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.BaseName;

}

// #################################################################
// ベース値の取り出し
double cio_Unit::GetBaseValue(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は0.0を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return 0.0;
  }

  //ベース値を返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.BaseValue;

}

// #################################################################
// Diff Name の取り出し
std::string cio_Unit::GetDiffName(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は空白を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return "";
  }

  //DiffNameを返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.DiffName;
}

// #################################################################
// Diff Valueの取り出し
double cio_Unit::GetDiffValue(const std::string Name, int &ret)
{
  map<std::string,cio_UnitElem>::iterator it;

  //Nameをキーにしてcio_UnitElemを検索
  it=UnitList.find(Name);

  //見つからなかった場合は0.0を返す
  if( it == UnitList.end() ) {
    ret = CIO::E_CIO_WARN_GETUNIT;
    return 0.0;
  }

  //Diff Valueを返す
  ret = CIO::E_CIO_SUCCESS;
  return (*it).second.DiffValue;

}

// #################################################################
// Unitの出力
CIO::E_CIO_ERRORCODE
cio_Unit::Write(FILE* fp, 
                const unsigned tab) 
{

  fprintf(fp, "Unit {\n");
  fprintf(fp, "\n");

  map<std::string,cio_UnitElem>::iterator it;
  for( it=UnitList.begin(); it!=UnitList.end(); it++ ) {
    if( (*it).second.Write(fp,tab+1) != CIO::E_CIO_SUCCESS ) return CIO::E_CIO_ERROR;
  }

  fprintf(fp, "\n");
  fprintf(fp, "}\n");
  fprintf(fp, "\n");

  return CIO::E_CIO_SUCCESS;

}

