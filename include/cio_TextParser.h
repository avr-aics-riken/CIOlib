#ifndef _CIO_TEXTPARSER_H_
#define _CIO_TEXTPARSER_H_

/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/**
 * @file   cio_TextParser.h 
 * @brief  TextParser Control class Header
 * @author kero
 */

#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "string.h"

//#include "cpm_ParaManager.h"

#include "cio_Define.h"
#include "TextParser.h"

using namespace std;

/*
class MediumTableInfo {
public:
  int type;
  string label;
  map<string, REAL_TYPE> m_fval;

public:
  MediumTableInfo() 
  {
    type = -1;
  }
  
  ~MediumTableInfo() {}
  
};
*/

class cio_TextParser {

private:
	TextParser* tp;  ///< テキストパーサ

public:
  /** コンストラクタ */
	cio_TextParser(){};

  /**　デストラクタ */
	~cio_TextParser(){};
	

  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（整数型）
   * @param [int] label 取得するベクトルのラベル（絶対パス）
   * @param [out] vec   ベクトル格納配列ポインタ
   * @param [in]  nvec  ベクトルサイズ
   */
	bool GetVector(const string label, int *vec, const int nvec);
  
  
  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（実数型）
   * @param [in]  label  取得するベクトルのラベル（絶対パス）
   * @param [out] vec    ベクトル格納配列ポインタ
   * @param [in]  nvec   ベクトルサイズ
   */
	//bool GetVector(const string label, REAL_TYPE *vec, const int nvec);
	bool GetVector(const string label, double *vec, const int nvec);

  
  /**
   * @brief TextParser入力ファイルからベクトル値を取得する（文字列型）
   * @param [in]  label  取得するベクトルのラベル（絶対パス）
   * @param [out] vec    ベクトル格納配列ポインタ
   * @param [in]  nvec   ベクトルサイズ
   */
	bool GetVector(const string label, string *vec, const int nvec);
  
  
  /**
   * @brief TextParser入力ファイルから変数を取得する（整数型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, int *ct);
  
  
  /**
   * @brief TextParser入力ファイルから変数を取得する（実数型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [out] ct    変数格納ポインタ
   */
	//bool GetValue(const string label, REAL_TYPE *ct);
	bool GetValue(const string label, double *ct);
  
  
  /**
   * @brief TextParser入力ファイルから変数を取得する（文字列型）
   * @param [in]  label 取得する変数のラベル（絶対パス）
   * @param [out] ct    変数格納ポインタ
   */
	bool GetValue(const string label, string *ct);

	
  /**
   * @brief ラベルの有無をチェック
   * @param [in] label チェックするラベル（絶対パス）
   */
	bool chkLabel(const string label);
  
  
  /**
   * @brief ノードの有無をチェック
   * @param [in] label チェックするノード（絶対パス）
   */
	bool chkNode(const string label);
  
  
  /**
   * @brief ノード以下のnnode番目の文字列を取得する
   * @param [in]  label ノードの絶対パス
   * @param [in]  nnode 取得する文字列が現れる順番
   * @param [out] ct    取得した文字列
   */
	bool GetNodeStr(const string label, const int nnode, string *ct);
  
  
  /**
   * @brief ノード以下のラベルの数を数える
   * @param [in] label ラベルを数えるノードの絶対パス
   * @retval ラベルの数（エラー、もしくはない場合は-1を返す）
   */
  int countLabels(const string label);
  
  /**
   * @brief TextParserLibraryのインスタンス生成
   * @retrun エラーコード
   */
  void getTPinstance();
  
  
  /**
   * @brief TextParserオブジェクトに入力ファイルをセットする
   * @param [in] filename 入力ファイル名
   * @retval エラーコード
   */
	int readTPfile(const string filename);
  
  
  /** テキストパーサーの内容を破棄 */
  int remove()
  {
    return tp->remove();
  }

};

#endif // _CIO_TXETPARSER_H_
