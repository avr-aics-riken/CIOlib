/*
 *
 * fmr( File Rank Mapper)  
 *
 * Copyright (c) RIKEN AICS, Japan. All right reserved. 2013
 *
 *
*/

/**
 * @file   main.C
 * @brief  Staging_Utility main関数
 * @author kero
 */

#include "Staging_Utility.h"

//void print_dfi(Staging STG);
void print_info(Staging STG);
void set_arg(Staging &STG, int argc, char **argv);

int main( int argc, char **argv )
{

  //MPI Initialize
  if( MPI_Init(&argc,&argv) != MPI_SUCCESS )
  {
      std::cerr << "MPI_Init error." << std::endl;
      return false;
  }

  Staging STG;

  //引数の取り出し＆セット
  set_arg(STG,argc,argv);

  CIO::E_CIO_ERRORCODE ret;
  vector<int> readRankList; ///< 読込みランクリスト

  //読込みDFIファイルのループ
  for( int i=0; i<STG.m_dfi_fname.size(); i++ ) {

    //初期化、ファイルの読込み、DFIのインスタンス  
    if( !STG.Initial(STG.m_infofile, STG.m_dfi_fname[i]) ) {
      printf("ERROR Initial()\n");
      return 0;
    }

    //DFIのdirectory path get
    STG.m_inPath = CIO::cioPath_DirName(STG.m_dfi_fname[i]);

    bool isSameDiv = true;        ///< 分割数フラグ　true:1x1 false:MxN
    bool isSame = true;           ///< 粗密フラグ true:密 false:粗
    CIO::E_CIO_READTYPE readflag; ///< 読込み判定フラグ

    // 分割数フラグの設定
    for(int j=0; j<3; j++ ) {
      if( STG.m_Gdiv[j] != STG.dfi_Domain->GlobalDivision[j] ) 
      {
        isSameDiv = false;
      }
    }

    // 粗密フラグの設定
    if( STG.CheckGlobalVoxel(STG.m_GVoxel,(int *)STG.dfi_Domain->GlobalVoxel) == STG_E_GV_SAME ) 
    {
      isSame = true;
    } 
    else if( STG.CheckGlobalVoxel(STG.m_GVoxel,(int *)STG.dfi_Domain->GlobalVoxel) == STG_E_GVX2_SAME )
    {
      isSame = false;
    } else {
     printf("ERROR Dimension size : %d %d %d\n",STG.m_GVoxel[0],STG.m_GVoxel[1],STG.m_GVoxel[2]);
     return 0;
    }

    //読込み判定フラグの設定
    if( isSameDiv == true )
    {
      if( isSame == true ) 
      {
        readflag = CIO::E_CIO_SAMEDIV_SAMERES;
      } 
      else 
      {
        readflag = CIO::E_CIO_SAMEDIV_REFINEMENT;
      }
    } 
    else 
    {
      if( isSame == true )
      {
        readflag = CIO::E_CIO_DIFFDIV_SAMERES;
      }
      else
      {
        readflag = CIO::E_CIO_DIFFDIV_REFINEMENT;
      }
    }

    int numRank = STG.m_GRankInfo.size();

    //ランクマップの生成
    int* rankMap = STG.CreateRankMap();

    //STG.m_GRankInfoの生成
    STG.m_HeadTail=NULL;
    if( numRank>0 ) {
      if( numRank != STG.m_NumberOfRank ) {
        printf("ERROR MISMATCH NumberOfRank\n");
        return 0;
      }
      STG.m_HeadTail = new int[numRank][6];
      //ランク毎のXYZ方向のheadとtaileテーブルの生成
      if( !STG.CreateHeadTail(rankMap,STG.m_GRankInfo) ) {
        printf("ERROR CreateHeadTail()\n");
        return 0;
      }
    } else {
      size_t ndiv = STG.m_Gdiv[0]*STG.m_Gdiv[1]*STG.m_Gdiv[2];
      STG.m_HeadTail = new int[ndiv][6];
      //ランク毎のXYZ方向のheadとtaileテーブルの生成
      if( !STG.CreateHeadTail(rankMap) ) {
        printf("ERROR CreateHeadTail()\n");
        return 0;
      }
    }

    //並列数のセット
    numRank=STG.m_GRankInfo.size();

    //並列数のループ
    char tmp[20];
    if( STG.m_outPath == "" ) STG.m_outPath=".";
    int len = STG.m_outPath.size()+7;

    for(int j=0; j<numRank; j++) {
       //読込みランクリストの生成
       readRankList.clear();

       cio_Domain domain;
       for(int k=0; k<3; k++ ) {
         domain.GlobalOrigin[k] = STG.dfi_Domain->GlobalOrigin[k];
         domain.GlobalRegion[k] = STG.dfi_Domain->GlobalRegion[k];
         domain.GlobalVoxel[k]  = STG.dfi_Domain->GlobalVoxel[k];
         domain.GlobalDivision[k]  = STG.dfi_Domain->GlobalDivision[k];
       }
       domain.ActiveSubdomainFile = STG.dfi_Domain->ActiveSubdomainFile;

       ret=STG.dfi_Process->CheckReadRank(domain,
           (const int *)STG.m_GRankInfo[j].HeadIndex,
           (const int *)STG.m_GRankInfo[j].TailIndex,readflag,readRankList);

       if( ret != CIO::E_CIO_SUCCESS ) return 0;

       //ファイルのコピー
       STG.FileCopy(readRankList,STG.m_GRankInfo[j].RankID);

    }

    //dfiファイル出力
    if( !STG.OutputDFI(STG.m_dfi_fname[i], rankMap) ) {
      printf("EEROR OutputDFI()\n");
      return 0;
    }
  }

  printf("####### normal end ######\n");

  return 1;
}

void set_arg(Staging &STG, int argc, char **argv)
{
  STG.m_outPath="";
  char *p;
  for(int i=1; i<argc; i++ ) {
    p=argv[i];
    if( *p=='-' ) {
      i++;
      p++;
      switch(*p) {
        case 'i' :
         STG.m_infofile = argv[i++];
         break;
        case 's' :
         STG.m_step = atoi(argv[i++]);
         break;
        case 'o' :
         STG.m_outPath = argv[i++];
         break;
      }
      i--;
    } else {
      STG.m_dfi_fname.push_back(argv[i]);
    }
  }

  printf("input proc file : %s\n",STG.m_infofile.c_str());
  printf("m_step          : %d\n",STG.m_step);
  printf("m_outPath       : %s\n",STG.m_outPath.c_str());
  printf("dfi_fname size  : %d\n",(int)STG.m_dfi_fname.size());
  for(int i=0;i<STG.m_dfi_fname.size();i++ ) printf("dfi_fname       : %s\n",STG.m_dfi_fname[i].c_str());

}
