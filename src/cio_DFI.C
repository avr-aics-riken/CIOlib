/*
 * CIO - Cartesian Input / Output library
 *
 * Copyright (c) RIKEN AICS, Japan. All right reserved. 2013
 *
 */

/** 
 * @file   cio_DFI.C
 * @brief  cio_DFI Class
 * @author kero    
 */

#include <unistd.h> // for gethostname() of FX10/K

#include "cio_DFI.h"
#include "cio_DFI_SPH.h"
#include "cio_DFI_BOV.h"

// #################################################################
// コンストラクタ
cio_DFI::cio_DFI()
{

 m_dfi_mng = 0;
 m_outSlice = false;
 m_start_type = 0;
 m_RankID = 0;

}


// #################################################################
// デストラクタ
cio_DFI::~cio_DFI()
{

}

// #################################################################
//
cio_DFI* cio_DFI::ReadInit(MPI_Comm comm, string DfiName)
{

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  cio_TextParser tpCntl;

  //index.dfi read
  //TPインスタンス
  tpCntl.getTPinstance();

  FILE*fp = NULL;
  if( !(fp=fopen(DfiName.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",DfiName.c_str());
    return false;
  }
  fclose(fp);

  //入力ファイル index.dfi をセット
  int ierror = 0;
  ierror = tpCntl.readTPfile(DfiName);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",DfiName.c_str());
    return NULL;
  }

  cio_FileInfo F_info;
  if( readFileInfo(DfiName,tpCntl,F_info) != CIO_SUCCESS ) 
  {
    printf("\tFileInfo Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  cio_FilePath F_path;
  if( readFilePath(DfiName,tpCntl,F_path) != CIO_SUCCESS )
  {
    printf("\tFilePath Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  cio_Unit unit;
  if( readUnit(DfiName,tpCntl,unit) != CIO_SUCCESS )
  {
    printf("\tUnit Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  vector<cio_Slice> TimeSlice;
  cio_Slice slice;
  if( readSlice(DfiName,tpCntl,TimeSlice,slice) != CIO_SUCCESS )
  {
    printf("\tTimeSlice Data Read error %s\n",DfiName.c_str());
    return NULL;
  }

  //TextParserの破棄
  tpCntl.remove();

  //proc file name set
  string procfile = F_path.Process;

  //proc.dfi read
  //TPインスタンス
  tpCntl.getTPinstance();

  fp = NULL;
  if( !(fp=fopen(procfile.c_str(),"rb")) ) {
    printf("Can't open file. (%s)\n",procfile.c_str());
    return false;
  }
  fclose(fp);

  //入力ファイル proc.dfi をセット
  ierror = tpCntl.readTPfile(procfile);
  if ( ierror )
  {
    printf("\tinput file not found '%s'\n",procfile.c_str());
    return NULL;
  }

  cio_Domain domain;
  if( readDomain(procfile,tpCntl,domain) != CIO_SUCCESS ) 
  {
    printf("\tDomain Data Read error %s\n",procfile.c_str());
    return NULL;
  }


  cio_MPI mpi;
  if( readMPI(procfile,tpCntl,domain,mpi) != CIO_SUCCESS )
  {
    printf("\tMPI Data Read error %s\n",procfile.c_str());
    return NULL;
  }

  vector<cio_Rank> RankInfo;
  cio_Rank rank;
  if( readRank(procfile,tpCntl,RankInfo,rank) != CIO_SUCCESS )
  {
    printf("\tProcess Data Read error %s\n",procfile.c_str());
    return NULL;
  }

  //TextParserの破棄
  tpCntl.remove();

  //std::string fmt = F_info.FileFormat;

  cio_DFI *dfi = NULL;
  if( F_info.FileFormat == E_CIO_FMT_SPH ) {
    dfi = new cio_DFI_SPH(F_info, F_path, unit, domain, mpi, TimeSlice, RankInfo);
  } else if( F_info.FileFormat == E_CIO_FMT_BOV ) {
    dfi = new cio_DFI_BOV(F_info, F_path, unit, domain, mpi, TimeSlice, RankInfo);
  } else {
    dfi = NULL;
    return dfi;
  }

  dfi->m_comm = comm;
  dfi->m_indexDfiName = DfiName;
  dfi->m_RankID = RankID;

  if( !strcasecmp(dfi->DFI_Finfo.TimeSliceDir.c_str() , "on" ) ) dfi->m_outSlice=true;
  else dfi->m_outSlice=false;

  //printf("TimeSliceDirectory : %s\n",dfi->DFI_Finfo.TimeSliceDir.c_str());

  int div[3];
  for(int i=0;i<3;i++) div[i]=dfi->DFI_Domain.GlobalDivision[i];
  int nRank = dfi->RankInfo.size();

  //ActiveSubdomainファイルの読込み処理
  if( !domain.ActiveSubdomain.empty() && nRank<=0 ) {
    int divSudomain[3] = {0,0,0};
    cio_ErrorCode ret = ReadActiveSubdomainFile( domain.ActiveSubdomain, 
                        dfi->DFI_subDomainInfo, divSudomain);
    if( ret != CIO_SUCCESS ) return NULL;
  } else {
    if( CheckData(nRank,div,dfi->DFI_subDomainInfo) != CIO_SUCCESS ) 
    return NULL;
  }

  dfi->DFI_HeadTail = NULL;
  if( nRank>0 ) {
    dfi->DFI_HeadTail = new int[nRank][6];
  } else {
    int *rankMap = CreateRankMap(div,dfi->DFI_subDomainInfo);
    size_t ndiv = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
    dfi->DFI_HeadTail = new int[ndiv][6];
    CreateHeadTail(div, rankMap, dfi->DFI_Domain.GlobalVoxel, dfi->DFI_HeadTail,
                   dfi->RankInfo,rank); 
  }

  SetHeadTail(dfi->RankInfo, dfi->DFI_HeadTail, 
              dfi->DFI_headX, dfi->DFI_headY, dfi->DFI_headZ);

  CreateHeadMap(dfi->DFI_headX, dfi->DFI_mapX);
  CreateHeadMap(dfi->DFI_headY, dfi->DFI_mapY);
  CreateHeadMap(dfi->DFI_headZ, dfi->DFI_mapZ);

  dfi->DFI_rankMap = CreateActiveRankMap(dfi->RankInfo,dfi->DFI_HeadTail,
                     dfi->DFI_mapX,dfi->DFI_mapY,dfi->DFI_mapZ);

#ifdef _CIO_DEBUG_
  if( RankID==0 ) {

    printf("FileFormat : %d\n",dfi->DFI_Finfo.FileFormat);

    int ndiv=dfi->DFI_mapX.size()*dfi->DFI_mapY.size()*dfi->DFI_mapZ.size();
    for(int i=0;i<ndiv;i++) printf("rankMap[%d] : %d\n",i,dfi->DFI_rankMap[i]);
  }
#endif

  return dfi;

}

// #################################################################
// 
cio_DFI* cio_DFI::get_dfi(cio_FileInfo out_F_info, cio_FilePath out_F_path, cio_Unit out_unit,
                          cio_Domain out_domain, cio_MPI out_mpi, 
                          vector<cio_Slice> out_TSlice, vector<cio_Rank> out_RankInfo)
{
  if( out_F_info.FileFormat == E_CIO_FMT_SPH ) {
    return  new cio_DFI_SPH(out_F_info, out_F_path, out_unit, out_domain, out_mpi, 
                          out_TSlice, out_RankInfo);
  } else if( out_F_info.FileFormat == E_CIO_FMT_BOV ) {

    return  new cio_DFI_BOV(out_F_info, out_F_path, out_unit, out_domain, out_mpi,
                          out_TSlice, out_RankInfo);
  }
  return NULL;
}
// #################################################################
// 初期化
void cio_DFI::InitDFI()
{
}
//

// #################################################################
//
int cio_DFI::readFileInfo(string dfifile, cio_TextParser tpCntl, cio_FileInfo &finfo)
{

  string str;
  string label,label_base,label_leaf;
  int ct;
  //REAL_TYPE dt;
  double dt;

  //Directorypath
  label = "/FileInfo/DirectoryPath";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.DirectoryPath=str;

  //TimeSilceDirectory
  label = "/FileInfo/TimeSliceDirectory";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.TimeSliceDir=str;

  //Prefix
  label = "/FileInfo/Prefix";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.Prefix=str;

  //FileFormat
  label = "/FileInfo/FileFormat";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  if( !strcasecmp(str.c_str(),"sph" ) ) {
    finfo.FileFormat=E_CIO_FMT_SPH;
  }
  else if( !strcasecmp(str.c_str(),"bov" ) ) {
    finfo.FileFormat=E_CIO_FMT_BOV;
  }
  else finfo.FileFormat=E_CIO_FMT_UNKNOWN;

  //GuidCell
  label = "/FileInfo/GuideCell";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.GuideCell=ct;

  //DataType
  label = "/FileInfo/DataType";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.DataType=str;

  //Endian
  label = "/FileInfo/Endian";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.Endian=str;

  //ArrayShape  
  label = "/FileInfo/ArrayShape";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.ArrayShape=str;

  //Componet  
  label = "/FileInfo/Component";
  if ( !(tpCntl.GetValue(label, &ct )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  finfo.Component=ct;

  return CIO_SUCCESS;
}

// #################################################################
//
int cio_DFI::readFilePath(string dfifile, cio_TextParser tpCntl, cio_FilePath &fpath)
{

  string str;
  string label;

  //Process
  label = "/FilePath/Process";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  fpath.Process=str;

  return CIO_SUCCESS;

}
// #################################################################
//
int cio_DFI::readUnit(string dfifile, cio_TextParser tpCntl, cio_Unit &unit)
{

  string str;
  string label;
  int ct;
  //REAL_TYPE dt;
  double dt;

  //Length
  label = "/Unit/Length";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
/*
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
*/
    str="NonDimensional"; 
  }
  unit.Length=str;

  //L0
  label = "/Unit/L0";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
/*
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
*/
    dt=1.0;
  }
  unit.L0=dt;

  //Velocity
  label = "/Unit/Velocity";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
/*
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
*/
    str="NonDimensional";
  }
  unit.Velocity=str;

  //L0
  label = "/Unit/V0";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
/*
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
*/
    dt=1.0; 
  }
  unit.V0=dt;

  //Pressure
  label = "/Unit/Pressure";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
/*
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
*/
    str="NonDimensional"; 
  }
  unit.Pressure=str;

  //P0
  label = "/Unit/P0";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
/*
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR; 
*/
    dt=0.0;
  }
  unit.P0=dt;

  //DiffPrs
  label = "/Unit/DiffPrs";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
/*
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
*/
    dt=0.0; 
  }
  unit.DiffPrs=dt;

  //Temperatur
  label = "/Unit/Temperatur";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR; 
    unit.Temperatur="";
  } else {
    unit.Temperatur=str;
  }

  //BaseTemp
  label = "/Unit/BaseTemp";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR; 
    unit.BaseTemp=0.0;
  } else {
    unit.BaseTemp=dt;
  }

  //DiffTemp
  label = "/Unit/DiffTemp";
  if ( !(tpCntl.GetValue(label, &dt )) )
  {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR; 
    unit.DiffTemp = 0.0;
  } else {
    unit.DiffTemp=dt;
  }

  return CIO_SUCCESS; 

}

// #################################################################
//
int cio_DFI::readSlice(string dfifile, cio_TextParser tpCntl, vector<cio_Slice> &TimeSlice, cio_Slice  slice)
{

  string str;
  string label,label_base,label_leaf,label_leaf_leaf;
  int ct;
  //REAL_TYPE dt;
  double dt;
  int nnode=0;

  //TimeSlice
  nnode=0;
  label_base = "/TimeSlice";
  if ( tpCntl.chkNode(label_base) )  //があれば
  {
    nnode = tpCntl.countLabels(label_base);
  }


  for (int i=0; i<nnode; i++) {

    int ncnt=0;

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      return CIO_ERROR;
    }
    if( strcasecmp(str.substr(0,5).c_str(), "Slice") ) continue;
    label_leaf=label_base+"/"+str;

    //Step
    label = label_leaf + "/Step";
    if ( !(tpCntl.GetValue(label, &ct )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    else {
      slice.step=ct;
    }

    ncnt++;

    //Time
    label = label_leaf + "/Time";
    if ( !(tpCntl.GetValue(label, &dt )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    else {
      slice.time= dt;
    }

    ncnt++;

    //AveragedStep
    label = label_leaf + "/AveragedStep";
    if ( !(tpCntl.GetValue(label, &ct )) ) {
      slice.AveragedStep=-1;
    }
    else {
      slice.AveragedStep= ct;
      ncnt++;
    }


    //AveragedTime
    label = label_leaf + "/AveragedTime";
    if ( !(tpCntl.GetValue(label, &dt )) ) {
      slice.AveragedTime=0.0;
    }
    else {
      slice.AveragedTime= dt;
      ncnt++;
    }


    //MinMax
    int ncomp=0;
    label_leaf_leaf = label_leaf + "/MinMax";
    if ( tpCntl.chkNode(label_leaf_leaf) )  //があれば
    {
      ncomp = tpCntl.countLabels(label_leaf_leaf);
    }
   
    ncnt++;
 
    slice.Min.clear();
    slice.Max.clear();

    for ( int j=0; j<ncomp; j++ ) {

      if(!tpCntl.GetNodeStr(label_leaf,j+ncnt,&str))
      {
        printf("\tParsing error : No Elem name\n");
        return CIO_ERROR;
      } 
      if( strcasecmp(str.substr(0,6).c_str(), "minmax") ) continue;
      label_leaf_leaf = label_leaf+"/"+str;

      label = label_leaf_leaf + "/Min";
      if ( !(tpCntl.GetValue(label, &dt )) ) {
        printf("\tParsing error : fail to get '%s'\n",label.c_str());
        return CIO_ERROR;
      }
      else {
        slice.Min.push_back(dt);
      }

      label = label_leaf_leaf + "/Max";
      if ( !(tpCntl.GetValue(label, &dt )) ) {
        printf("\tParsing error : fail to get '%s'\n",label.c_str());
        return CIO_ERROR;
      }
      else {
        slice.Max.push_back(dt);
      }

    }

   TimeSlice.push_back(slice); 

  }

  return CIO_SUCCESS;

}

// #################################################################
//
int cio_DFI::readDomain(string dfifile, cio_TextParser tpCntl, cio_Domain &domain)
{

  string str;
  string label;
  int ct;
  //REAL_TYPE dt;
  //REAL_TYPE v[3];
  double dt;
  double v[3];

  //GlobalOrign
  label = "/Domain/GlobalOrigin";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalOrigin[0]=v[0];
  domain.GlobalOrigin[1]=v[1];
  domain.GlobalOrigin[2]=v[2];

  //GlobalRegion
  label = "/Domain/GlobalRegion";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalRegion[0]=v[0];
  domain.GlobalRegion[1]=v[1];
  domain.GlobalRegion[2]=v[2];

  //Global_Voxel
  label = "/Domain/GlobalVoxel";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalVoxel[0]=v[0];
  domain.GlobalVoxel[1]=v[1];
  domain.GlobalVoxel[2]=v[2];

  //Global_Division
  label = "/Domain/GlobalDivision";
  for (int n=0; n<3; n++) v[n]=0.0;
  if ( !(tpCntl.GetVector(label, v, 3 )) ) 
  {
    printf("\tParsing error : fail to get '%s'\n",label.c_str());
    return CIO_ERROR;
  }
  domain.GlobalDivision[0]=v[0];
  domain.GlobalDivision[1]=v[1];
  domain.GlobalDivision[2]=v[2];

  //ActiveSubdomain
  label = "/Domain/ActiveSubdomain";
  if ( !(tpCntl.GetValue(label, &str )) )
  {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR;
    str="";
  }
  domain.ActiveSubdomain=str;

  return CIO_SUCCESS;

}

// #################################################################
//
int cio_DFI::readMPI(string dfifile, cio_TextParser tpCntl, cio_Domain domain, cio_MPI &mpi)
{

  string str;
  string label;
  int ct;

  //NumberOfRank
  label = "/MPI/NumberOfRank";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //printf("\tParsing warning : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR;
    ct = domain.GlobalDivision[0]*domain.GlobalDivision[1]*domain.GlobalDivision[2];
  }
  else {
    mpi.NumberOfRank = ct;
  }

  //NumberOfGroup
  label = "/MPI/NumberOfGroup";
  if ( !(tpCntl.GetValue(label, &ct )) ) {
    //printf("\tParsing error : fail to get '%s'\n",label.c_str());
    //printf("\tParsing warning : fail to get '%s'\n",label.c_str());
    //return CIO_ERROR;
    ct = 1;
  }
  else {
    mpi.NumberOfGroup = ct;
  }

  return CIO_SUCCESS;

}

// #################################################################
//
int cio_DFI::readRank(string dfifile, cio_TextParser tpCntl, vector<cio_Rank> &RankInfo, cio_Rank rank)
{

  string str;
  string label,label_base,label_leaf;
  int ct;
  //REAL_TYPE dt;
  //REAL_TYPE v[3];
  double dt;
  double v[3];
  int nnode=0;

  //Process 
  nnode=0;
  label_base = "/Process";
  if ( tpCntl.chkNode(label_base) )  //nodeがあれば
  {
    nnode = tpCntl.countLabels(label_base);
  }

  for (int i=0; i<nnode; i++) {

    if(!tpCntl.GetNodeStr(label_base,i+1,&str))
    {
      printf("\tParsing error : No Elem name\n");
      return CIO_ERROR;
    }
    if( strcasecmp(str.substr(0,4).c_str(), "Rank") ) continue;
    label_leaf=label_base+"/"+str;

    //ID
    label = label_leaf + "/ID";
    if ( !(tpCntl.GetValue(label, &ct )) ) {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    else {
      rank.RankID= ct;
    }

    //HostName
    label = label_leaf + "/HostName";
    if ( !(tpCntl.GetValue(label, &str )) ) {
      printf("\tParsing error : fail to get '%s'\n", label.c_str());
      return CIO_ERROR;
    }
    rank.HostName= str;

    //VoxelSize
    label = label_leaf + "/VoxelSize";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
    {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    rank.VoxelSize[0]=v[0];
    rank.VoxelSize[1]=v[1];
    rank.VoxelSize[2]=v[2];

    //HeadIndex
    label = label_leaf + "/HeadIndex";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
    {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    rank.HeadIndex[0]=v[0];
    rank.HeadIndex[1]=v[1];
    rank.HeadIndex[2]=v[2];

    //TailIndex
    label = label_leaf + "/TailIndex";
    for (int n=0; n<3; n++) v[n]=0.0;
    if ( !(tpCntl.GetVector(label, v, 3 )) ) 
    {
      printf("\tParsing error : fail to get '%s'\n",label.c_str());
      return CIO_ERROR;
    }
    rank.TailIndex[0]=v[0];
    rank.TailIndex[1]=v[1];
    rank.TailIndex[2]=v[2];

    RankInfo.push_back(rank); 

  }

  return CIO_SUCCESS;

}

// #################################################################
// 配列形状の取り出し
std::string cio_DFI::getArrayShape()
{
  return DFI_Finfo.ArrayShape;
}

// #################################################################
// データタイプの取り出し
std::string cio_DFI::getDataType()
{
  return DFI_Finfo.DataType;
}

// #################################################################
// 成分数の取り出し
int cio_DFI::getComponent()
{
  return DFI_Finfo.Component;
}

// #################################################################
// データタイプの取り出し
cio_DFI::E_CIO_DTYPE cio_DFI::get_cio_Datatype(string datatype)
{

  if     ( !strcasecmp(datatype.c_str(),"Int8"   ) ) return E_CIO_INT8;
  else if( !strcasecmp(datatype.c_str(),"Int16"  ) ) return E_CIO_INT16;
  else if( !strcasecmp(datatype.c_str(),"Int32"  ) ) return E_CIO_INT32;
  else if( !strcasecmp(datatype.c_str(),"Int64"  ) ) return E_CIO_INT64;
  else if( !strcasecmp(datatype.c_str(),"UInt8"  ) ) return E_CIO_UINT8;
  else if( !strcasecmp(datatype.c_str(),"UInt16" ) ) return E_CIO_UINT16;
  else if( !strcasecmp(datatype.c_str(),"UInt32" ) ) return E_CIO_UINT32;
  else if( !strcasecmp(datatype.c_str(),"UInt64" ) ) return E_CIO_UINT64;
  else if( !strcasecmp(datatype.c_str(),"Float32") ) return E_CIO_FLOAT32;
  else if( !strcasecmp(datatype.c_str(),"Float64") ) return E_CIO_FLOAT64;

  return E_CIO_DUMMY;
}

int cio_DFI::get_cio_Datasize(E_CIO_DTYPE Dtype)
{

  if     ( Dtype == E_CIO_INT8    ) return sizeof(char);
  else if( Dtype == E_CIO_INT16   ) return sizeof(short);
  else if( Dtype == E_CIO_INT32   ) return sizeof(int);
  else if( Dtype == E_CIO_INT64   ) return sizeof(long);
  else if( Dtype == E_CIO_UINT8   ) return sizeof(unsigned char);
  else if( Dtype == E_CIO_UINT16  ) return sizeof(unsigned short);
  else if( Dtype == E_CIO_UINT32  ) return sizeof(unsigned int);
  else if( Dtype == E_CIO_UINT64  ) return sizeof(unsigned long long);
  else if( Dtype == E_CIO_FLOAT32 ) return sizeof(float);
  else if( Dtype == E_CIO_FLOAT64 ) return sizeof(double);
  else return 0;

}
// #################################################################
// Create Domain & Process  
void cio_DFI::cio_Create_Domain(MPI_Comm comm, 
                                int G_voxel[3], int G_division[3],
                                int head[3], int tail[3],
                                cio_Domain &G_domain, vector<cio_Rank> &G_RankInfo,
                                cio_Rank G_Rank)
{

  for(int i=0; i<3; i++ ) {
    G_domain.GlobalVoxel[i]   =G_voxel[i];
    G_domain.GlobalDivision[i]=G_division[i];
  }

  int RankID;
  MPI_Comm_rank( comm, &RankID );

  int nrank;
  MPI_Comm_size( comm, &nrank );

  int *head_All  = new int[3*nrank];
  int *tail_All  = new int[3*nrank];
  for(int i=0; i<3; i++) {
     head_All[RankID*3+i] =head[i];
     tail_All[RankID*3+i] =tail[i];
  }

  if( nrank > 1 ) {
     MPI_Allgather(head,3,MPI_INT,head_All,3,MPI_INT,comm);
     MPI_Allgather(tail,3,MPI_INT,tail_All,3,MPI_INT,comm);
  }

  for(int i=0; i<nrank; i++) {
     G_Rank.RankID=i;
     for(int j=0; j<3; j++) G_Rank.HeadIndex[j]=head_All[i*3+j];
     for(int j=0; j<3; j++) G_Rank.TailIndex[j]=tail_All[i*3+j];
     for(int j=0; j<3; j++) G_Rank.VoxelSize[j]=G_Rank.TailIndex[j]-G_Rank.HeadIndex[j]+1;
     G_RankInfo.push_back(G_Rank);
  }

  if( head_All ) delete [] head_All;
  if( tail_All ) delete [] tail_All;

}

// #################################################################
// 粗密判定
cio_EGlobalVoxel cio_DFI::CheckGlobalVoxel(int Gvoxel[3], int DFI_Gvoxel[3])
{

  if( Gvoxel[0] == DFI_Gvoxel[0]   &&
      Gvoxel[1] == DFI_Gvoxel[1]   &&
      Gvoxel[2] == DFI_Gvoxel[2]   ) return CIO_E_GV_SAME;

  if( Gvoxel[0] == DFI_Gvoxel[0]*2 &&
      Gvoxel[1] == DFI_Gvoxel[1]*2 &&
      Gvoxel[2] == DFI_Gvoxel[2]*2 ) return CIO_E_GVX2_SAME;
  
  return CIO_E_OTHER;
}

// #################################################################
// 粗密ファイルのMxN判定
bool cio_DFI::CheckMxN(vector<int> &rankList, int head[3], int tail[3], int gc, int dfi_gc)
{

  int dfi_head[3],dfi_tail[3];

  for(int i=0; i<rankList.size(); i++) {
    int n=rankList[i];
    for(int j=0; j<3; j++) {
      dfi_head[j]=RankInfo[n].HeadIndex[j]*2-1;
      dfi_tail[j]=RankInfo[n].TailIndex[j]*2;
    }

  if( head[0] != dfi_head[0] || head[1] != dfi_head[1] || head[2] != dfi_head[2] ||
      tail[0] != dfi_tail[0] || tail[1] != dfi_tail[1] || tail[2] != dfi_tail[2] ||
      gc != dfi_gc ) return false;
  }
  return true;
}

// #################################################################
// 読込み範囲を求める
bool cio_DFI::CheckReadArea(int head[3], int tail[3], int gc, int DFI_head[3],
           int DFI_tail[3], int DFI_gc, cio_EGlobalVoxel readflag, int sta[3], int end[3])
{

  int cal_gc=0;
  if( DFI_gc > 0 || gc > 0 ) {
    if( gc <= DFI_gc ) {
      cal_gc = gc;
    } else {
      cal_gc = DFI_gc;
    }
  } 

  for(int i=0; i<3; i++) {
   sta[i]=0;
   end[i]=0;
  }

  for( int i=0; i<3; i++ ) {
    sta[i] = max(head[i],DFI_head[i]) -1;
    end[i] = min(tail[i],DFI_tail[i]) -1;
    if( sta[i] == 0 ) sta[i] -= cal_gc;
    else if( head[i]>DFI_head[i] ) sta[i] -= gc;

    if( readflag == CIO_E_GV_SAME ) {
      if( (end[i]+1) == DFI_Domain.GlobalVoxel[i] ) end[i] += cal_gc;
      else if( tail[i]<DFI_tail[i] ) end[i] += gc;
    } else {
      if( (end[i]+1) == DFI_Domain.GlobalVoxel[i]*2 ) end[i] += cal_gc;
      else if( tail[i]<DFI_tail[i] ) end[i] += gc;
    }
  }

  if( head[0] == DFI_head[0] && head[1] == DFI_head[1] && head[2] == DFI_head[2] &&
      tail[0] == DFI_tail[0] && tail[1] == DFI_tail[1] && tail[2] == DFI_tail[2] &&
      gc == DFI_gc ) return true;

  return false;

}

// #################################################################
// 読込みランクファイルリストの作成 
void cio_DFI::CreateRankList(int head[3], int tail[3], int gc, cio_EGlobalVoxel readflag,
                             vector<int> &rankList) 
{
  int StartEnd[6];
  int ndiv = DFI_Domain.GlobalDivision[0]*DFI_Domain.GlobalDivision[1]*DFI_Domain.GlobalDivision[2];

  int (*HeadTail)[6] = new int[ndiv][6];
  std::set<int>headx,heady,headz;

  if( readflag == CIO_E_GV_SAME ) {
    for(int i=0; i<RankInfo.size(); i++) {
      HeadTail[i][0]=RankInfo[i].HeadIndex[0];
      HeadTail[i][1]=RankInfo[i].HeadIndex[1];
      HeadTail[i][2]=RankInfo[i].HeadIndex[2];
      HeadTail[i][3]=RankInfo[i].TailIndex[0];
      HeadTail[i][4]=RankInfo[i].TailIndex[1];
      HeadTail[i][5]=RankInfo[i].TailIndex[2];

      headx.insert(HeadTail[i][0]);
      heady.insert(HeadTail[i][1]);
      headz.insert(HeadTail[i][2]);

    }
  }else if( readflag == CIO_E_GVX2_SAME ) {
    for(int i=0; i<RankInfo.size(); i++) {
      HeadTail[i][0]=RankInfo[i].HeadIndex[0]*2-1;
      HeadTail[i][1]=RankInfo[i].HeadIndex[1]*2-1;
      HeadTail[i][2]=RankInfo[i].HeadIndex[2]*2-1;
      HeadTail[i][3]=RankInfo[i].TailIndex[0]*2;
      HeadTail[i][4]=RankInfo[i].TailIndex[1]*2;
      HeadTail[i][5]=RankInfo[i].TailIndex[2]*2;

      headx.insert(HeadTail[i][0]);
      heady.insert(HeadTail[i][1]);
      headz.insert(HeadTail[i][2]);

    }
  }

  headT mapx,mapy,mapz;
  CreateHeadMap(headx,mapx);
  CreateHeadMap(heady,mapy);
  CreateHeadMap(headz,mapz);
/*
  for(int i=0; i<DFI_Domain.GlobalDivision[0]; i++) {
    int rank = DFI_rankMap[_IDX_IJK(i,0,0,DFI_Domain.GlobalDivision[0],
                                          DFI_Domain.GlobalDivision[1],
                                          DFI_Domain.GlobalDivision[2],0)];
    if( rank<0 ) continue;
    if( head[0] >= HeadTail[rank][0] ) StartEnd[0]=i;
    if( tail[0] >= HeadTail[rank][0] ) StartEnd[3]=i;
  }

  for(int j=0; j<DFI_Domain.GlobalDivision[1]; j++) {
    int rank = DFI_rankMap[_IDX_IJK(0,j,0,DFI_Domain.GlobalDivision[0],
                                          DFI_Domain.GlobalDivision[1],
                                          DFI_Domain.GlobalDivision[2],0)];
    if( rank<0 ) continue;
    if( head[1] >= HeadTail[rank][1] ) StartEnd[1]=j;
    if( tail[1] >= HeadTail[rank][1] ) StartEnd[4]=j;
  }

  for(int k=0; k<DFI_Domain.GlobalDivision[2]; k++) {
    int rank = DFI_rankMap[_IDX_IJK(0,0,k,DFI_Domain.GlobalDivision[0],
                                          DFI_Domain.GlobalDivision[1],
                                          DFI_Domain.GlobalDivision[2],0)];
    if( rank<0 ) continue;
    if( head[2] >= HeadTail[rank][2] ) StartEnd[2]=k;
    if( tail[2] >= HeadTail[rank][2] ) StartEnd[5]=k;
  }
*/

//x方向の絞り込み
  for( std::set<int>::iterator it=headx.begin();it!=headx.end();it++ )
  {
    if( head[0] >= *it ) StartEnd[0] = mapx.find(*it)->second; 
    if( tail[0] >= *it ) StartEnd[3] = mapx.find(*it)->second;
  }

//y方向の絞り込み
  for( std::set<int>::iterator it=heady.begin();it!=heady.end();it++ )
  {
    if( head[1] >= *it ) StartEnd[1] = mapy.find(*it)->second; 
    if( tail[1] >= *it ) StartEnd[4] = mapy.find(*it)->second;
  }

//z方向の絞り込み
  for( std::set<int>::iterator it=headz.begin();it!=headz.end();it++ )
  {
    if( head[2] >= *it ) StartEnd[2] = mapz.find(*it)->second; 
    if( tail[2] >= *it ) StartEnd[5] = mapz.find(*it)->second;
  }

  rankList.clear();

  //int dfi_head[3],dfi_tail[3];
  int sta_x,end_x,sta_y,end_y,sta_z,end_z;

  for(int k=StartEnd[2]; k<=StartEnd[5]; k++) {
  for(int j=StartEnd[1]; j<=StartEnd[4]; j++) {
  for(int i=StartEnd[0]; i<=StartEnd[3]; i++) {
    int rank = DFI_rankMap[_IDX_IJK(i,j,k,DFI_Domain.GlobalDivision[0],
                                          DFI_Domain.GlobalDivision[1],
                                          DFI_Domain.GlobalDivision[2],0)];
    if( rank<0 ) continue;

    sta_x = max(head[0],HeadTail[rank][0]);
    end_x = min(tail[0],HeadTail[rank][3]);
    sta_y = max(head[1],HeadTail[rank][1]);
    end_y = min(tail[1],HeadTail[rank][4]);
    sta_z = max(head[2],HeadTail[rank][2]);
    end_z = min(tail[2],HeadTail[rank][5]);

    if( sta_x <= end_x && sta_y <= end_y && sta_z <= end_z ) rankList.push_back(rank);
  }}}

  if( HeadTail ) delete [] HeadTail;

  if( headx.size()>0 ) headx.clear();
  if( heady.size()>0 ) heady.clear();
  if( headz.size()>0 ) headz.clear();
  if( mapx.size()>0 ) mapx.clear();
  if( mapy.size()>0 ) mapy.clear();
  if( mapz.size()>0 ) mapz.clear();


/*
  for(int i=0; i<RankInfo.size(); i++) {

    if( readflag == CIO_E_GV_SAME ) {
      for(int n=0; n<3; n++) {
        dfi_head[n]=RankInfo[i].HeadIndex[n];
        dfi_tail[n]=RankInfo[i].TailIndex[n];
      }
    }else if( readflag == CIO_E_GVX2_SAME ) {
      for(int n=0; n<3; n++) {
        dfi_head[n]=RankInfo[i].HeadIndex[n]*2-1;
        dfi_tail[n]=RankInfo[i].TailIndex[n]*2;
      }
    }

    //x 方向のスタートエンド
    sta_x=max(head[0],dfi_head[0]);
    end_x=min(tail[0],dfi_tail[0]);

    //y 方向のスタートエンド
    sta_y=max(head[1],dfi_head[1]);
    end_y=min(tail[1],dfi_tail[1]);

    //z 方向のスタートエンド
    sta_z=max(head[2],dfi_head[2]);
    end_z=min(tail[2],dfi_tail[2]);

    if( sta_x <= end_x && sta_y <= end_y && sta_z <= end_z ) rankList.push_back(i);
  }
*/

}

// #################################################################
// ファイル名を作成
std::string cio_DFI::Generate_FileName(int RankID, int step, const bool mio)
{

  if( DFI_Finfo.DirectoryPath.empty() ) return NULL;
  if( DFI_Finfo.Prefix.empty() ) return NULL;

  string fmt;
  if( DFI_Finfo.FileFormat == E_CIO_FMT_SPH ) {
    fmt=D_CIO_EXT_SPH;
  } else if( DFI_Finfo.FileFormat == E_CIO_FMT_BOV ) {
    fmt=D_CIO_EXT_BOV;
  }

  int len = DFI_Finfo.DirectoryPath.size() + DFI_Finfo.Prefix.size() + fmt.size() + 25; 
  // id(6) + step(10) + 1(\0) + "_"(2) + "."(1)+"id"(2)
  if( m_outSlice ) len += 11;

  char* tmp = new char[len];
  memset(tmp, 0, sizeof(char)*len);

  if( mio ) {
    if( m_outSlice ) {
      sprintf(tmp, "%s/%010d/%s_%010d_id%06d.%s",DFI_Finfo.DirectoryPath.c_str(),step,DFI_Finfo.Prefix.c_str(), 
            step,RankID,fmt.c_str());
    } else {
      sprintf(tmp, "%s/%s_%010d_id%06d.%s",DFI_Finfo.DirectoryPath.c_str(),DFI_Finfo.Prefix.c_str(), 
            step,RankID,fmt.c_str());
    }
  } else {
    if( m_outSlice ) {
      sprintf(tmp, "%s/%010d/%s_%010d.%s",DFI_Finfo.DirectoryPath.c_str(),step,DFI_Finfo.Prefix.c_str(), 
            step,fmt.c_str());
    } else {
      sprintf(tmp, "%s/%s_%010d.%s",DFI_Finfo.DirectoryPath.c_str(),DFI_Finfo.Prefix.c_str(), 
            step,fmt.c_str());
    }
  }
  
  std::string fname(tmp);
  if( tmp ) delete [] tmp;

  return fname;
}


// #################################################################
// ディレクトリがなければ作成、既存なら何もしない
int cio_DFI::MakeDirectory(string path)
{
  // 標準ライブラリ呼び出し
  // パーミッションはumaskで指定
  umask(022);
  int ret;

  ret = mkdir(path.c_str(), 0777); // rwx

  if ( 0 != ret )
  {
    // 既存以外のエラー
    if ( EEXIST != errno )
    {
      printf( "\tError(errno)=[%s]\n", strerror(errno) );
      return 0;
    }
  }

  return 1;
}

// #################################################################
// ディレクトリがなければ作成、既存なら何もしない
int cio_DFI::MakeDirectory(string path, int step)
{
  // 標準ライブラリ呼び出し
  // パーミッションはumaskで指定
  umask(022);
  int ret;

  if( m_outSlice ) {
    int len = DFI_Finfo.DirectoryPath.size() + 11;
    char* tmp = new char[len];
    memset(tmp, 0, sizeof(char)*len);
    sprintf(tmp, "%s/%010d",DFI_Finfo.DirectoryPath.c_str(),step);
    ret = mkdir(tmp, 0777); // rwx
  } else { 
    ret = mkdir(path.c_str(), 0777); // rwx
  }

  if ( 0 != ret )
  {
    // 既存以外のエラー
    if ( EEXIST != errno )
    {
      printf( "\tError(errno)=[%s]\n", strerror(errno) );
      return 0;
    }
  }

  return 1;
}

// #################################################################
// ActiveSubdomainファイルのエンディアンチェック
cio_EMatchType cio_DFI::isMatchEndianSbdmMagick( int ident )
{
  char magick_c[] = "SBDM";
  int  magick_i=0;

  //cheak match
  magick_i = (magick_c[3]<<24) + (magick_c[2]<<16) + (magick_c[1]<<8) + magick_c[0];
  if( magick_i == ident )
  {
    return CIO_Match;
  }

  //chack unmatch
  magick_i = (magick_c[0]<<24) + (magick_c[1]<<16) + (magick_c[2]<<8) + magick_c[3];
  if( magick_i == ident )
  {
    return CIO_UnMatch;
  }

  //unknown format
  return CIO_UnKnown;
}

// #################################################################
// ActiveSubdomainファイルの読み込み(static関数)
cio_ErrorCode cio_DFI::ReadActiveSubdomainFile( std::string subDomainFile,
    std::vector<cio_ActiveSubDomain>& subDomainInfo, int div[3] )
{
  if( subDomainFile.empty() ) return CIO_ERROR_OPEN_SBDM;

  // ファイルオープン
  FILE*fp = fopen( subDomainFile.c_str(), "rb" );
  if( !fp ) return CIO_ERROR_OPEN_SBDM;

  //エンディアン識別子
  int ident;
    if( fread( &ident, sizeof(int), 1, fp ) != 1 )
  {
    fclose(fp);
    return CIO_ERROR_READ_SBDM_HEADER;
  }

  cio_EMatchType endian = isMatchEndianSbdmMagick( ident );
  if( endian == CIO_UnKnown )
  {
    fclose(fp);
    return CIO_ERROR_READ_SBDM_FORMAT;
  }

  // 領域分割数
  if( fread( div, sizeof(int), 3, fp ) != 3 )
  {
    fclose(fp);
    return CIO_ERROR_READ_SBDM_DIV;
  }
  if( endian == CIO_UnMatch ) BSWAPVEC(div,3);

  // contents
  size_t nc = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
  unsigned char *contents = new unsigned char[nc];
  if( fread( contents, sizeof(unsigned char), nc, fp ) != nc )
  {
    delete [] contents;
    fclose(fp);
    return CIO_ERROR_READ_SBDM_CONTENTS;
  }

  // ファイルクローズ
  fclose(fp);

  size_t ptr = 0;
  // 活性ドメイン情報の生成
  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    if( contents[ptr] == 0x01 )
    {
      int pos[3] = {i,j,k};
      cio_ActiveSubDomain dom( pos );
      subDomainInfo.push_back(dom);
    }
    ptr++;
  }}}

  // contentsのdelete
  delete [] contents;

  // 活性ドメインの数をチェック
  if( subDomainInfo.size() == 0 )
  {
    return CIO_ERROR_SBDM_NUMDOMAIN_ZERO;
  }

  return CIO_SUCCESS;
}

// #################################################################
// 領域情報のチェック
cio_ErrorCode cio_DFI::CheckData( int nRank , int div[3] ,
              std::vector<cio_ActiveSubDomain>& subDomainInfo)
{

  // 領域分割数
  if( div[0] <= 0 || div[1] <= 0 || div[2] <= 0 ) 
    return CIO_ERROR_INVALID_DIVNUM;

  //活性サブドメイン情報
  int ndom = subDomainInfo.size();
  if( ndom == 0 ) {
    //活性サブドメイン情報が空のとき、全領域を活性サブドメインとする
    if( nRank != div[0]*div[1]*div[2] )
    {
      //return CIO_ERROR_MISMATCH_NP_SUBDOMAIN;
    }
    for( int k=0;k<div[2];k++ ){
    for( int j=0;j<div[1];j++ ){
    for( int i=0;i<div[0];i++ ){
      int pos[3] = {i,j,k};
      cio_ActiveSubDomain dom( pos );
      subDomainInfo.push_back(dom);
    }}}
  } else {
    if( nRank != ndom )
    {
      return CIO_ERROR_MISMATCH_NP_SUBDOMAIN;
    }
  }


  return CIO_SUCCESS;
}

// #################################################################
// ランクマップを生成 （非活性を含む）
int* cio_DFI::CreateRankMap(int div[3],std::vector<cio_ActiveSubDomain> subDomainInfo)
{
  
  size_t ndiv = size_t(div[0]) * size_t(div[1]) * size_t(div[2]);
  int *rankMap = new int[ndiv];
  if( !rankMap ) return NULL;

  for( size_t i=0;i<ndiv;i++ ) rankMap[i] = -1;

  // 活性サブドメイン情報配置位置に0をセット
  for( int i=0; i<subDomainInfo.size(); i++ )
  {
    //サブドメイン情報
    const cio_ActiveSubDomain* dom = &subDomainInfo[i];
    if( !dom ) 
    {
      delete [] rankMap;
      return NULL;
    }

    //位置を取得
    const int *pos = dom->GetPos();
    if( !pos )
    {
      delete [] rankMap;
      return NULL;
    }

    //0をセット
    rankMap[_IDX_IJK(pos[0],pos[1],pos[2],div[0],div[1],div[2],0)] = 0;
  }

 //i->j->kの優先順で活性サブドメインにランク番号をセット
  int rankCount = 0;
  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    if( rankMap[_IDX_IJK(i,j,k,div[0],div[1],div[2],0)] == 0 )
    {
      rankMap[_IDX_IJK(i,j,k,div[0],div[1],div[2],0)] = rankCount;
      //rankCount++;
    }
    rankCount++;
  }}}

  //for(int i=0;i<ndiv;i++) printf("rankMap[%d] : %d\n",i,rankMap[i]);

  return rankMap;
}

// #################################################################
// ランクマップを生成 
int* cio_DFI::CreateActiveRankMap(vector<cio_Rank> &RankInfo, int (*HeadTail)[6],
                                  headT mapX, headT mapY, headT mapZ)
{

  int i,j,k;

  int div[3];
  div[0] = mapX.size();
  div[1] = mapY.size();
  div[2] = mapZ.size();

  size_t ndiv = div[0]*div[1]*div[2];

  //printf("ndiv : %d\n",ndiv);

  int *rankMap = new int[ndiv];
  for(int i=0; i<ndiv; i++) rankMap[i]=-1;

  headT::iterator it;

  for(int n=0; n<RankInfo.size(); n++)
  {
    it=mapX.find(RankInfo[n].HeadIndex[0]);
    i=it->second;

    it=mapY.find(RankInfo[n].HeadIndex[1]);
    j=it->second;

    it=mapZ.find(RankInfo[n].HeadIndex[2]);
    k=it->second;

    int rnkPos=_IDX_IJK(i,j,k,div[0],div[1],div[2],0);

    //printf("i : %d j : %d k : %d rnkPos : %d RankID : %d\n",i,j,k,rnkPos,RankInfo[n].RankID);

    //rankMap[_IDX_IJK(i,j,k,div[0],div[1],div[2],0)] = RankInfo[n].RankID;    
    rankMap[_IDX_IJK(i,j,k,div[0],div[1],div[2],0)] = n;    
  }

  return rankMap;

}

// #################################################################
// head&tail の生成
bool cio_DFI::SetHeadTail( vector<cio_Rank> RankInfo, int (*HeadTail)[6], 
                              std::set<int>&headx, std::set<int>&heady, std::set<int>&headz ) 
{

  for(int i=0; i<RankInfo.size(); i++ ) {
   HeadTail[i][0] = RankInfo[i].HeadIndex[0];
   HeadTail[i][1] = RankInfo[i].HeadIndex[1];
   HeadTail[i][2] = RankInfo[i].HeadIndex[2];
   HeadTail[i][3] = RankInfo[i].TailIndex[0];
   HeadTail[i][4] = RankInfo[i].TailIndex[1];
   HeadTail[i][5] = RankInfo[i].TailIndex[2];

   headx.insert(HeadTail[i][0]);
   heady.insert(HeadTail[i][1]);
   headz.insert(HeadTail[i][2]);

  }

  return true;

}


// #################################################################
// head&tail の生成
bool cio_DFI::CreateHeadTail(int div[3], int* rankMap, int gvox[3], int (*HeadTail)[6], vector<cio_Rank> &RankInfo, cio_Rank rank) 
{

  if( !rankMap ) return false;
 
  //ローカルのVOXEL数
  int *nvX = new int[div[0]];
  int *nvY = new int[div[1]];
  int *nvZ = new int[div[2]];
  int *nv[3] = {nvX,nvY,nvZ};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    //基準のボクセル数
    int nbase = gvox[n] / div[n];

    //余り
    int amari = gvox[n] % div[n];

    //ボクセル数をセット
    for( int i=0;i<div[n];i++ )
    {
      nvd[i] = nbase;
      if( i<amari ) nvd[i]++;
    }
  }

  //head
  int *headX = new int[div[0]];
  int *headY = new int[div[1]];
  int *headZ = new int[div[2]];
  int *head[3] = {headX,headY,headZ};
  for( int n=0;n<3;n++ )
  {
    int *nvd = nv[n];
    int *hd = head[n];
    hd[0] = 1;

    for( int i=1;i<div[n];i++ )
    {
      hd[i] = hd[i-1]+nvd[i-1];
    }
  }

  for( int k=0;k<div[2];k++ ){
  for( int j=0;j<div[1];j++ ){
  for( int i=0;i<div[0];i++ ){
    int rankNo = rankMap[_IDX_IJK(i,j,k,div[0],div[1],div[2],0)];
    if( rankNo < 0 ) continue;
    HeadTail[rankNo][0] = headX[i];
    HeadTail[rankNo][1] = headY[j];
    HeadTail[rankNo][2] = headZ[k];
    HeadTail[rankNo][3] = headX[i]+nvX[i]-1;
    HeadTail[rankNo][4] = headY[j]+nvY[j]-1;
    HeadTail[rankNo][5] = headZ[k]+nvZ[k]-1;

    rank.RankID = rankNo;
    rank.VoxelSize[0]=HeadTail[rankNo][3]-HeadTail[rankNo][0]+1;
    rank.VoxelSize[1]=HeadTail[rankNo][4]-HeadTail[rankNo][1]+1;
    rank.VoxelSize[2]=HeadTail[rankNo][5]-HeadTail[rankNo][2]+1;
    rank.HeadIndex[0]=HeadTail[rankNo][0];
    rank.HeadIndex[1]=HeadTail[rankNo][1];
    rank.HeadIndex[2]=HeadTail[rankNo][2];
    rank.TailIndex[0]=HeadTail[rankNo][3];
    rank.TailIndex[1]=HeadTail[rankNo][4];
    rank.TailIndex[2]=HeadTail[rankNo][5];
    RankInfo.push_back(rank); 

  }}}

  delete [] nvX;
  delete [] nvY;
  delete [] nvZ;
  delete [] headX;
  delete [] headY;
  delete [] headZ; 

  return true;

}

// #################################################################
// head map の生成
bool cio_DFI::CreateHeadMap(std::set<int>head, headT &map)
{

  int cnt=0;
  for(std::set<int>::iterator it=head.begin();it!=head.end();it++)
  {
    int key=*it;
    map.insert(headT::value_type(key,cnt));
    cnt++;
  }
  return true;
}

// #################################################################
// set Unit
void cio_DFI::SetUnitLength(bool out_length, string Length, double L0) 
{

  DFI_Unit.out_Length = out_length;
  DFI_Unit.Length = Length;
  DFI_Unit.L0 = L0;

}

void cio_DFI::SetUnitVelo(bool out_Velocity, string Velocity, double V0)
{

  DFI_Unit.out_Velocity = out_Velocity;
  DFI_Unit.Velocity = Velocity;
  DFI_Unit.V0 = V0;

}
void cio_DFI::SetUnitPres(bool out_Pressure, string Pressure, double P0, double DiffPrs)
{

  DFI_Unit.out_Pressure = out_Pressure;
  DFI_Unit.Pressure =Pressure;
  DFI_Unit.P0 = P0;
  DFI_Unit.DiffPrs = DiffPrs;

}
void cio_DFI::SetUnitTemp(bool out_Temp, string Temp, double Btemp, double DiffTemp)
{

  DFI_Unit.out_Temperature = out_Temp;
  DFI_Unit.Temperatur = Temp;
  DFI_Unit.BaseTemp = Btemp;
  DFI_Unit.DiffTemp = DiffTemp;

}

// #################################################################
// debug write index.dfi & proc.dfi
bool cio_DFI::dbwrite(int RankID)
{

  if( RankID != 0 ) return true;

  printf("**** input index.dfi file infomation ****\n"); 
  printf("FileInfo {\n");
  printf("  DirectoryPath               :%s\n",DFI_Finfo.DirectoryPath.c_str()); 
  printf("  Prefix                      :%s\n",DFI_Finfo.Prefix.c_str()); 
  printf("  FileFormat                  :%d\n",DFI_Finfo.FileFormat); 
  printf("  GuideCell                   :%d\n",DFI_Finfo.GuideCell); 
  printf("  DataType                    :%s\n",DFI_Finfo.DataType.c_str()); 
  printf("  Endian                      :%s\n",DFI_Finfo.Endian.c_str()); 
  printf("  ArrayShape                  :%s\n",DFI_Finfo.ArrayShape.c_str()); 
  printf("  Component                   :%d\n",DFI_Finfo.Component); 
  printf("}\n");
  printf("\n");

  printf("FilePath {\n");
  printf("  Process                     :%s\n",DFI_Fpath.Process.c_str());
  printf("}\n");
  printf("\n");

  printf("Unit {\n");
  printf("  Length                      :%s\n",DFI_Unit.Length.c_str());
  printf("  L0                          :%e\n",DFI_Unit.L0);
  printf("  Velocity                    :%s\n",DFI_Unit.Velocity.c_str());
  printf("  V0                          :%e\n",DFI_Unit.V0);
  printf("  Pressure                    :%s\n",DFI_Unit.Pressure.c_str());
  printf("  P0                          :%e\n",DFI_Unit.P0);
  printf("  DiffPrs                     :%e\n",DFI_Unit.DiffPrs);
  printf("  Temperatur                  :%s\n",DFI_Unit.Temperatur.c_str());
  printf("  BaseTemp                    :%e\n",DFI_Unit.BaseTemp);
  printf("  DiffTemp                    :%e\n",DFI_Unit.DiffTemp);
  printf("}\n");
  printf("\n");

  printf("TimeSlice {\n");
  for(int i=0; i<TimeSlice.size(); i++) {
    printf("  Slice[%d] {\n",i);
    printf("    Step                        :%d\n",TimeSlice[i].step);
    printf("    Time                        :%e\n",TimeSlice[i].time);
    for(int j=0; j<TimeSlice[i].Min.size(); j++) {
      printf("    MinMax[%d] {\n",j+1);
      printf("      Min                         :%e\n",TimeSlice[i].Min[j]);
      printf("      Max                         :%e\n",TimeSlice[i].Max[j]);
      printf("    }\n");
    }
    printf("  }\n");
  }
  printf("}\n");
  printf("\n");

  printf("**** input proc.dfi file infomation ****\n"); 
  printf("Domian {\n");
  printf("  GlovalOrigin                :%e %e %e\n",DFI_Domain.GlobalOrigin[0],
                                                     DFI_Domain.GlobalOrigin[1],
                                                     DFI_Domain.GlobalOrigin[2]); 
  printf("  GlovalRegion                :%e %e %e\n",DFI_Domain.GlobalRegion[0],
                                                     DFI_Domain.GlobalRegion[1],
                                                     DFI_Domain.GlobalRegion[2]); 
  printf("  GlobalVoxel                 :%d %d %d\n",DFI_Domain.GlobalVoxel[0],
                                                     DFI_Domain.GlobalVoxel[1],
                                                     DFI_Domain.GlobalVoxel[2]); 
  printf("  GlobalDivision              :%d %d %d\n",DFI_Domain.GlobalDivision[0],
                                                     DFI_Domain.GlobalDivision[1],
                                                     DFI_Domain.GlobalDivision[2]); 
  printf("  ActiveSubdomain             :%s\n",      DFI_Domain.ActiveSubdomain.c_str());
  printf("}\n");
  printf("\n");

  printf("MPI {\n");
  printf("  NumberofRank                :%d\n",DFI_MPI.NumberOfRank); 
  printf("  NumberofGroup               :%d\n",DFI_MPI.NumberOfGroup); 
  printf("}\n");
  printf("\n");

  printf("Process {\n");
  for (int i=0; i<RankInfo.size(); i++) {
    printf("  Rank[%d] {\n",i);
    printf("    ID                          :%d\n",RankInfo[i].RankID);
    printf("    HostName                    :%s\n",RankInfo[i].HostName.c_str()); 
    printf("    VoxelSize                   :%d %d %d\n",RankInfo[i].VoxelSize[0],RankInfo[i].VoxelSize[1],RankInfo[i].VoxelSize[2]); 
    printf("    HeadIndex                   :%d %d %d\n",RankInfo[i].HeadIndex[0],RankInfo[i].HeadIndex[1],RankInfo[i].HeadIndex[2]); 
    printf("    TailIndex                   :%d %d %d\n",RankInfo[i].TailIndex[0],RankInfo[i].TailIndex[1],RankInfo[i].TailIndex[2]); 
    printf("  }\n");
  }
  printf("}\n");

  return true;

}

