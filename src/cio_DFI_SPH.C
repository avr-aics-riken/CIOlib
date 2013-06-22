/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_SPH.C
 * @brief  cio_DFI_SPH Class
 * @author kero    
 */

#include "cio_DFI.h"
#include "cio_DFI_SPH.h"

// #################################################################
// コンストラクタ
cio_DFI_SPH::cio_DFI_SPH()
{

}


// #################################################################
// デストラクタ
cio_DFI_SPH::~cio_DFI_SPH()
{

}

// #################################################################
// 
void* cio_DFI_SPH::ReadData(int step, int gc,
                           int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                           double &time,
                           const bool mode, unsigned &step_avr, double &time_avr)
{

  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);

  int sz[3];
  for(int i=0; i<3; i++) sz[i]=tail[i]-head[i]+1+2*gc;
  int val_size = sz[0]*sz[1]*sz[2];

  void* val;

  if       ( Dtype == E_CIO_FLOAT32 ) {
      val = new float[val_size];  
      memset(val, 0, sizeof(float)*val_size);
  } else if( Dtype == E_CIO_FLOAT64 ) { 
      val = new double[val_size];  
      memset(val, 0, sizeof(double)*val_size);
  } else return NULL;

  ReadData(step, gc, Gvoxel, Gdivision, head, tail, val, time, mode, step_avr, time_avr);
  return val;
}

// #################################################################
// 
void cio_DFI_SPH::ReadData(int step, int gc, 
                           int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                           void *val, double &time,
                           const bool mode, unsigned &step_avr, double &time_avr)
{
  bool mio = false;
  
  if( RankInfo.size() > 1 ) mio = true;

/*
  cio_Domain G_Domain;
  vector<cio_Rank> G_RankInfo;
  cio_Rank G_Rank;
  cio_Create_Domain(m_comm,Gvoxel, Gdivision, head, tail,G_Domain,G_RankInfo,G_Rank);
*/

  cio_EGlobalVoxel readflag;  //CIO_E_GV_SAME:密　CIO_E_GVX2_SAME:粗 CIO_E_OTHER:その他
  readflag = CheckGlobalVoxel(Gvoxel, DFI_Domain.GlobalVoxel);

//読込みランクリストの生成
  vector<int> RList;
  CreateRankList(head, tail, gc, readflag, RList);

//密データの読込み処理
  if( readflag == CIO_E_GV_SAME ) { 
    for(int i=0; i<RList.size(); i++) {
      int n = RList[i];
      int ID= RankInfo[n].RankID;
      std::string fname = Generate_FileName(ID,step,mio);

      //printf("fname : %s\n",fname.c_str());

      int sta[3],end[3];
      if( CheckReadArea(head, tail, gc, RankInfo[n].HeadIndex, RankInfo[n].TailIndex,
                        DFI_Finfo.GuideCell, readflag, sta, end) )
      { 
        //1対1の読込み
        if( !read_S3D(fname.c_str(), step, DFI_Finfo.GuideCell, val, time, 
                      mode, step_avr, time_avr ) ) {
          val = NULL;
          printf("**** error read s3d ****");
          return;
        }

      } else {
        //1対多の読込み
        if( !read_MxN(fname.c_str(), step, gc, head, tail, RankInfo[n].HeadIndex,
                    RankInfo[n].TailIndex, sta, end, val, time, n, mode, step_avr, time_avr) ) {
          val = NULL;
          printf("**** error read mxn ID %d dfi_ID : %d fname : %s\n",m_RankID,n,
                  fname.c_str());
          return;
        }
      }
      if( m_RankID == 0 ) {
        printf("\t[%s] has read :\tstep=%d  time=%e ]\n",fname.c_str(), step, time);
      }
    }
    return;
//粗い粒子の処理
  }else if( readflag == CIO_E_GVX2_SAME ) { 

    for(int i=0; i<RList.size(); i++) {
      std::string fname = Generate_FileName(RList[i],step,mio);

      //printf("fname : %s\n",fname.c_str());

      int sta[3],end[3],dfi_head[3],dfi_tail[3];

      int n = RList[i];

      for(int j=0; j<3; j++){
        dfi_head[j]=RankInfo[n].HeadIndex[j]*2-1;
        dfi_tail[j]=RankInfo[n].TailIndex[j]*2;
      }
      //1対1の読込み&補間
      if( CheckReadArea(head, tail, gc, dfi_head, dfi_tail,DFI_Finfo.GuideCell, readflag, 
                        sta, end)) 
      {
        if( !read_Coarse(fname.c_str(), step, gc, head, tail, Gvoxel, RankInfo[n].HeadIndex,
                          RankInfo[n].TailIndex, RankInfo[n].VoxelSize, DFI_Finfo.GuideCell,
                          DFI_Finfo.Component, sta, end, val, time, mode, step_avr, time_avr )) {
           val = NULL;
           return ;
        }
      } else {
      //1対多の読込み&補間
        if( !read_Coarse_MxN(fname.c_str(), step, gc, head, tail, Gvoxel, RankInfo[n].HeadIndex,
                          RankInfo[n].TailIndex, RankInfo[n].VoxelSize, DFI_Finfo.GuideCell,
                          DFI_Finfo.Component, sta, end, val, time, n, mode, step_avr, time_avr ) ) {
           val = NULL;
           return ;
        }
      }
      if( m_RankID == 0 ) {
        printf("\t[%s] has read :\tstep=%d  time=%e ]\n",fname.c_str(), step, time);
      }
    }
  } else {
    printf("Error field data format\n");
    val = NULL;
    return;
  }

}

// #################################################################
// 
bool cio_DFI_SPH::read_Head(FILE* fp, int Etype, int step, int voxsize[3], double &time, 
                           RealType &real_type)
{

  unsigned int dmy,type_dmy;
  
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != 8 ) { fclose(fp); return false; }

  DataDims data_dims;
  if( fread(&data_dims, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(data_dims); 
  if( data_dims == _SCALAR && DFI_Finfo.Component != 1 ) { fclose(fp); return false; } 
  if( data_dims == _VECTOR && DFI_Finfo.Component <= 1 ) { fclose(fp); return false; } 

  //RealType real_type;
  if( fread(&real_type, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(real_type); 
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != 8 ) { fclose(fp); return false; }

  if( real_type == _FLOAT ) type_dmy=12;
  if( real_type == _DOUBLE) type_dmy=24;

//ボクセルサイズ
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) {
    if( fread(voxsize, sizeof(int), 3, fp) != 3 ){fclose(fp);return false;}
    if( Etype == CIO_UnMatch ) {
      BSWAP32(voxsize[0]); 
      BSWAP32(voxsize[1]); 
      BSWAP32(voxsize[2]);
    }
  } else if( real_type == _DOUBLE ) {
    long long tmp[3];
    if( fread(tmp, sizeof(long long), 3, fp) != 3 ){fclose(fp);return false;}
    if( Etype == CIO_UnMatch ) {
      BSWAP64(tmp[0]); 
      BSWAP64(tmp[1]); 
      BSWAP64(tmp[2]);
    }
    voxsize[0]=(int)tmp[0];
    voxsize[1]=(int)tmp[1];
    voxsize[2]=(int)tmp[2];
  } 
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }

//原点座標
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) {
    float voxorg[3];
    if( fread(voxorg, sizeof(float), 3, fp) != 3 ){fclose(fp);return false;}
    if( Etype == CIO_UnMatch ) {
      BSWAP32(voxorg[0]); 
      BSWAP32(voxorg[1]); 
      BSWAP32(voxorg[2]); 
    }
  } else if( real_type == _DOUBLE ) {
    double voxorg[3];
    if( fread(voxorg, sizeof(double), 3, fp) != 3 ){fclose(fp);return false;}
    if( Etype == CIO_UnMatch ) {
      BSWAP64(voxorg[0]); 
      BSWAP64(voxorg[1]); 
      BSWAP64(voxorg[2]); 
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }
  
//pit
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) {
    float voxpit[3];
    if( fread(voxpit, sizeof(float), 3, fp) != 3 ){fclose(fp);return false;}
    if( Etype == CIO_UnMatch ) {
      BSWAP32(voxpit[0]); 
      BSWAP32(voxpit[1]); 
      BSWAP32(voxpit[2]); 
    }
  } else if( real_type == _DOUBLE ) {
    double voxpit[3];
    if( fread(voxpit, sizeof(double), 3, fp) != 3 ){fclose(fp);return false;}
    if( Etype == CIO_UnMatch ) {
      BSWAP64(voxpit[0]); 
      BSWAP64(voxpit[1]); 
      BSWAP64(voxpit[2]); 
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }

//step,time
  if( real_type == _FLOAT ) type_dmy = 8;
  if( real_type == _DOUBLE) type_dmy = 16;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }
  if( real_type == _FLOAT ) { 
    int r_step;
    if( fread(&r_step, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP32(r_step);
    if( r_step != step ) { fclose(fp); return false; }
  } else if( real_type == _DOUBLE ) {
    long long r_step;
    if( fread(&r_step, sizeof(long long), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP64(r_step);
    if( r_step != step ) { fclose(fp); return false; }
  }
  if( real_type == _FLOAT ) {
    float r_time;
    if( fread(&r_time, sizeof(float), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP32(r_time); 
    time = r_time;
  } else if( real_type == _DOUBLE ) {
    double r_time;
    if( fread(&r_time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP64(r_time); 
    time = r_time;
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy); 
  if( dmy != type_dmy ) { fclose(fp); return false; }


return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_S3D(const char* fname, int step, int gc, void*val, 
                           double &time, 
                           const bool mode, unsigned &step_avr, double &time_avr)
{

  if( !fname || !DFI_Finfo.Component ) return false;

  RealType real_type;

  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
  int Dsize = get_cio_Datasize(Dtype);
  //printf("**** Dtype : %d Dsize : %d\n",Dtype,Dsize);

  //Endian check....
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  string Endian="";
  if( cdumy[0] == 0x01 ) Endian = "little";
  if( cdumy[0] == 0x00 ) Endian = "big";

  int Etype = CIO_UnKnown;
  if( DFI_Finfo.Endian == Endian ) {
    Etype = CIO_Match;
  }else {
    Etype = CIO_UnMatch;
  }

  if( Etype == CIO_UnKnown ) return false;

  FILE* fp;
  if( !(fp=fopen(fname,"rb")) ) {
    printf("Can't open file. (%s)\n",fname);
    //return false;
    exit(0);
  }


  int voxsize[3];
  if( !read_Head(fp, Etype, step, voxsize, time, real_type) ) { fclose(fp); return false; }

  unsigned int dmy, type_dmy;
//data 
  int arraylen = voxsize[0]*voxsize[1]*voxsize[2]*DFI_Finfo.Component;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy);
  if( dmy != arraylen*Dsize ) { fclose(fp); return false; }
  if( fread(val, Dsize, arraylen, fp) != arraylen ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) {
    if( Dtype == E_CIO_FLOAT32) {
      float *val2 = (float *)val;
      BSWAPVEC(val2,arraylen);
    } else {
      double *val2 = (double *)val;
      DBSWAPVEC(val2,arraylen);
    }
  }
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy);
 
//read average step & time
  if( !mode ) {
    if( real_type == _FLOAT ) type_dmy = 8;
    if( real_type == _DOUBLE) type_dmy = 16;
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP32(dmy);
    if( dmy != type_dmy ) { fclose(fp); return false; }
    if( real_type == _FLOAT ) {
      int r_step;
      if( fread(&r_step, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP32(r_step);
      step_avr = (unsigned)r_step;
    } else if( real_type == _DOUBLE ) {
      long long r_step;
      if( fread(&r_step, sizeof(long long), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP64(r_step);
      step_avr = (unsigned)r_step;
    }
    if( real_type == _FLOAT ) {
      float r_time;
      if( fread(&r_time, sizeof(float), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP32(r_time);
      time_avr = (double)r_time;
    } else if( real_type == _DOUBLE ) {
      double r_time;
      if( fread(&r_time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP64(r_time);
      time_avr = r_time;
    }
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP32(dmy);
    if( dmy != type_dmy ) { fclose(fp); return false; }
  }

 // printf("***********************\n");
 // printf("*** read sph s3d ok ***\n");
 // printf("***********************\n");

  fclose(fp);
  return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_MxN(const char* fname, int step, int gc, int head[3], int tail[3], 
                           int R_head[3], int R_tail[3], int sta[3], int end[3], void* val, 
                           double &time, int n_dfi, 
                           const bool mode, unsigned &step_avr, double &time_avr)
{

  if( !fname || !DFI_Finfo.Component ) return false;

  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
  int Dsize = get_cio_Datasize(Dtype);
  //printf("**** Dtype : %d Dsize : %d\n",Dtype,Dsize);

  //Endian check....
  int idumy = 1;
  char* cdumy = (char*)(&idumy);
  string Endian="";
  if( cdumy[0] == 0x01 ) Endian = "little";
  if( cdumy[0] == 0x00 ) Endian = "big";

  int Etype = CIO_UnKnown;
  if( DFI_Finfo.Endian == Endian ) {
    Etype = CIO_Match;
  }else {
    Etype = CIO_UnMatch;
  }

  if( Etype == CIO_UnKnown ) return false;

  FILE* fp;
  if( !(fp=fopen(fname,"rb")) ) {
    printf("Can't open file. (%s)\n",fname);
    //return false;
    exit(0);
  }

  RealType real_type;
  int voxsize[3];
  if( !read_Head(fp, Etype, step, voxsize, time, real_type) ) { fclose(fp); return false; }

  int off_set;
  long off_setleng;
//ijkn
  if( !strcasecmp( DFI_Finfo.ArrayShape.c_str(), "ijkn" ) || DFI_Finfo.Component == 1 ) { 
    off_set = (head[2]-gc)-(R_head[2]-DFI_Finfo.GuideCell);
    off_setleng = (voxsize[0])*(voxsize[1])*(off_set);
  } else if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "nijk" ) ) {
//nijk
    off_set = (head[2]-gc)-(R_head[2]-DFI_Finfo.GuideCell);
    off_setleng = (voxsize[0])*(voxsize[1])*(off_set)*DFI_Finfo.Component;
  }

  int nkrec = end[2]-sta[2]+1;

  unsigned int dmy, type_dmy;
  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }

  long long idx,idx2;
  int imax,jmax,kmax;

  imax = tail[0]-head[0]+1;
  jmax = tail[1]-head[1]+1;
  kmax = tail[2]-head[2]+1;


  // ijkn
  if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || DFI_Finfo.Component == 1 ) {
    int arraylen = (voxsize[0])*(voxsize[1]);
    void *rbuff;
    if( Dtype == E_CIO_FLOAT32 ) {
      rbuff = new float[arraylen];
    } else {
      rbuff = new double[arraylen];
    }

    for(int ncomp=0; ncomp<DFI_Finfo.Component; ncomp++) {
      if( off_set > 0 ) {
        if( fseek(fp,off_setleng*Dsize,SEEK_CUR) != 0 ) 
        {
          printf("*** Error fseek length offset length : %d\n", off_set*Dsize);
          fclose(fp); return false; 
        }
      }

      int iimax,jjmax;
      iimax=voxsize[0]-2*DFI_Finfo.GuideCell;
      jjmax=voxsize[1]-2*DFI_Finfo.GuideCell;
      if( sta[2]>0 && R_head[2]>1) {
        for(int n=0; n<DFI_Finfo.GuideCell; n++) {
           if( fread(rbuff, Dsize, arraylen, fp) != arraylen ) { fclose(fp); return false; }
           if( Etype == CIO_UnMatch ) {
             if( Dtype == E_CIO_FLOAT32 ) {
               float *rbuff2 = (float *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else {
               double *rbuff2 = (double *)rbuff;
               DBSWAPVEC(rbuff2,arraylen);
             }
           }
        }
        nkrec-=DFI_Finfo.GuideCell;
      }

      for(int n=0; n<nkrec; n++) {
        if( fread(rbuff, Dsize, arraylen, fp) != arraylen ) { fclose(fp); return false; }
        if( Etype == CIO_UnMatch ) {
          if( Dtype == E_CIO_FLOAT32 ) {
            float *rbuff2 = (float *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else {
            double *rbuff2 = (double *)rbuff;
            DBSWAPVEC(rbuff2,arraylen);
          }
        }
        int k=sta[2]+n;
        for( int j=sta[1]; j<=end[1]; j++ ) {
        for( int i=sta[0]; i<=end[0]; i++ ) {
           idx = _IDX_IJKN(i-(head[0]-1),j-(head[1]-1),k-(head[2]-1),ncomp,imax,jmax,kmax,gc);
           idx2= _IDX_IJ(i-(R_head[0]-1),j-(R_head[1]-1),iimax,jjmax,DFI_Finfo.GuideCell);
           //val[idx]=rbuff[idx2];
           if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)rbuff,idx2);
           if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)rbuff,idx2);
        }}
      }
    }

    delete [] rbuff;

  // nijk
  } else if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "nijk" ) ) {

    int arraylen = (voxsize[0])*(voxsize[1])*DFI_Finfo.Component;
    void *rbuff;
    if( Dtype == E_CIO_FLOAT32 ) {
      rbuff = new float[arraylen];
    } else {
      rbuff = new double[arraylen];
    }

    if( off_set > 0 ) {
      if( fseek(fp,off_setleng*Dsize,SEEK_CUR) != 0 ) 
      {
        printf("*** Error fseek length offset length : %d\n", off_set*Dsize);
        fclose(fp); return false; 
      }
    }

    int iimax,jjmax;
    iimax=voxsize[0]-2*DFI_Finfo.GuideCell;
    jjmax=voxsize[1]-2*DFI_Finfo.GuideCell;
    if( sta[2]>0 && R_head[2]>1 ) {
      for(int n=0; n<DFI_Finfo.GuideCell; n++) {
         if( fread(rbuff, Dsize, arraylen, fp) != arraylen ) { fclose(fp); return false; }
         if( Etype == CIO_UnMatch ) {
           if( Dtype == E_CIO_FLOAT32 ) {
             float* rbuff2 = (float *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else {
             double* rbuff2 = (double *)rbuff;
             DBSWAPVEC(rbuff2,arraylen);
           }
         }
      }
      nkrec-=DFI_Finfo.GuideCell;
    }

    for(int n=0; n<nkrec; n++) {
      if( fread(rbuff, Dsize, arraylen, fp) != arraylen ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) {
        if( Dtype == E_CIO_FLOAT32 ) {
          float* rbuff2 = (float *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else {
          double* rbuff2 = (double *)rbuff;
          DBSWAPVEC(rbuff2,arraylen);
        }
      }
      int k=sta[2]+n;
      for( int j=sta[1]; j<=end[1]; j++ ) {
      for( int i=sta[0]; i<=end[0]; i++ ) {
      for(int ncomp=0; ncomp<DFI_Finfo.Component; ncomp++) {
         idx = _IDX_NIJK(ncomp,i-(head[0]-1),j-(head[1]-1),k-(head[2]-1),DFI_Finfo.Component,imax,jmax,kmax,gc);
         idx2= _IDX_NIJ(ncomp,i-(R_head[0]-1),j-(R_head[1]-1),iimax,jjmax,DFI_Finfo.Component,
                          DFI_Finfo.GuideCell);
         //val[idx]=rbuff[idx2];
         if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)rbuff,idx2);
         if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)rbuff,idx2);

      }}}
    }

    delete [] rbuff;

  }

  if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) BSWAP32(dmy);

//read average step & time
  if( !mode ) {
    if( real_type == _FLOAT ) type_dmy = 8;
    if( real_type == _DOUBLE) type_dmy = 16;
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP32(dmy);
    if( dmy != type_dmy ) { fclose(fp); return false; }
    if( real_type == _FLOAT ) {
      int r_step;
      if( fread(&r_step, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP32(r_step);
      step_avr=(unsigned)r_step;
    } else if( real_type == _DOUBLE ) {
      long long r_step;
      if( fread(&r_step, sizeof(long long), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP64(r_step);
      step_avr=(unsigned)r_step;
    }
    if( real_type == _FLOAT ) {
      float r_time;
      if( fread(&r_time, sizeof(float), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP32(r_time);
      time_avr = (double)r_time;
    } else if( real_type == _DOUBLE ) {
      double r_time;
      if( fread(&r_time, sizeof(double), 1, fp) != 1 ) { fclose(fp); return false; }
      if( Etype == CIO_UnMatch ) BSWAP64(r_time);
      time_avr = r_time;
    }
    if( fread(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return false; }
    if( Etype == CIO_UnMatch ) BSWAP32(dmy);
    if( dmy != type_dmy ) { fclose(fp); return false; }
  }
  //printf("***********************\n");
  //printf("*** read sph MxN ok ***\n");
  //printf("***********************\n");

  fclose(fp);
  //fclose(xfp);

  return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_Coarse(const char* fname, int step, int gc, int head[3], int tail[3], 
                           int G_Voxelsize[3], int R_head[3], int R_tail[3], int Voxelsize[3], 
                           int dfi_gc, int ncomp, int sta[3], int end[3], void* val, 
                           double &time, 
                           const bool mode, unsigned &step_avr, double &time_avr)
{

  int size_val = (Voxelsize[0]+2*dfi_gc)*(Voxelsize[1]+2*dfi_gc)*(Voxelsize[2]+2*dfi_gc)*ncomp;

  void *src;
  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
  if     ( Dtype == E_CIO_FLOAT32 ) src = new float[size_val]; 
  else if( Dtype == E_CIO_FLOAT64 ) src = new double[size_val];

  if( !read_S3D(fname, step, dfi_gc, src, time, mode, step_avr, time_avr) ) {
    delete [] src;
    val = NULL;
    printf("**** error read s3d ****");
    return false;
  }

  int idx,idx2,dg2,ii,jj,kk,i2,j2,k2;
  dg2 = (gc+1)/2*2;

  int szval[3],szdfi[3];
  for(int i=0; i<3; i++) szdfi[i] = Voxelsize[i];
  for(int i=0; i<3; i++) szval[i] = tail[i]-head[i]+1;

  if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || ncomp == 1 ) {
    //cio_dfi_coarse_ijkn_(val, head, tail, &gc, src, R_head, R_tail, &dfi_gc, &ncomp, sta, end, &m_RankID);
    for(int n=0; n<DFI_Finfo.Component; n++) {
    for(int k=sta[2]+1; k<=end[2]+1; k++) {
       kk = (k+dg2+1)/2-dg2/2 - R_head[2];
       k2 = k- head[2];
    for(int j=sta[1]+1; j<=end[1]+1; j++) {
       jj = (j+dg2+1)/2-dg2/2 - R_head[1];
       j2 = j- head[1];
    for(int i=sta[0]+1; i<=end[0]+1; i++) {
       ii = (i+dg2+1)/2-dg2/2 - R_head[0];
       i2 = i- head[0];

       idx  = _IDX_IJKN(i2,j2,k2,n,szval[0],szval[1],szval[2],gc);
       idx2 = _IDX_IJKN(ii,jj,kk,n,szdfi[0],szdfi[1],szdfi[2],DFI_Finfo.GuideCell);

       if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

    }}}}
  } else {
    //cio_dfi_coarse_nijk_(val, head, tail, &gc, src, R_head, R_tail, &dfi_gc, &ncomp, sta, end, &m_RankID);
    for(int k=sta[2]+1; k<=end[2]+1; k++) {
       kk = (k+dg2+1)/2-dg2/2 - R_head[2];
       k2 = k- head[2];
    for(int j=sta[1]+1; j<=end[1]+1; j++) {
       jj = (j+dg2+1)/2-dg2/2 - R_head[1];
       j2 = j- head[1];
    for(int i=sta[0]+1; i<=end[0]+1; i++) {
       ii = (i+dg2+1)/2-dg2/2 - R_head[0];
       i2 = i- head[0];
    for(int n=0; n<DFI_Finfo.Component; n++) {

       idx  = _IDX_NIJK(n,i2,j2,k2,DFI_Finfo.Component,szval[0],szval[1],szval[2],gc);
       idx2 = _IDX_NIJK(n,ii,jj,kk,DFI_Finfo.Component,szdfi[0],szdfi[1],szdfi[2],DFI_Finfo.GuideCell);

       if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

    }}}}
  }

  delete [] src;

  return true;

}

// #################################################################
// 
bool cio_DFI_SPH::read_Coarse_MxN(const char* fname, int step, int gc, int head[3], int tail[3], 
                           int G_Voxelsize[3], int R_head[3], int R_tail[3], int Voxelsize[3], 
                           int dfi_gc, int ncomp, int sta[3], int end[3], void *val, 
                           double &time, int n_dfi, 
                           const bool mode, unsigned &step_avr, double &time_avr)
{

  int read_s[3],read_e[3];
  int temp_head[3],temp_tail[3],t_voxsize[3];
  for(int i=0; i<3; i++ ) {
    temp_head[i]=(head[i]-1)/2+1;
    temp_tail[i]=(tail[i]-1)/2+1;
    if( temp_head[i] > R_head[i] ) temp_head[i] = R_head[i];
    if( temp_tail[i] < R_tail[i] ) temp_tail[i] = R_tail[i];
    t_voxsize[i]=(temp_tail[i]-temp_head[i])+1;
  }

  cio_EGlobalVoxel readflag = CIO_E_GV_SAME;
  CheckReadArea(temp_head,temp_tail,dfi_gc,R_head,R_tail,dfi_gc,readflag,read_s,read_e);

  int size_val = (t_voxsize[0]+2*dfi_gc)*(t_voxsize[1]+2*dfi_gc)*(t_voxsize[2]+2*dfi_gc)*ncomp;
  void *src;
  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
  if     ( Dtype == E_CIO_FLOAT32 ) src = new float[size_val]; 
  else if( Dtype == E_CIO_FLOAT64 ) src = new double[size_val];

  if( !read_MxN(fname,step,dfi_gc,temp_head,temp_tail,R_head,R_tail,read_s,read_e,src,
                time,n_dfi,mode,step_avr,time_avr) ) { 
    delete [] src;
    val = NULL;
    printf("**** error read Mxn ****");
    return false;
  }

  int idx,idx2,dg2,ii,jj,kk,i2,j2,k2;
  dg2 = (gc+1)/2*2;

  int szval[3],szdfi[3];
  //for(int i=0; i<3; i++) szdfi[i] = Voxelsize[i];
  for(int i=0; i<3; i++) szdfi[i] = temp_tail[i]-temp_head[i]+1;
  for(int i=0; i<3; i++) szval[i] = tail[i]-head[i]+1;

  if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || ncomp == 1 ) {
    //cio_dfi_coarse_ijkn_(val, head, tail, &gc, src, temp_head, temp_tail, &dfi_gc, &ncomp, sta, end, &m_RankID);

    for(int n=0; n<DFI_Finfo.Component; n++) {
    for(int k=sta[2]+1; k<=end[2]+1; k++) {
       kk = (k+dg2+1)/2-dg2/2 - temp_head[2];
       k2 = k- head[2];
    for(int j=sta[1]+1; j<=end[1]+1; j++) {
       jj = (j+dg2+1)/2-dg2/2 - temp_head[1];
       j2 = j- head[1];
    for(int i=sta[0]+1; i<=end[0]+1; i++) {
       ii = (i+dg2+1)/2-dg2/2 - temp_head[0];
       i2 = i- head[0];
       idx  = _IDX_IJKN(i2,j2,k2,n,szval[0],szval[1],szval[2],gc);
       //idx2 = _IDX_IJKN(ii,jj,kk,n,szdfi[0],szdfi[1],szdfi[2],DFI_Finfo.GuideCell);
       idx2 = _IDX_IJKN(ii,jj,kk,n,szdfi[0],szdfi[1],szdfi[2],dfi_gc);

       //int iunt=RankID+10;
       //dbidx_(&iunt,&ii,&jj,&kk,&i2,&j2,&k2,&idx,&idx2);

       if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

    }}}}
  } else {
    //cio_dfi_coarse_nijk_(val, head, tail, &gc, src, temp_head, temp_tail, &dfi_gc, &ncomp, sta, end, &m_RankID);
    for(int k=sta[2]+1; k<=end[2]+1; k++) {
       kk = (k+dg2+1)/2-dg2/2 - temp_head[2];
       k2 = k- head[2];
    for(int j=sta[1]+1; j<=end[1]+1; j++) {
       jj = (j+dg2+1)/2-dg2/2 - temp_head[1];
       j2 = j- head[1];
    for(int i=sta[0]+1; i<=end[0]+1; i++) {
       ii = (i+dg2+1)/2-dg2/2 - temp_head[0];
       i2 = i- head[0];
    for(int n=0; n<DFI_Finfo.Component; n++) {

       idx  = _IDX_NIJK(n,i2,j2,k2,DFI_Finfo.Component,szval[0],szval[1],szval[2],gc);
       //idx2 = _IDX_NIJK(n,ii,jj,kk,DFI_Finfo.Component,szdfi[0],szdfi[1],szdfi[2],DFI_Finfo.GuideCell);
       idx2 = _IDX_NIJK(n,ii,jj,kk,DFI_Finfo.Component,szdfi[0],szdfi[1],szdfi[2],dfi_gc);

       if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

    }}}}
  }

  delete [] src;

  return true;

}

// #################################################################
// 
void cio_DFI_SPH::WriteData(int step, int gc, void* time,
                           void *val, void *minmax, int interval,
                           const bool air_mode, const unsigned step_avr, const double time_avr,
                           bool force)
{

  if( !force  && interval>0 && step%interval != 0 ) return; 

  bool mio=false;
  if( DFI_MPI.NumberOfRank > 1 ) mio=true;

  //m_outSlice = true;

  std::string outFile = Generate_FileName(m_RankID,step,mio);

  if( MakeDirectory(DFI_Finfo.DirectoryPath,step) != 1 ) return;

  FILE* fp;
  if( (fp = fopen(outFile.c_str(),"wb")) == NULL ) {
    fprintf(stderr,"Can't open file.(%s)\n",outFile.c_str());
    return;
  }

  if( !write_head(fp, step, time, m_RankID) ) { fclose(fp); return; }

  if( gc == DFI_Finfo.GuideCell ) { 
    if( !write_data(fp, val, gc, m_RankID) ) { fclose(fp); return; }
  } else {
    int sz[3];
    for(int i=0; i<3; i++ ) sz[i] = RankInfo[m_RankID].VoxelSize[i]+2*DFI_Finfo.GuideCell;

    int size_val = sz[0]*sz[1]*sz[2]*DFI_Finfo.Component;
    void *src;
    E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
    if     ( Dtype == E_CIO_FLOAT32 ) src = new float[size_val]; 
    else if( Dtype == E_CIO_FLOAT64 ) src = new double[size_val];

    int gc2,idx,idx2;
    gc2 = min(DFI_Finfo.GuideCell,gc);
    if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || DFI_Finfo.Component == 1 ) {
      //cio_dfi_copy_ijkn_(src,sz,&DFI_Finfo.GuideCell,val,&gc,&DFI_Finfo.Component); 
      for(int n=0; n<DFI_Finfo.Component; n++) {
      for(int k=0; k<sz[2]+2*gc2; k++) {
      for(int j=0; j<sz[1]+2*gc2; j++) {
      for(int i=0; i<sz[0]+2*gc2; i++) {
         idx = _IDX_IJKN(i,j,k,n,sz[0],sz[1],sz[2],gc);
         idx2 = _IDX_IJKN(i,j,k,n,sz[0],sz[1],sz[2],DFI_Finfo.GuideCell);
         if( Dtype == E_CIO_FLOAT32 ) setval((float*)src,idx2,(float*)val,idx);
         if( Dtype == E_CIO_FLOAT64 ) setval((double*)src,idx2,(double*)val,idx);
      }}}}
    } else {
      //cio_dfi_copy_nijk_(src,sz,&DFI_Finfo.GuideCell,val,&gc,&DFI_Finfo.Component); 
      for(int k=0; k<sz[2]+2*gc2; k++) {
      for(int j=0; j<sz[1]+2*gc2; j++) {
      for(int i=0; i<sz[0]+2*gc2; i++) {
      for(int n=0; n<DFI_Finfo.Component; n++) {
         idx = _IDX_NIJK(n,i,j,k,DFI_Finfo.Component,sz[0],sz[1],sz[2],gc);
         idx2 = _IDX_NIJK(n,i,j,k,DFI_Finfo.Component,sz[0],sz[1],sz[2],DFI_Finfo.GuideCell);
         if( Dtype == E_CIO_FLOAT32 ) setval((float*)src,idx2,(float*)val,idx);
         if( Dtype == E_CIO_FLOAT64 ) setval((double*)src,idx2,(double*)val,idx);
      }}}}
    }
    if( !write_data(fp, src, DFI_Finfo.GuideCell, m_RankID) ) { 
      fclose(fp); 
      delete [] src;
      return;
    }

    delete [] src;
  }

//average
  if( !air_mode ) {
    int dType = 0;
    if( !strcasecmp(DFI_Finfo.DataType.c_str(),"Float32")) dType = 1;
    if( !strcasecmp(DFI_Finfo.DataType.c_str(),"Float64")) dType = 2;

    unsigned int dmy;
    int Int_size,Real_size;
    if ( dType == 1 ) {
      dmy = 8;
      Int_size  = sizeof(int);
      Real_size = sizeof(float);
    }else{
      dmy = 16;
      Int_size  = sizeof(long long);
      Real_size = sizeof(double);
    }
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return ; }
    if( dType == 1 ){
      int istep = (int)step_avr;
      float ttime = (float)time_avr;
      if( fwrite(&istep, Int_size, 1, fp) != 1 ) { fclose(fp); return ; }
      if( fwrite(&ttime, Real_size, 1, fp) != 1 ) { fclose(fp); return ; }
    } else {
      long long dstep = (long long)step_avr;
      double ttime = (double)time_avr;
      if( fwrite(&dstep, Int_size, 1, fp) != 1 ) { fclose(fp); return ; }
      if( fwrite(&ttime, Real_size, 1, fp) != 1 ) { fclose(fp); return ; }
    }
    if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) { fclose(fp); return ; }
  }

  fclose(fp);

  WriteIndexDfiFile(m_indexDfiName,m_RankID,DFI_Finfo.Prefix,step,time,minmax,true,
                    air_mode, step_avr, time_avr);
}

// #################################################################
// 
bool cio_DFI_SPH::write_head(FILE* fp, int step, void* time, int n)
{

  int svType = 0;
  if( DFI_Finfo.Component == 1 ) svType = 1;
  if( DFI_Finfo.Component > 1  ) svType = 2;
  if( svType == 0 ) return false;

  int dType = 0;
  if( !strcasecmp(DFI_Finfo.DataType.c_str(),"Float32")) dType = 1; 
  if( !strcasecmp(DFI_Finfo.DataType.c_str(),"Float64")) dType = 2; 
  if( dType == 0 ) return false;

  unsigned int dmy;
  dmy = 8;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&svType, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dType, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;


  if( dType == 1 ) dmy = 12; //float
  else             dmy = 24; //double

  //voxel size
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( dType == 1 ) {
    int size[3];
    for(int i=0; i<3; i++ ) size[i] = (int)RankInfo[n].VoxelSize[i]+(int)(2*DFI_Finfo.GuideCell);
    if( fwrite(size, sizeof(int), 3, fp) !=3 ) return false;
  } else {
    long long  size[3];
    for(int i=0; i<3; i++ ) size[i] = (long long)RankInfo[n].VoxelSize[i]+(long long)(2*DFI_Finfo.GuideCell);
    if( fwrite(size, sizeof(long long), 3, fp) !=3 ) return false;
  }

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  //origin
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( dType == 1 ) {
    float org[3];
    for(int i=0; i<3; i++ ) org[i]=(float)DFI_Domain.GlobalOrigin[i];
    if( fwrite(org, sizeof(float), 3, fp) !=3 ) return false;
  } else {
    double org[3];
    for(int i=0; i<3; i++ ) org[i]=(double)DFI_Domain.GlobalOrigin[i];
    if( fwrite(org, sizeof(double), 3, fp) !=3 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  //pitch
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( dType == 1 ) {
    float pch[3];
    for(int i=0; i<3; i++ ) pch[i]=(float)DFI_Domain.GlobalRegion[i]/DFI_Domain.GlobalVoxel[i];
    if( fwrite(pch, sizeof(float), 3, fp) !=3 ) return false;
  } else {
    double pch[3];
    for(int i=0; i<3; i++ ) pch[i]=(double)DFI_Domain.GlobalRegion[i]/DFI_Domain.GlobalVoxel[i];
    if( fwrite(pch, sizeof(double), 3, fp) !=3 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;

  //step&time
  int Int_size,Real_size;
  if ( dType == 1 ) {
    dmy = 8;
    Int_size  = sizeof(int);
    Real_size = sizeof(float);
  }else{
    dmy = 16;
    Int_size  = sizeof(long long);
    Real_size = sizeof(double);
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( dType == 1 ){
    float *ttime = (float *)time;
    if( fwrite(&step, Int_size, 1, fp) != 1 ) return false;
    if( fwrite(ttime, Real_size, 1, fp) != 1 ) return false;
  } else {
    long long dstep = (long long)step;
    double *ttime = (double *)time;
    if( fwrite(&dstep, Int_size, 1, fp) != 1 ) return false;
    if( fwrite(ttime, Real_size, 1, fp) != 1 ) return false;
  }
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  
  return true;
}

// #################################################################
// 
bool cio_DFI_SPH::write_data(FILE* fp, void* val, int gc, int n)
{

  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
  int Real_size = get_cio_Datasize(Dtype);

  int size[3];
  for(int i=0; i<3; i++ ) size[i] = (int)RankInfo[n].VoxelSize[i]+(int)(2*gc);

  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.Component > 1 ) dLen *= 3;  

  unsigned int dmy = dLen * Real_size;

  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  if( fwrite(val, Real_size, dLen, fp) != dLen ) return false; 
  if( fwrite(&dmy, sizeof(int), 1, fp) != 1 ) return false;
  return true;
}
