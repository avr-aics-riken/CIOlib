/*
 * CIOlib - Cartesian Input / Output library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

/** 
 * @file   cio_DFI_BOV.C
 * @brief  cio_DFI_BOV Class
 * @author kero    
 */

#include "cio_DFI.h"
#include "cio_DFI_BOV.h"

// #################################################################
// コンストラクタ
cio_DFI_BOV::cio_DFI_BOV()
{

}


// #################################################################
// デストラクタ
cio_DFI_BOV::~cio_DFI_BOV()
{

}


// #################################################################
//
void* cio_DFI_BOV::ReadData(int step, int gc,
                           int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                           double &time,
                           const bool mode, unsigned &idummy, double &fdummy)
{

  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);

  int sz[3];
  for(int i=0; i<3; i++) sz[i]=tail[i]-head[i]+1+2*gc;
  int val_size = sz[0]*sz[1]*sz[2];
  //printf("ID : %d val_size %d\n",RankID,val_size);

  void* val;

  if       ( Dtype == E_CIO_INT8    ) {                           //int8
      val = new char[val_size];
      memset(val, 0, sizeof(char)*val_size);
  }else  if( Dtype == E_CIO_INT16   ) {                           //int16
      val = new short[val_size];
      memset(val, 0, sizeof(short)*val_size);
  }else  if( Dtype == E_CIO_INT32   ) {                           //int32
      val = new int[val_size];
      memset(val, 0, sizeof(int)*val_size);
  }else  if( Dtype == E_CIO_INT64   ) {                           //int64
      val = new long long[val_size];
      memset(val, 0, sizeof(long long)*val_size);
  }else  if( Dtype == E_CIO_UINT8   ) {                           //unit8
      val = new unsigned char[val_size];
      memset(val, 0, sizeof(unsigned char)*val_size);
  }else  if( Dtype == E_CIO_UINT16  ) {                           //unit16
      val = new unsigned short[val_size];
      memset(val, 0, sizeof(unsigned short)*val_size);
  }else  if( Dtype == E_CIO_UINT32  ) {                           //uint32
      val = new unsigned int[val_size];
      memset(val, 0, sizeof(unsigned int)*val_size);
  }else  if( Dtype == E_CIO_UINT64  ) {                           //uint64
      val = new unsigned long long[val_size];
      memset(val, 0, sizeof(unsigned long long)*val_size);
  }else  if( Dtype == E_CIO_FLOAT32 ) {                           //float32
      val = new float[val_size];
      memset(val, 0, sizeof(float)*val_size);
  } else if( Dtype == E_CIO_FLOAT64 ) {                           //float64
      val = new double[val_size];
      memset(val, 0, sizeof(double)*val_size);
  }

  ReadData(step, gc, Gvoxel, Gdivision, head, tail, val, time, mode, idummy, fdummy);
  return val;
}

// #################################################################
//
void cio_DFI_BOV::ReadData(int step, int gc,
                           int Gvoxel[3], int Gdivision[3], int head[3], int tail[3],
                           void *val, double &time,
                           const bool mode, unsigned &idummy, double &fdummy)
{

  //printf("###### BOV READ START #######\n");

  bool mio = false;

  if( RankInfo.size() > 1 ) mio = true;

/*
  cio_Domain G_Domain;
  vector<cio_Rank> G_RankInfo;
  cio_Rank G_Rank;
  cio_Create_Domain(m_comm,Gvoxel, Gdivision, head, tail,G_Domain,G_RankInfo,G_Rank);
*/

  cio_EGlobalVoxel readflag; //CIO_E_GV_SAME:密　CIO_E_GVX2_SAME:粗 CIO_E_OTHER:その他
  readflag = CheckGlobalVoxel(Gvoxel, DFI_Domain.GlobalVoxel);

  vector<int> RList;
  CreateRankList(head, tail, gc, readflag, RList);

  if( readflag == CIO_E_GV_SAME ) {
    for(int i=0; i<RList.size(); i++) {
      int n = RList[i];
      int ID= RankInfo[n].RankID;
      std::string fname = Generate_FileName(ID,step,mio);

      int sta[3],end[3];

      int voxsize[3];
      for(int i=0; i<3; i++) voxsize[i]=RankInfo[n].VoxelSize[i]+(2*DFI_Finfo.GuideCell);

      if( CheckReadArea(head, tail, gc, RankInfo[n].HeadIndex, RankInfo[n].TailIndex,
                        DFI_Finfo.GuideCell, readflag, sta, end) )
      {
        if( !read_S3D(fname.c_str(), step, DFI_Finfo.GuideCell, val, voxsize, time ) ) {
          val = NULL;
          printf("**** error read s3d ****");
          return;
        }

      } else {
        if( !read_MxN(fname.c_str(), step, gc, head, tail, RankInfo[n].HeadIndex,
                    RankInfo[n].TailIndex, sta, end, val, time, n, voxsize ) ) {
          val = NULL;
          printf("**** error read mxn ID %d dfi_ID : %d fname : %s\n",m_RankID,n,
                  fname.c_str());
          return;
        }
      }
      if( m_RankID == 0 ) {
        printf("\t[%s] has read :\tstep=%d ]\n",fname.c_str(), step);
      }

    }
  } else if( readflag == CIO_E_GVX2_SAME ) {

    for(int i=0; i<RList.size(); i++) {
      std::string fname = Generate_FileName(RList[i],step,mio);
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
                          DFI_Finfo.Component, sta, end, val, time) ) {
           val = NULL;
           return ;
        }
      } else {
      //1対多の読込み&補間
        if( !read_Coarse_MxN(fname.c_str(), step, gc, head, tail, Gvoxel, RankInfo[n].HeadIndex,
                          RankInfo[n].TailIndex, RankInfo[n].VoxelSize, DFI_Finfo.GuideCell,
                          DFI_Finfo.Component, sta, end, val, time, n) ) {
           val = NULL;
           return ;
        }
      }
      if( m_RankID == 0 ) {
        printf("\t[%s] has read :\tstep=%d ]\n",fname.c_str(), step);
      }
    }
  } else {
   printf("Error filed data format\n");
   return;
  }

  time=0.0;
  for(int i=0; i<TimeSlice.size(); i++) {
     if( TimeSlice[i].step == step ) {
       time=(double)TimeSlice[i].time;
       return;
     }
  }

}

// #################################################################
//
bool cio_DFI_BOV::read_S3D(const char* fname, int step, int gc, void*val, int voxsize[3],
                           double &time)
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

  unsigned int dmy;
//data
  int arraylen = voxsize[0]*voxsize[1]*voxsize[2]*DFI_Finfo.Component;
  if( fread(val, Dsize, arraylen, fp) != arraylen ) { fclose(fp); return false; }
  if( Etype == CIO_UnMatch ) {
    if(        Dtype == E_CIO_INT8 ) {
      char *val2 = (char *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_INT16 ) {
      short *val2 = (short *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_INT32 ) {
      int *val2 = (int *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_INT64 ) {
      long long *val2 = (long long *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_UINT8 ) {
      unsigned char *val2 = (unsigned char *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_UINT16 ) {
      unsigned short *val2 = (unsigned short *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_UINT32 ) {
      unsigned int *val2 = (unsigned int *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_UINT64 ) {
      unsigned long long *val2 = (unsigned long long *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_FLOAT32 ) {
      float *val2 = (float *)val;
      BSWAPVEC(val2,arraylen);
    } else if( Dtype == E_CIO_FLOAT64 ) {
      double *val2 = (double *)val;
      DBSWAPVEC(val2,arraylen);
    }
  }

  fclose(fp);

  //printf("**** BOV read s3d ok!!!!\n");

  return true;

}

// #################################################################
//
bool cio_DFI_BOV::read_MxN(const char* fname, int step, int gc, int head[3], int tail[3],
                           int R_head[3], int R_tail[3], int sta[3], int end[3], void* val,
                           double &time, int n_dfi, int voxsize[3])
{

  //printf("###### BOV MxN READ START #######\n");

  if( !fname || !DFI_Finfo.Component ) {
    printf("*** error fname %s DFI_Finfo.Component %d\n",fname,
            DFI_Finfo.Component);
    return false;
  }

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

  if( Etype == CIO_UnKnown ) {
    printf("*** error Endian type unknown\n");
    return false;
  }

  FILE* fp;
  if( !(fp=fopen(fname,"rb")) ) {
    printf("Can't open file. (%s)\n",fname);
    //return false;
    exit(0);
  }

  //printf("fname : %s\n",fname);

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

  long long idx,idx2;
  int imax,jmax,kmax;

  imax = tail[0]-head[0]+1;
  jmax = tail[1]-head[1]+1;
  kmax = tail[2]-head[2]+1;

  if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || DFI_Finfo.Component == 1 ) {

    int arraylen = (voxsize[0])*(voxsize[1]);
    void *rbuff;
    if(        Dtype == E_CIO_INT8 ) {
      rbuff = new char[arraylen];
    } else if( Dtype == E_CIO_INT16 ) {
      rbuff = new short[arraylen];
    } else if( Dtype == E_CIO_INT32 ) {
      rbuff = new int[arraylen];
    } else if( Dtype == E_CIO_INT64 ) {
      rbuff = new long long[arraylen];
    } else if( Dtype == E_CIO_UINT8 ) {
      rbuff = new unsigned char[arraylen];
    } else if( Dtype == E_CIO_UINT16 ) {
      rbuff = new unsigned short[arraylen];
    } else if( Dtype == E_CIO_UINT32 ) {
      rbuff = new unsigned int[arraylen];
    } else if( Dtype == E_CIO_UINT64 ) {
      rbuff = new unsigned long long[arraylen];
    } else if( Dtype == E_CIO_FLOAT32 ) {
      rbuff = new float[arraylen];
    } else if( Dtype == E_CIO_FLOAT64 ) {
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
           if( fread(rbuff, Dsize, arraylen, fp) != arraylen ) { fclose(fp); return false;}
           if( Etype == CIO_UnMatch ) {
             if(        Dtype == E_CIO_INT8 ) {
               char *rbuff2 = (char *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_INT16 ) {
               short *rbuff2 = (short *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_INT32 ) {
               int *rbuff2 = (int *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_INT64 ) {
               long long *rbuff2 = (long long *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_UINT8 ) {
               unsigned char *rbuff2 = (unsigned char *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_UINT16 ) {
               unsigned short *rbuff2 = (unsigned short *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_UINT32 ) {
               unsigned int *rbuff2 = (unsigned int *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_UINT64 ) {
               unsigned long long *rbuff2 = (unsigned long long *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_FLOAT32 ) {
               float *rbuff2 = (float *)rbuff;
               BSWAPVEC(rbuff2,arraylen);
             } else if( Dtype == E_CIO_FLOAT64 ) {
               double *rbuff2 = (double *)rbuff;
               DBSWAPVEC(rbuff2,arraylen);
             }
           }
        }
        nkrec-=DFI_Finfo.GuideCell;
      }
      for(int n=0; n<nkrec; n++) {
        if( fread(rbuff, Dsize, arraylen, fp) != arraylen ) {fclose(fp); return false;}
        if( Etype == CIO_UnMatch ) {
          if(        Dtype == E_CIO_INT8 ) {
            char *rbuff2 = (char *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_INT16 ) {
            short *rbuff2 = (short *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_INT32 ) {
            int *rbuff2 = (int *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_INT64 ) {
            long long *rbuff2 = (long long *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_UINT8 ) {
            unsigned char *rbuff2 = (unsigned char *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_UINT16 ) {
            unsigned short *rbuff2 = (unsigned short *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_UINT32 ) {
            unsigned int *rbuff2 = (unsigned int *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_UINT64 ) {
            unsigned long long *rbuff2 = (unsigned long long *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_FLOAT32 ) {
            float *rbuff2 = (float *)rbuff;
            BSWAPVEC(rbuff2,arraylen);
          } else if( Dtype == E_CIO_FLOAT64 ) {
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
           if( Dtype == E_CIO_INT8    ) setval((char *)val,idx,(char *)rbuff,idx2);
           if( Dtype == E_CIO_INT16   ) setval((short *)val,idx,(short *)rbuff,idx2);
           if( Dtype == E_CIO_INT32   ) setval((int *)val,idx,(int *)rbuff,idx2);
           if( Dtype == E_CIO_INT64   ) setval((long long *)val,idx,(long long *)rbuff,idx2);
           if( Dtype == E_CIO_UINT8   ) setval((unsigned char *)val,idx,(unsigned char *)rbuff,idx2);
           if( Dtype == E_CIO_UINT16  ) setval((unsigned short *)val,idx,(unsigned short *)rbuff,idx2);
           if( Dtype == E_CIO_UINT32  ) setval((unsigned int *)val,idx,(unsigned int *)rbuff,idx2);
           if( Dtype == E_CIO_UINT64  ) setval((unsigned long long *)val,idx,(unsigned long long *)rbuff,idx2);
           if( Dtype == E_CIO_FLOAT32 ) setval((float *)val,idx,(float *)rbuff,idx2);
           if( Dtype == E_CIO_FLOAT64 ) setval((double *)val,idx,(double *)rbuff,idx2);
        }}
      }
    }

    delete [] rbuff;

//nijk
  } else if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "nijk" ) ) {

    int arraylen = (voxsize[0])*(voxsize[1])*DFI_Finfo.Component;
    void *rbuff;
    if(        Dtype == E_CIO_INT8 ) {
      rbuff = new char[arraylen];
    } else if( Dtype == E_CIO_INT16 ) {
      rbuff = new short[arraylen];
    } else if( Dtype == E_CIO_INT32 ) {
      rbuff = new int[arraylen];
    } else if( Dtype == E_CIO_INT64 ) {
      rbuff = new long long[arraylen];
    } else if( Dtype == E_CIO_UINT8 ) {
      rbuff = new unsigned char[arraylen];
    } else if( Dtype == E_CIO_UINT16 ) {
      rbuff = new unsigned short[arraylen];
    } else if( Dtype == E_CIO_UINT32 ) {
      rbuff = new unsigned int[arraylen];
    } else if( Dtype == E_CIO_UINT64 ) {
      rbuff = new unsigned long long[arraylen];
    } else if( Dtype == E_CIO_FLOAT32 ) {
      rbuff = new float[arraylen];
    } else if( Dtype == E_CIO_FLOAT64 ) {
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
           if(        Dtype == E_CIO_INT8 ) {
             char* rbuff2 = (char *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_INT16 ) {
             short* rbuff2 = (short *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_INT32 ) {
             int* rbuff2 = (int *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_INT64 ) {
             long long* rbuff2 = (long long *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_UINT8 ) {
             unsigned char* rbuff2 = (unsigned char *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_UINT16 ) {
             unsigned short* rbuff2 = (unsigned short *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_UINT32 ) {
             unsigned int* rbuff2 = (unsigned int *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_UINT64 ) {
             unsigned long long* rbuff2 = (unsigned long long *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_FLOAT32 ) {
             float* rbuff2 = (float *)rbuff;
             BSWAPVEC(rbuff2,arraylen);
           } else if( Dtype == E_CIO_FLOAT64 ) {
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
        if(        Dtype == E_CIO_INT8 ) {
          char* rbuff2 = (char *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_INT16 ) {
          short* rbuff2 = (short *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_INT32 ) {
          int* rbuff2 = (int *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_INT64 ) {
          long long* rbuff2 = (long long *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_UINT8 ) {
          unsigned char* rbuff2 = (unsigned char *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_UINT16 ) {
          unsigned short* rbuff2 = (unsigned short *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_UINT32 ) {
          unsigned int* rbuff2 = (unsigned int *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_UINT64 ) {
          unsigned long long* rbuff2 = (unsigned long long *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_FLOAT32 ) {
          float* rbuff2 = (float *)rbuff;
          BSWAPVEC(rbuff2,arraylen);
        } else if( Dtype == E_CIO_FLOAT64 ) {
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
         if( Dtype == E_CIO_INT8    ) setval((char *)val,idx,(char *)rbuff,idx2);
         if( Dtype == E_CIO_INT16   ) setval((short *)val,idx,(short *)rbuff,idx2);
         if( Dtype == E_CIO_INT32   ) setval((int *)val,idx,(int *)rbuff,idx2);
         if( Dtype == E_CIO_INT64   ) setval((long long *)val,idx,(long long *)rbuff,idx2);
         if( Dtype == E_CIO_UINT8   ) setval((unsigned char *)val,idx,(unsigned char *)rbuff,idx2);
         if( Dtype == E_CIO_UINT16  ) setval((unsigned short *)val,idx,(unsigned short *)rbuff,idx2);
         if( Dtype == E_CIO_UINT32  ) setval((unsigned int *)val,idx,(unsigned int *)rbuff,idx2);
         if( Dtype == E_CIO_UINT64  ) setval((unsigned long long *)val,idx,(unsigned long long *)rbuff,idx2);
         if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)rbuff,idx2);
         if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)rbuff,idx2);
      }}}
    }

    delete [] rbuff;

  }

  fclose(fp);

  return true;
}
// #################################################################
//
bool cio_DFI_BOV::read_Coarse(const char* fname, int step, int gc, int head[3], int tail[3],
                           int G_Voxelsize[3], int R_head[3], int R_tail[3], int Voxelsize[3],
                           int dfi_gc, int ncomp, int sta[3], int end[3], void* val,
                           double &time)
{

  //printf("**** BOV read_Coarse in *****\n");

  int size_val = (Voxelsize[0]+2*dfi_gc)*(Voxelsize[1]+2*dfi_gc)*(Voxelsize[2]+2*dfi_gc)*ncomp;

  void *src;
  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
  if     ( Dtype == E_CIO_INT8    ) src = new char[size_val];
  else if( Dtype == E_CIO_INT16   ) src = new short[size_val];
  else if( Dtype == E_CIO_INT32   ) src = new int[size_val];
  else if( Dtype == E_CIO_INT64   ) src = new long long[size_val];
  else if( Dtype == E_CIO_UINT8   ) src = new unsigned char[size_val];
  else if( Dtype == E_CIO_UINT16  ) src = new unsigned short[size_val];
  else if( Dtype == E_CIO_UINT32  ) src = new unsigned int[size_val];
  else if( Dtype == E_CIO_UINT64  ) src = new unsigned long long[size_val];
  else if( Dtype == E_CIO_FLOAT32 ) src = new float[size_val];
  else if( Dtype == E_CIO_FLOAT64 ) src = new double[size_val];

  int vsize[3];
  for( int i=0; i<3; i++ ) vsize[i]=Voxelsize[i]+2*dfi_gc;

  if( !read_S3D(fname, step, dfi_gc, src, vsize, time) ) {
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

       if     ( Dtype == E_CIO_INT8    ) setval((char*)val,idx,(char*)src,idx2);
       else if( Dtype == E_CIO_INT16   ) setval((short*)val,idx,(short*)src,idx2);
       else if( Dtype == E_CIO_INT32   ) setval((int*)val,idx,(int*)src,idx2);
       else if( Dtype == E_CIO_INT64   ) setval((long long*)val,idx,(long long*)src,idx2);
       else if( Dtype == E_CIO_UINT8   ) setval((unsigned char*)val,idx,(unsigned char*)src,idx2);
       else if( Dtype == E_CIO_UINT16  ) setval((unsigned short*)val,idx,(unsigned short*)src,idx2);
       else if( Dtype == E_CIO_UINT32  ) setval((unsigned int*)val,idx,(unsigned int*)src,idx2);
       else if( Dtype == E_CIO_UINT64  ) setval((unsigned long long*)val,idx,(unsigned long long*)src,idx2);
       else if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       else if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

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

       if     ( Dtype == E_CIO_INT8    ) setval((char*)val,idx,(char*)src,idx2);
       else if( Dtype == E_CIO_INT16   ) setval((short*)val,idx,(short*)src,idx2);
       else if( Dtype == E_CIO_INT32   ) setval((int*)val,idx,(int*)src,idx2);
       else if( Dtype == E_CIO_INT64   ) setval((long long*)val,idx,(long long*)src,idx2);
       else if( Dtype == E_CIO_UINT8   ) setval((unsigned char*)val,idx,(unsigned char*)src,idx2);
       else if( Dtype == E_CIO_UINT16  ) setval((unsigned short*)val,idx,(unsigned short*)src,idx2);
       else if( Dtype == E_CIO_UINT32  ) setval((unsigned int*)val,idx,(unsigned int*)src,idx2);
       else if( Dtype == E_CIO_UINT64  ) setval((unsigned long long*)val,idx,(unsigned long long*)src,idx2);
       else if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       else if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

    }}}}
  }
  delete [] src;

  return true;
}

// #################################################################
//
bool cio_DFI_BOV::read_Coarse_MxN(const char* fname, int step, int gc, int head[3], int tail[3],
                           int G_Voxelsize[3], int R_head[3], int R_tail[3], int Voxelsize[3],
                           int dfi_gc, int ncomp, int sta[3], int end[3], void *val,
                           double &time, int n_dfi)
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
  if     ( Dtype == E_CIO_INT8    ) src = new char[size_val];
  else if( Dtype == E_CIO_INT16   ) src = new short[size_val];
  else if( Dtype == E_CIO_INT32   ) src = new int[size_val];
  else if( Dtype == E_CIO_INT64   ) src = new long long[size_val];
  else if( Dtype == E_CIO_UINT8   ) src = new unsigned char[size_val];
  else if( Dtype == E_CIO_UINT16  ) src = new unsigned short[size_val];
  else if( Dtype == E_CIO_UINT32  ) src = new unsigned int[size_val];
  else if( Dtype == E_CIO_UINT64  ) src = new unsigned long long[size_val];
  else if( Dtype == E_CIO_FLOAT32 ) src = new float[size_val];
  else if( Dtype == E_CIO_FLOAT64 ) src = new double[size_val];

  int vsize[3];
  //for( int i=0; i<3; i++ ) vsize[i]=t_voxsize[i]+2*dfi_gc;
  for( int i=0; i<3; i++ ) vsize[i]=(R_tail[i]-R_head[i]+1)+2*dfi_gc;

  if( !read_MxN(fname,step,dfi_gc,temp_head,temp_tail,R_head,R_tail,read_s,read_e,src,time,n_dfi,
                vsize ) ) {
    delete [] src;
    val = NULL;
    printf("**** error read Mxn ****\n");
    return false;
  }

  int idx,idx2,dg2,ii,jj,kk,i2,j2,k2;
  dg2 = (gc+1)/2*2;

  int szval[3],szdfi[3];
  //for(int i=0; i<3; i++) szdfi[i] = Voxelsize[i];
  for(int i=0; i<3; i++) szdfi[i] = Voxelsize[i];
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

       if     ( Dtype == E_CIO_INT8    ) setval((char*)val,idx,(char*)src,idx2);
       else if( Dtype == E_CIO_INT16   ) setval((short*)val,idx,(short*)src,idx2);
       else if( Dtype == E_CIO_INT32   ) setval((int*)val,idx,(int*)src,idx2);
       else if( Dtype == E_CIO_INT64   ) setval((long long*)val,idx,(long long*)src,idx2);
       else if( Dtype == E_CIO_UINT8   ) setval((unsigned char*)val,idx,(unsigned char*)src,idx2);
       else if( Dtype == E_CIO_UINT16  ) setval((unsigned short*)val,idx,(unsigned short*)src,idx2);
       else if( Dtype == E_CIO_UINT32  ) setval((unsigned int*)val,idx,(unsigned int*)src,idx2);
       else if( Dtype == E_CIO_UINT64  ) setval((unsigned long long*)val,idx,(unsigned long long*)src,idx2);
       else if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       else if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

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
       idx2 = _IDX_NIJK(n,ii,jj,kk,DFI_Finfo.Component,szdfi[0],szdfi[1],szdfi[2],dfi_gc);

       if     ( Dtype == E_CIO_INT8    ) setval((char*)val,idx,(char*)src,idx2);
       else if( Dtype == E_CIO_INT16   ) setval((short*)val,idx,(short*)src,idx2);
       else if( Dtype == E_CIO_INT32   ) setval((int*)val,idx,(int*)src,idx2);
       else if( Dtype == E_CIO_INT64   ) setval((long long*)val,idx,(long long*)src,idx2);
       else if( Dtype == E_CIO_UINT8   ) setval((unsigned char*)val,idx,(unsigned char*)src,idx2);
       else if( Dtype == E_CIO_UINT16  ) setval((unsigned short*)val,idx,(unsigned short*)src,idx2);
       else if( Dtype == E_CIO_UINT32  ) setval((unsigned int*)val,idx,(unsigned int*)src,idx2);
       else if( Dtype == E_CIO_UINT64  ) setval((unsigned long long*)val,idx,(unsigned long long*)src,idx2);
       else if( Dtype == E_CIO_FLOAT32 ) setval((float*)val,idx,(float*)src,idx2);
       else if( Dtype == E_CIO_FLOAT64 ) setval((double*)val,idx,(double*)src,idx2);

    }}}}
  }

  delete [] src;

  return true;
}

// #################################################################
//
void cio_DFI_BOV::WriteData(int step, int gc, void* time,
                           void *val, void *minmax, int interval, 
                           const bool air_mode, const unsigned step_avr, const double time_avr,
                           bool force)
{

  if( !force  && interval>0 && step%interval != 0 ) return;

  //printf("###### BOV WRITE START #######\n");

  bool mio=false;
  if( DFI_MPI.NumberOfRank > 1 ) mio=true;

  std::string outFile = Generate_FileName(m_RankID,step,mio);
  if( MakeDirectory(DFI_Finfo.DirectoryPath) != 1 ) return;

  FILE* fp;
  if( (fp = fopen(outFile.c_str(),"wb")) == NULL ) {
    fprintf(stderr,"Can't open file.(%s)\n",outFile.c_str());
    return;
  }

  if( gc == DFI_Finfo.GuideCell ) {
    if( !write_data(fp, val, gc, m_RankID) ) { fclose(fp); return; }
  } else {

    int gc2,idx,idx2;
    gc2 = min(DFI_Finfo.GuideCell,gc);

    int sz[3];
    for(int i=0; i<3; i++ ) sz[i] = RankInfo[m_RankID].VoxelSize[i]+2*DFI_Finfo.GuideCell;

    int size_val = sz[0]*sz[1]*sz[2]*DFI_Finfo.Component;
    void *src;
    E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
    if     ( Dtype == E_CIO_INT8    ) src = new char[size_val];
    else if( Dtype == E_CIO_INT16   ) src = new short[size_val];
    else if( Dtype == E_CIO_INT32   ) src = new int[size_val];
    else if( Dtype == E_CIO_INT64   ) src = new long long[size_val];
    else if( Dtype == E_CIO_UINT8   ) src = new unsigned char[size_val];
    else if( Dtype == E_CIO_UINT16  ) src = new unsigned short[size_val];
    else if( Dtype == E_CIO_UINT32  ) src = new unsigned int[size_val];
    else if( Dtype == E_CIO_UINT64  ) src = new unsigned long long[size_val];
    else if( Dtype == E_CIO_FLOAT32 ) src = new float[size_val];
    else if( Dtype == E_CIO_FLOAT64 ) src = new double[size_val];
    else return;

    if( !strcasecmp(DFI_Finfo.ArrayShape.c_str(), "ijkn") || DFI_Finfo.Component == 1 ) {
      //cio_dfi_copy_ijkn_(src,sz,&DFI_Finfo.GuideCell,val,&gc,&DFI_Finfo.Component);
      for(int n=0; n<DFI_Finfo.Component; n++) {
      for(int k=0; k<sz[2]+2*gc2; k++) {
      for(int j=0; j<sz[1]+2*gc2; j++) {
      for(int i=0; i<sz[0]+2*gc2; i++) {
         idx = _IDX_IJKN(i,j,k,n,sz[0],sz[1],sz[2],gc);
         idx2 = _IDX_IJKN(i,j,k,n,sz[0],sz[1],sz[2],DFI_Finfo.GuideCell);
         if     ( Dtype == E_CIO_INT8    ) setval((char*)src,               idx2,(char*)val,idx);
         else if( Dtype == E_CIO_INT16   ) setval((short*)src,              idx2,(short*)val,idx);
         else if( Dtype == E_CIO_INT32   ) setval((int*)src,                idx2,(int*)val,idx);
         else if( Dtype == E_CIO_INT64   ) setval((long long *)src,         idx2,(long long*)val,idx);
         else if( Dtype == E_CIO_UINT8   ) setval((unsigned char *)src,     idx2,(unsigned char*)val,idx);
         else if( Dtype == E_CIO_UINT16  ) setval((unsigned short *)src,    idx2,(unsigned short*)val,idx);
         else if( Dtype == E_CIO_UINT32  ) setval((unsigned int *)src,      idx2,(unsigned int*)val,idx);
         else if( Dtype == E_CIO_UINT64  ) setval((unsigned long long *)src,idx2,(unsigned long long*)val,idx);
         else if( Dtype == E_CIO_FLOAT32 ) setval((float*)src,              idx2,(float*)val,idx);
         else if( Dtype == E_CIO_FLOAT64 ) setval((double*)src,             idx2,(double*)val,idx);
      }}}}
    } else {
      //cio_dfi_copy_nijk_(src,sz,&DFI_Finfo.GuideCell,val,&gc,&DFI_Finfo.Component);
      for(int k=0; k<sz[2]+2*gc2; k++) {
      for(int j=0; j<sz[1]+2*gc2; j++) {
      for(int i=0; i<sz[0]+2*gc2; i++) {
      for(int n=0; n<DFI_Finfo.Component; n++) {
         idx = _IDX_NIJK(n,i,j,k,DFI_Finfo.Component,sz[0],sz[1],sz[2],gc);
         idx2 = _IDX_NIJK(n,i,j,k,DFI_Finfo.Component,sz[0],sz[1],sz[2],DFI_Finfo.GuideCell);
         if     ( Dtype == E_CIO_INT8    ) setval((char*)src,               idx2,(char*)val,idx);
         else if( Dtype == E_CIO_INT16   ) setval((short*)src,              idx2,(short*)val,idx);
         else if( Dtype == E_CIO_INT32   ) setval((int*)src,                idx2,(int*)val,idx);
         else if( Dtype == E_CIO_INT64   ) setval((long long *)src,         idx2,(long long*)val,idx);
         else if( Dtype == E_CIO_UINT8   ) setval((unsigned char *)src,     idx2,(unsigned char*)val,idx);
         else if( Dtype == E_CIO_UINT16  ) setval((unsigned short *)src,    idx2,(unsigned short*)val,idx);
         else if( Dtype == E_CIO_UINT32  ) setval((unsigned int *)src,      idx2,(unsigned int*)val,idx);
         else if( Dtype == E_CIO_UINT64  ) setval((unsigned long long *)src,idx2,(unsigned long long*)val,idx);
         else if( Dtype == E_CIO_FLOAT32 ) setval((float*)src,              idx2,(float*)val,idx);
         else if( Dtype == E_CIO_FLOAT64 ) setval((double*)src,             idx2,(double*)val,idx);
      }}}}
    }
    if( !write_data(fp, src, DFI_Finfo.GuideCell, m_RankID) ) {
      fclose(fp);
      delete [] src;
      return;
    }
    delete [] src;
  }

  fclose(fp);

  WriteIndexDfiFile(m_indexDfiName,m_RankID,DFI_Finfo.Prefix,step,time,minmax,true,
                    air_mode, step_avr, time_avr);

}

// #################################################################
//
bool cio_DFI_BOV::write_data(FILE* fp, void* val, int gc, int n)
{

  E_CIO_DTYPE Dtype = get_cio_Datatype(DFI_Finfo.DataType);
  int Real_size = get_cio_Datasize(Dtype);

  int size[3];
  for(int i=0; i<3; i++ ) size[i] = (int)RankInfo[n].VoxelSize[i]+(int)(2*gc);

  size_t dLen = (size_t)(size[0] * size[1] * size[2]);
  if( DFI_Finfo.Component > 1 ) dLen *= 3;

  unsigned int dmy = dLen * Real_size;

  if( fwrite(val, Real_size, dLen, fp) != dLen ) return false;

  return true;
}
