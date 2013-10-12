!########################################################################
!
! CIOlib - Cartesian Input / Output library
!
! Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
! All rights reserved.
!
!########################################################################

!! *********************************************************************
!! szS srcの実ボクセル数
!! gcS srcの仮想セル数
!! szD dstの実ボクセル数(=szS*2)
!! gcD dstの仮想セル数(=gcS*2)
!! nc  成分数
!! src 粗データの配列
!! dst 密データの配列

  do n=1,nc
  do k=1-gcS,szS(3)+gcS
    kk=(k-1)*2+1
  do j=1-gcS,szS(2)+gcS
    jj=(j-1)*2+1
  do i=1-gcS,szS(1)+gcS
    ii=(i-1)*2+1

    q = src(i,j,k,n)

    dst(ii  ,jj  ,kk  ,n) = q
    dst(ii+1,jj  ,kk  ,n) = q
    dst(ii  ,jj+1,kk  ,n) = q
    dst(ii+1,jj+1,kk  ,n) = q
    dst(ii  ,jj  ,kk+1,n) = q
    dst(ii+1,jj  ,kk+1,n) = q
    dst(ii  ,jj+1,kk+1,n) = q
    dst(ii+1,jj+1,kk+1,n) = q

  enddo
  enddo
  enddo
  enddo

