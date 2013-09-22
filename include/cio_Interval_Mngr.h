#ifndef _CIO_INTERVAL_MNGR_H_
#define _CIO_INTERVAL_MNGR_H_


// 使用例
// 1. 初期化処理
//   (1) インスタンス
//     cio_Interval_Mngr intvl;
//   (2) モードをセット
//     intvl.setMode( cio_Interval_Mngr::By_step );
//   (3) インターバルステップ(時刻)のセット
//     intvl.setInterval( (double)stepInterval );
//   (4) セッション開始ステップ(時刻)のセット(未指定のとき0)
//     intvl.setStart( (double)startStep );
//   (5) 最終ステップ(時刻)のセット(指定の場合、計算の最終ステップは必ず出力する)
//     intvl.setLast( (double)lastStep );
//   (6) トリガーの初期化
//     intvl.initTrigger( step, time, dt );
// 2.ループ処理
//   (1) 出力ステップ(時刻)かどうか判定
//     bool out = intvl.isTriggerd( step, time );

class cio_Interval_Mngr
{
// enum
public:
  /// モード
  enum type_IO_spec
  {
    noset  =-1,
    By_step= 1,
    By_time= 2,
  };

// メンバ変数
protected:
  type_IO_spec m_mode; ///< モード(noset:未定義、By_step:ステップ間隔指定、By_time:時刻間隔指定)

  int    m_base_step;  ///< インターバルの基点となるステップ(By_stepのとき有効)
  int    m_intvl_step; ///< ステップ間隔(By_stepのとき有効)
  int    m_start_step; ///< セッションを開始するステップ(By_stepのとき有効)
  int    m_last_step;  ///< 最終ステップ(By_stepのとき有効)

  double m_base_time;  ///< インターバルの基点となる時刻(By_timeのとき有効)
  double m_intvl_time; ///< 時刻間隔(By_timeのとき有効)
  double m_start_time; ///< セッションを開始する時刻(By_timeのとき有効)
  double m_last_time;  ///< 最終時刻(By_timeのとき有効)

  double m_dt;         ///< 計算の時間間隔Δt

// メンバ関数
public:
  /// コンストラクタ
  cio_Interval_Mngr()
  {
    m_mode = noset;
    m_base_step  = 0;
    m_intvl_step = 0;
    m_start_step = 0;
    m_last_step  = -1;
    m_base_time  = 0.0;
    m_intvl_time = 0.0;
    m_start_time = 0.0;
    m_last_time  = -1.0;
    m_dt         = 0.0;
  }

  /// デストラクタ
  ~cio_Interval_Mngr()
  {
  }

  /// モードをセット
  void setMode( const type_IO_spec mode )
  {
    m_mode = mode;
  }

  /// モードを取得
  type_IO_spec getMode()
  {
    return m_mode;
  }

  /// インターバル値のセット
  bool setInterval( const double interval )
  {
    if( m_mode == By_step )
    {
      m_intvl_step = (int)interval;
    }
    else if( m_mode == By_time )
    {
      m_intvl_time = interval;
    }
    else
    {
      return false;
    }
    return true;
  }

  /// インターバル(ステップ)の取得
  int getIntervalStep()
  {
    return m_intvl_step;
  }

  /// インターバル(時刻)の取得
  double getIntervalTime()
  {
    return m_intvl_time;
  }

  /// 最終ステップ(時刻)のセット
  bool setLast( const double last )
  {
    if( m_mode == By_step )
    {
      m_last_step = (int)last;
    }
    else if( m_mode == By_time )
    {
      m_last_time = last;
    }
    else
    {
      return false;
    }
    return true;
  }

  /// セッション開始ステップ(時刻)のセット
  bool setStart( const double start )
  {
    if( m_mode == By_step )
    {
      m_start_step = (int)start;
    }
    else if( m_mode == By_time )
    {
      m_start_time = start;
    }
    else
    {
      return false;
    }
    return true;
  }

  /// セッション開始ステップの取得
  int getStartStep()
  {
    return m_start_step;
  }

  /// セッション開始時刻の取得
  double getStartTime()
  {
    return m_start_time;
  }

  /// トリガーの初期化
  bool initTrigger( const int step, const double time, const double dt, bool d_flag=false )
  {
    m_dt = dt;
    if( m_mode == By_step )
    {
      m_base_step = step;
      if( d_flag )
      {
        int next_step = calcNextStep(step);
        if( next_step == m_base_step ) next_step += m_intvl_step;
        printf("cio_Interval_Mngr : mode=By_step, base_step=%d, interval=%d, start=%d, first_step=%d"
              , m_base_step, m_intvl_step, m_start_step, next_step);
        if( m_last_step > 0 )
          printf(", last_step=%d\n", m_last_step);
        else
          printf("\n");
      }
    }
    else if( m_mode == By_time )
    {
      m_base_time = time;
      if( d_flag )
      {
        double next_time = calcNextTime(time);
        if( next_time == m_base_time ) next_time += m_intvl_time;
        printf("cio_Interval_Mngr : mode=By_time, base_time=%e, delta_t=%e, interval=%e, start=%e, first_time=%e"
              , m_base_time, m_dt, m_intvl_time, m_start_time, next_time);
        if( m_last_time > 0.0 )
          printf(", last_time=%e\n", m_last_time);
        else
          printf("\n");
      }
    }
    else
    {
      return false;
    }
    return true;
  }

  /// インターバルのチェック
  /// 出力ステップ(時刻)であっても、同じステップ(時刻)に既にisTriggeredがコールされている場合はfalseを返す(重複出力対応)
  bool isTriggered( const int step, const double time, bool forced_out=false, bool d_flag=false )
  {
    // 強制出力フラグ
    if( forced_out )
    {
      return true;
    }

    // nosetのときは必ずtrue(上位プログラムで判定しているものとする)
    if( m_mode == noset )
    {
      return true;
    }

    // セッションが開始しているか
    if( !isStarted(step, time) )
    {
      return false;
    }

    // ステップ指定のとき
    if( m_mode == By_step )
    {
      // intervalが正常か
      if( m_intvl_step <= 0 )
      {
        return false;
      }
      // 次の出力ステップかどうか
      if( step == calcNextStep( step ) )
      {
        if( d_flag )
        {
          printf("cio_Interval_Mngr::isTriggerd : step=%d\n", step);
        }
        return true;
      }
      // 最終ステップかどうか
      if( isLastStep( step ) )
      {
        if( d_flag )
        {
          printf("cio_Interval_Mngr::isTriggerd : last step=%d\n", step);
        }
        return true;
      }
    }
    // 時刻指定のとき
    else if( m_mode == By_time )
    {
      // intervalが正常か
      if( m_intvl_time <= 0.0 )
      {
        return false;
      }
      // 次の出力時刻
      double next_time = calcNextTime( time );
      // 次の出力時刻かどうか
      if( time <= next_time && next_time < time + m_dt )
      {
        if( d_flag )
        {
          printf("cio_Interval_Mngr::isTriggerd : time=%e\n", time);
        }
        return true;
      }
      // 最終時刻かどうか
      if( isLastTime( time ) )
      {
        if( d_flag )
        {
          printf("cio_Interval_Mngr::isTriggerd : last time=%e\n", time);
        }
        return true;
      }
    }
    return false;
  }

  /// セッションが開始しているかをチェック
  bool isStarted( const int step, const double time )
  {
    if( m_mode == By_step )
    {
      if( m_start_step <= step ) return true;
    }
    else if( m_mode == By_time )
    {
      if( m_start_time <= time ) return true;
    }
    return false;
  }

  // 時刻を指定のスケールで無次元化する
  bool normalizeTime( const double scale )
  {
    if( m_mode != By_time )
    {
      return false;
    }
    normalizeBaseTime(scale);
    normalizeIntervalTime(scale);
    normalizeStartTime(scale);
    normalizeLastTime(scale);
    normalizeDelteT(scale);
    return true;
  }
  void normalizeBaseTime(const double scale)
  {
    m_base_time /= scale;
  }
  void normalizeIntervalTime(const double scale)
  {
    m_intvl_time /= scale;
  }
  void normalizeStartTime(const double scale)
  {
    m_start_time /= scale;
  }
  void normalizeLastTime(const double scale)
  {
    m_last_time /= scale;
  }
  void normalizeDelteT(const double scale)
  {
    m_dt /= scale;
  }

protected:
  /// 次の出力ステップを計算
  int calcNextStep( const int step )
  {
    int s_step = (m_start_step > step) ? m_start_step : step;
    int inc    = ((s_step-m_base_step)%m_intvl_step==0) ? 0 : 1;
    int next_step = int((s_step-m_base_step)/m_intvl_step+inc) * m_intvl_step + m_base_step;
    return next_step;
  }

  /// 次の出力時刻を計算
  double calcNextTime( const double time )
  {
    //int s_time = (m_start_time > time) ? m_start_time : time;
    double s_time = (m_start_time > time) ? m_start_time : time;
    int inc    = (dmod(s_time-m_base_time, m_intvl_time)==0.0) ? 0 : 1;
    double next_time = int((s_time-m_base_time)/m_intvl_time+inc) * m_intvl_time + m_base_time;
    return next_time;
  }

  /// 実数の余り(fortranのmodと同じ)
  static double dmod( double a, double b )
  {
    return a - int(a/b) * b;
  }

  /// 最終ステップかどうか
  bool isLastStep( const int step )
  {
    if( m_last_step <= 0 ) return false;
    if( m_last_step == step ) return true;
    return false;
  }

  /// 最終時刻かどうか
  bool isLastTime( const double time )
  {
    if( m_last_time <= 0 ) return false;
    if( time <= m_last_time && m_last_time < time + m_dt ) return true;
    return false;
  }

 
};

#endif  /* _CIO_INTERVAL_MNGR_H_ */

