# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

//# define DEBUG

int main ( void );
int explode ( double x, double y, int count_max );
int ppma_write ( char *out_filename, int xsize, int ysize, int *r, int *g, int *b );
int ppma_write_data ( FILE *file_out, int xsize, int ysize, int *r, int *g, int *b );
int ppma_write_header ( FILE *file_out, char *out_filename, int xsize, int ysize, int rgb_max );
void timestamp ( void );
int explodeParallel ( double x, double y, int count_max );

int main ( void )
{
    int *b;
    int c;
    int c_max;
    int *count;
    int count_max = 1000;
    char *filename = "mandelbrot.ppm";
    int *g;
    int i;
    int ierror;
    int j;
    int k;
    int n = 5501;
    int *r;
    double x;
    double x_max =   1.25;
    double x_min = - 2.25;
    double y;
    double y_max =   1.75;
    double y_min = - 1.75;
    double TinicioSeq;
    double TinicioParallelOMP;
    double TinicioParallelVarPrivada;
    double TinicioParallelRumetime;
    double TtotalSeq;
    double TtotalParallelOMP;
    double TtotalParallelVarPrivada;
    double TtotalParallelRuntime;
    double TinicioBloque2Omp = 0.0;
    double TinicioBloque2Runtime = 0.0;
    double TinicioBloque2Secuencial = 0.0;
    double TinicioBloque2VarPrivada = 0.0;
    double TtotalBloque2Omp = 0.0;
    double TtotalBloque2Runtime = 0.0;
    double TtotalBloque2Secuencial = 0.0;
    double TtotalBloque2VarPrivada = 0.0;

    timestamp ( );
    printf ( "\n" );
    printf ( "MANDELBROT\n" );
    printf ( "  C version\n" );
    printf ( "\n" );
    printf ( "  Create an ASCII PPM image of the Mandelbrot set.\n" );
    printf ( "\n" );
    printf ( "  For each point C = X + i*Y\n" );
    printf ( "  with X range [%f,%f]\n", x_min, x_max );
    printf ( "  and  Y range [%f,%f]\n", y_min, y_max );
    printf ( "  carry out %d iterations of the map\n", count_max );
    printf ( "  Z(n+1) = Z(n)^2 + C.\n" );
    printf ( "  If the iterates stay bounded (norm less than 2)\n" );
    printf ( "  then C is taken to be a member of the set.\n" );
    printf ( "\n" );
    printf ( "  An ASCII PPM image of the set is created using\n" );
    printf ( "    N = %d pixels in the X direction and\n", n );
    printf ( "    N = %d pixels in the Y direction.\n", n );
    /* -----------------------Sequencial---------------------------
     SETUP DATA
     */
    count = ( int * ) malloc ( n * n * sizeof ( int ) );
    r = ( int * ) malloc ( n * n * sizeof ( int ) );
    g = ( int * ) malloc ( n * n * sizeof ( int ) );
    b = ( int * ) malloc ( n * n * sizeof ( int ) );
    /*
    Carry out the iteration for each pixel, determining COUNT.
    */

    printf ( "--------Secuencial---------\n" );

    TinicioSeq = omp_get_wtime();

    for ( i = 0; i < n; i++ )
    {
    for ( j = 0; j < n; j++ )
    {
      x = ( ( double ) (     j     ) * x_max
          + ( double ) ( n - j - 1 ) * x_min )
          / ( double ) ( n     - 1 );

      y = ( ( double ) (     i     ) * y_max
          + ( double ) ( n - i - 1 ) * y_min )
          / ( double ) ( n     - 1 );

      count[i+j*n] = explode ( x, y, count_max );
    }
    }
    /*
    Determine the coloring of each pixel.
    */
    TinicioBloque2Secuencial = omp_get_wtime();
        c_max = 0;
        for ( j = 0; j < n; j++ )
        {
            for ( i = 0; i < n; i++ )
            {
                if ( c_max < count[i+j*n] )
                {
                c_max = count[i+j*n];
                }
            }
        }
    TtotalBloque2Secuencial = omp_get_wtime() - TinicioBloque2Secuencial;
    #ifdef DEBUG
            printf("c_max Secuencial = %d\n", c_max);
    #endif

    /*
    Set the image data.
    */

    for ( i = 0; i < n; i++ )
    {
    for ( j = 0; j < n; j++ )
    {
      if ( count[i+j*n] % 2 == 1 )
      {
        r[i+j*n] = 255;
        g[i+j*n] = 255;
        b[i+j*n] = 255;
      }
      else
      {
        c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt (
          ( ( double ) ( count[i+j*n] ) / ( double ) ( c_max ) ) ) ) ) );
        r[i+j*n] = 3 * c / 5;
        g[i+j*n] = 3 * c / 5;
        b[i+j*n] = c;
      }
    }
    }

    TtotalSeq = omp_get_wtime() - TinicioSeq;

    /*
    Write an image file.
    */
    ierror = ppma_write ( filename, n, n, r, g, b );

    printf ( "\n" );
    printf ( "  ASCII PPM image data stored in \"%s\".\n", filename );

    free ( b );
    free ( count );
    free ( g );
    free ( r );
    /*
    Terminate.
    */
    printf ( "\n" );
    printf ( "MANDELBROT\n" );
    printf ( "  Normal end of execution. Duracion Secuencial = %f\n", TtotalSeq );
    printf ( "\n" );


    /* -----------------------Paralelo OMP---------------------------
     SETUP DATA
     */
    count = ( int * ) malloc ( n * n * sizeof ( int ) );
    r = ( int * ) malloc ( n * n * sizeof ( int ) );
    g = ( int * ) malloc ( n * n * sizeof ( int ) );
    b = ( int * ) malloc ( n * n * sizeof ( int ) );
    filename = "mandelbrotParallel.ppm";
    /*
    Carry out the iteration for each pixel, determining COUNT.
    */

    printf ( "---------Paralelo OMP--------\n" );

    TinicioParallelOMP = omp_get_wtime();

    #pragma omp parallel firstprivate(n) private(j, i, x, y, c) shared(x_max, x_min, y_max, y_min, count, r, b, g)
    {
        #pragma omp for
            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    #pragma omp task
                    {
                        x = ( ( double ) (     j     ) * x_max
                          + ( double ) ( n - j - 1 ) * x_min )
                          / ( double ) ( n     - 1 );

                        y = ( ( double ) (     i     ) * y_max
                          + ( double ) ( n - i - 1 ) * y_min )
                          / ( double ) ( n     - 1 );

                        count[i+j*n] = explode ( x, y, count_max );
                    }
                }
            }



        /*
        Determine the coloring of each pixel.
        */
        #pragma omp single
        {
            TinicioBloque2Omp = omp_get_wtime();
            c_max = 0;
        }
        #pragma omp for
            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < n; i++ )
                {

                    if ( c_max < count[i+j*n] )
                    {
                        #pragma omp critical
                        {
                            c_max = count[i+j*n];
                        }
                    }

                }
            }
        #pragma omp single
        {
            TtotalBloque2Omp = omp_get_wtime() - TinicioBloque2Omp;
            #ifdef DEBUG
            printf("c_max OMP = %d\n", c_max);
            #endif
        }

        /*
        Set the image data.
        */
        #pragma omp for
        for ( i = 0; i < n; i++ )
        {
        for ( j = 0; j < n; j++ )
        {
          if ( count[i+j*n] % 2 == 1 )
          {
            r[i+j*n] = 255;
            g[i+j*n] = 255;
            b[i+j*n] = 255;
          }
          else
          {
            c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt (
              ( ( double ) ( count[i+j*n] ) / ( double ) ( c_max ) ) ) ) ) );
            r[i+j*n] = 3 * c / 5;
            g[i+j*n] = 3 * c / 5;
            b[i+j*n] = c;
          }
        }
        }
    }

    TtotalParallelOMP = omp_get_wtime() - TinicioParallelOMP;

    /*
    Write an image file.
    */
    ierror = ppma_write ( filename, n, n, r, g, b );

    printf ( "\n" );
    printf ( "  ASCII PPM image data stored in \"%s\".\n", filename );

    free ( b );
    free ( count );
    free ( g );
    free ( r );
    /*
    Terminate.
    */
    printf ( "\n" );
    printf ( "MANDELBROT\n" );
    printf ( "  Normal end of execution. Duracion Paralelo = %f\n  SpeedUp = %f\n", TtotalParallelOMP, TtotalSeq/TtotalParallelOMP);
    printf ( "\n" );



    /* -----------------------Paralelo Variables privadas---------------------------
     SETUP DATA
     */
    count = ( int * ) malloc ( n * n * sizeof ( int ) );
    r = ( int * ) malloc ( n * n * sizeof ( int ) );
    g = ( int * ) malloc ( n * n * sizeof ( int ) );
    b = ( int * ) malloc ( n * n * sizeof ( int ) );
    int * c_maxmulti;
    int numThreads;
    filename = "mandelbrotParallelVarPrivadas.ppm";
    /*
    Carry out the iteration for each pixel, determining COUNT.
    */

    printf ( "---------Paralelo Variables Privadas--------\n" );

    TinicioParallelVarPrivada = omp_get_wtime();

    #pragma omp parallel firstprivate(n) private(j, i, x, y, c) shared(x_max, x_min, y_max, y_min, count, r, b, g, c_maxmulti)
    {
        #pragma omp for
            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    #pragma omp task
                    {
                        x = ( ( double ) (     j     ) * x_max
                          + ( double ) ( n - j - 1 ) * x_min )
                          / ( double ) ( n     - 1 );

                        y = ( ( double ) (     i     ) * y_max
                          + ( double ) ( n - i - 1 ) * y_min )
                          / ( double ) ( n     - 1 );

                        count[i+j*n] = explode ( x, y, count_max );
                    }
                }
            }

        /*
        Determine the coloring of each pixel.
        */

        //Usando variables privadas
        #pragma omp single
        {
            TinicioBloque2VarPrivada = omp_get_wtime();
            numThreads = omp_get_max_threads();
            c_maxmulti = (int *) malloc(sizeof(int) * numThreads);
        }
        c_maxmulti[omp_get_thread_num()] = 0;
        #ifdef DEBUG
            printf ("soy el %d\n", omp_get_thread_num());
            printf("array inicializado %d = %d\n", omp_get_thread_num(), c_maxmulti[omp_get_thread_num()]);
        #endif
        #pragma omp for
            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < n; i++ )
                {
                    if ( c_maxmulti[omp_get_thread_num()] < count[i+j*n] )
                    {

                        c_maxmulti[omp_get_thread_num()] = count[i+j*n];
                    }
                }
            }

        #pragma omp barrier
        #pragma omp single
        {
            c_max = 0;
            for (i = 0; i < numThreads; i++){
                #ifdef DEBUG
                    printf("c_max Vprivada %d = %d\n", i, c_maxmulti[i]);
                #endif
                if ( c_max < c_maxmulti[i])
                    {
                        c_max = c_maxmulti[i];
                    }
            }


            TtotalBloque2VarPrivada = omp_get_wtime() - TinicioBloque2VarPrivada;
            #ifdef DEBUG
                printf("c_max Vprivada = %d\n", c_max);
            #endif
        }

        /*
        Set the image data.
        */
        #pragma omp for
        for ( i = 0; i < n; i++ )
        {
        for ( j = 0; j < n; j++ )
        {
          if ( count[i+j*n] % 2 == 1 )
          {
            r[i+j*n] = 255;
            g[i+j*n] = 255;
            b[i+j*n] = 255;
          }
          else
          {
            c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt (
              ( ( double ) ( count[i+j*n] ) / ( double ) ( c_max ) ) ) ) ) );
            r[i+j*n] = 3 * c / 5;
            g[i+j*n] = 3 * c / 5;
            b[i+j*n] = c;
          }
        }
        }
    }

    TtotalParallelVarPrivada = omp_get_wtime() - TinicioParallelVarPrivada;

    /*
    Write an image file.
    */
    ierror = ppma_write ( filename, n, n, r, g, b );

    printf ( "\n" );
    printf ( "  ASCII PPM image data stored in \"%s\".\n", filename );

    free ( b );
    free ( count );
    free ( g );
    free ( r );

    /*
    Terminate.
    */
    printf ( "\n" );
    printf ( "MANDELBROT\n" );
    printf ( "  Normal end of execution. Duracion Paralelo Variables privadas = %f\n  SpeedUp = %f\n", TtotalParallelVarPrivada, TtotalSeq/TtotalParallelVarPrivada);
    printf ( "\n" );


/* -----------------------Paralelo Funciones Runtime ---------------------------
     SETUP DATA
     */
    count = ( int * ) malloc ( n * n * sizeof ( int ) );
    r = ( int * ) malloc ( n * n * sizeof ( int ) );
    g = ( int * ) malloc ( n * n * sizeof ( int ) );
    b = ( int * ) malloc ( n * n * sizeof ( int ) );
    omp_lock_t mutex;
    filename = "mandelbrotParallelRuntime.ppm";
    /*
    Carry out the iteration for each pixel, determining COUNT.
    */

    printf ( "---------Paralelo Funciones Runtime--------\n" );

    TinicioParallelRumetime = omp_get_wtime();

    #pragma omp parallel firstprivate(n) private(j, i, x, y, c) shared(x_max, x_min, y_max, y_min, count, r, b, g, mutex)
    {
        #pragma omp for
            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    #pragma omp task
                    {
                        x = ( ( double ) (     j     ) * x_max
                          + ( double ) ( n - j - 1 ) * x_min )
                          / ( double ) ( n     - 1 );

                        y = ( ( double ) (     i     ) * y_max
                          + ( double ) ( n - i - 1 ) * y_min )
                          / ( double ) ( n     - 1 );

                        count[i+j*n] = explode ( x, y, count_max );
                    }
                }
            }



        /*
        Determine the coloring of each pixel.
        */
        #pragma omp single
        {
            TinicioBloque2Runtime = omp_get_wtime();
            c_max = 0;
            omp_init_lock(&mutex);

        }
        #pragma omp for
            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < n; i++ )
                {

                    if ( c_max < count[i+j*n] )
                    {
                        omp_set_lock(&mutex);
                        c_max = count[i+j*n];
                        omp_unset_lock(&mutex);
                    }
                }

            }

        #pragma omp single
        {
            omp_destroy_lock(&mutex);
            TtotalBloque2Runtime = omp_get_wtime() - TinicioBloque2Runtime;
            #ifdef DEBUG
                printf("c_max Runtime = %d\n", c_max);
            #endif
        }

        /*
        Set the image data.
        */
        #pragma omp for
        for ( i = 0; i < n; i++ )
        {
        for ( j = 0; j < n; j++ )
        {
          if ( count[i+j*n] % 2 == 1 )
          {
            r[i+j*n] = 255;
            g[i+j*n] = 255;
            b[i+j*n] = 255;
          }
          else
          {
            c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt (
              ( ( double ) ( count[i+j*n] ) / ( double ) ( c_max ) ) ) ) ) );
            r[i+j*n] = 3 * c / 5;
            g[i+j*n] = 3 * c / 5;
            b[i+j*n] = c;
          }
        }
        }
    }

    TtotalParallelRuntime = omp_get_wtime() - TinicioParallelRumetime;

    /*
    Write an image file.
    */
    ierror = ppma_write ( filename, n, n, r, g, b );

    printf ( "\n" );
    printf ( "  ASCII PPM image data stored in \"%s\".\n", filename );

    free ( b );
    free ( count );
    free ( g );
    free ( r );
    /*
    Terminate.
    */
    printf ( "\n" );
    printf ( "MANDELBROT\n" );
    printf ( "  Normal end of execution. Duracion Paralelo Rumetime = %f\n  SpeedUp = %f\n", TtotalParallelRuntime, TtotalSeq/TtotalParallelRuntime);
    printf ( "\n" );
    printf ( "\n" );
    printf ( "---------Resumen Resultados--------\n" );
    printf ( "\n" );
    printf ( "MANDELBROT\n" );
     printf ( "  Normal end of execution.\n ");
     printf ( "  Tiempos y Speedups \n ");
     printf("   Duracion Secuencial = %f  SpeedUp = %f\n", TtotalSeq, TtotalSeq/TtotalSeq);
     printf("    Duracion Paralelo OMP = %f  SpeedUp = %f\n", TtotalParallelOMP, TtotalSeq/TtotalParallelOMP);
     printf("    Duracion Paralelo Runtime = %f  SpeedUp = %f\n", TtotalParallelRuntime, TtotalSeq/TtotalParallelRuntime);
     printf("    Duracion Paralelo Variables Privadas = %f  SpeedUp = %f\n", TtotalParallelVarPrivada, TtotalSeq/TtotalParallelVarPrivada);
     printf ( "   Time to determine the coloring of each pixel.\n ");
     printf ("   Directivas Omp = %f\n    Funciones runtime = %f\n    Secuencial = %f\n    Con Variables privadas = %f\n", TtotalBloque2Omp, TtotalBloque2Runtime, TtotalBloque2Secuencial, TtotalBloque2VarPrivada);
    printf ( "\n" );






    timestamp ( );

    return 0;
}

/*int explodeParallel ( double x, double y, int count_max )
{
    int k;
    int value;
    double x1;
    double x2;
    double y1;
    double y2;

    //NO SE COMO PARALLELIZAR CON EL BRAKE
    value = 0;

    x1 = x;
    y1 = y;
    bool flag = true;

    #pragma omp parallel private(k,x2,x1,y1,y2) shared(flag)
    for ( k = 1; k <= count_max; k++ )
    {
        x2 = x1 * x1 - y1 * y1 + x;
        y2 = 2.0 * x1 * y1 + y;

        #pragma omp critical
        {
            if(flag) {
                if ( x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2 )
                {
                    flag = false;
                    value = k;

                }
            }
        }
        x1 = x2;
        y1 = y2;
    }
    return value;
}*/


int explode ( double x, double y, int count_max )
{
  int k;
  int value;
  double x1;
  double x2;
  double y1;
  double y2;

  value = 0;

  x1 = x;
  y1 = y;

  for ( k = 1; k <= count_max; k++ )
  {
    x2 = x1 * x1 - y1 * y1 + x;
    y2 = 2.0 * x1 * y1 + y;

    if ( x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2 )
    {
      value = k;
      break;
    }
    x1 = x2;
    y1 = y2;
  }
  return value;
}

int ppma_write ( char *filename, int xsize, int ysize, int *r, int *g, int *b )
{
  int *b_index;
  int error;
  FILE *file_out;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_max;
/*
  Open the output file.
*/
  file_out = fopen ( filename, "wt" );

  if ( !file_out )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  Cannot open the output file \"%s\".\n", filename );
    return 1;
  }
/*
  Compute the maximum.
*/
  rgb_max = 0;
  r_index = r;
  g_index = g;
  b_index = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( rgb_max < *r_index )
      {
        rgb_max = *r_index;
      }
      r_index = r_index + 1;

      if ( rgb_max < *g_index )
      {
        rgb_max = *g_index;
      }
      g_index = g_index + 1;

      if ( rgb_max < *b_index )
      {
        rgb_max = *b_index;
      }
      b_index = b_index + 1;
    }
  }
/*
  Write the header.
*/
  error = ppma_write_header ( file_out, filename, xsize, ysize, rgb_max );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  PPMA_WRITE_HEADER failed.\n" );
    return 1;
  }
/*
  Write the data.
*/
  error = ppma_write_data ( file_out, xsize, ysize, r, g, b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  PPMA_WRITE_DATA failed.\n" );
    return 1;
  }
/*
  Close the file.
*/
  fclose ( file_out );

  return 0;
}

int ppma_write_data ( FILE *file_out, int xsize, int ysize, int *r, int *g, int *b )
{
  int *b_index;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_num;

  r_index = r;
  g_index = g;
  b_index = b;
  rgb_num = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      fprintf ( file_out, "%d  %d  %d", *r_index, *g_index, *b_index );
      rgb_num = rgb_num + 3;
      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;

      if ( rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize )
      {
        fprintf ( file_out, "\n" );
      }
      else
      {
        fprintf ( file_out, " " );
      }
    }
  }
  return 0;
}

int ppma_write_header ( FILE *file_out, char *filename, int xsize, int ysize, int rgb_max )
{
  fprintf ( file_out, "P3\n" );
  fprintf ( file_out, "# %s created by ppma_write.c.\n", filename );
  fprintf ( file_out, "%d  %d\n", xsize, ysize );
  fprintf ( file_out, "%d\n", rgb_max );

  return 0;
}

void timestamp ( void )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
