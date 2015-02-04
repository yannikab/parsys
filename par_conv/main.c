/* 
 * File:   main.c
 * Author: John
 *
 * Created on December 24, 2014, 12:02 AM
 */

#include <stdio.h>
#include <stdlib.h>

int main_serial(int argc, char** argv);
int main_serial_omp(int argc, char** argv);

int main_sync(int argc, char** argv);
int main_sync_omp(int argc, char** argv);
int main_sync_omp_simple(int argc, char** argv);

int main_async(int argc, char** argv);
int main_async_nonper(int argc, char** argv);
int main_async_omp(int argc, char** argv);
int main_async_omp_simple(int argc, char** argv);

/*
 * 
 */
int main(int argc, char** argv)
{
    // return main_serial(argc, argv);
    // return main_serial_omp(argc, argv);

    return main_async(argc, argv);
    // return main_async_nonper(argc, argv);
    // return main_async_omp(argc, argv);
    // return main_async_omp_simple(argc, argv);

    // return main_sync(argc, argv);
    // return main_sync_omp(argc, argv);
    // return main_sync_omp_simple(argc, argv);
}
