/* 
 * File:   main.c
 * Author: John
 *
 * Created on December 24, 2014, 12:02 AM
 */

#include <stdio.h>
#include <stdlib.h>

/*
 * 
 */
int main(int argc, char** argv)
{
	int main_serial(int argc, char** argv);
	int main_sync(int argc, char** argv);
	int main_async(int argc, char** argv);

	// return main_serial(argc, argv);
	// return main_sync(argc, argv);
	return main_async(argc, argv);
}
