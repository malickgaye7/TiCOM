#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#ifndef __TIRADE_H__
#define __TIRADE_H__
#include "tirade.h"
#endif

unsigned char check_container_members(Tiradata*);
unsigned char flags_to_Tiradata(char*);
Tiradata* construct_Tiradata(char*, int*);
void destroy_Tiradata(Tiradata*);

extern struct unitFlags uFlags;
extern struct matprop *layers;

long double getAvg(long double*, int, int, int);
long double linearApprox(long double, int, long double*, int, int);
long double layerData(int, char, char*, long double, int, int, long double*);
unsigned char flagSearch(char*);
void stitchLine(char*, char*);
void setUnits(char*, char*);
void convertUnits(long double*, char*, int, int);
long double* readCSVData(int, int, FILE*, char*, char*);

#endif /* __UTILITY_H__ */