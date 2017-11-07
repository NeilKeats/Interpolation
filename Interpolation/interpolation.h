#pragma once
#ifndef INTERPOLATION
#define INTERPOLATION

#include<Windows.h>   

#define MODE_BICUBIC 0x00
#define MODE_SPLINE 0x0F

void gradientX(const float *Imdata, float * dx, float *dbuffer, DWORD weight, DWORD hight);

void gradientY(const float *Imdata, float * dx, float *dbuffer, DWORD weight, DWORD hight);

void bicubic_coeff(const float *Imdata,float * coeff, DWORD weight, DWORD hight);

float cal_bicubic(float *coeff,float s_x, float s_y);

void bicubic_spline_coeff(const float *f_data, float * coeff, DWORD s_weight, DWORD s_hight);

void fill_data(const char* s_data, float * f_data, DWORD s_weight, DWORD s_hight);

void interpolation(const char *s_data, char *d_data, DWORD s_weight, DWORD s_hight, int weight_scale, int hight_scale, int MODE);

#endif // !1
