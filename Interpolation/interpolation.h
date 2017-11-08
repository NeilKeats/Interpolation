#pragma once
#ifndef INTERPOLATION
#define INTERPOLATION

#include<Windows.h>   
#include<inttypes.h>

#define MODE_BICUBIC 0x00
#define MODE_SPLINE 0x0F
#define MODE_BC_KERNEL 0xF0
#define MODE_LANCZOS_KERNEL 0xFF
#define MODE_NEAREST_NEIGHBOUR 0x07

float bubic_conv_kernel(float a, float x);

float lanczos_conv_kernel(float a, float x);

void gradientX(const float *Imdata, float * dx, float *dbuffer, DWORD weight, DWORD hight);

void gradientY(const float *Imdata, float * dx, float *dbuffer, DWORD weight, DWORD hight);

void bicubic_coeff(const float *Imdata,float * coeff, DWORD weight, DWORD hight);

void bicubic_spline_coeff(const float *f_data, float * coeff, DWORD s_weight, DWORD s_hight);

float cal_bicubic(float *coeff, float s_x, float s_y, DWORD s_weight, DWORD s_hight);

float cal_bicubic_kernel(const float *f_data, float s_x, float s_y, DWORD s_weight, DWORD s_hight, float a, int MODE);

float nearest_neighbour(const float *f_data, float s_x, float s_y, DWORD s_weight, DWORD s_hight);

void fill_data(const uint8_t* s_data, float * f_data, DWORD s_weight, DWORD s_hight);

void interpolation(const uint8_t *s_data, uint8_t *d_data, DWORD s_weight, DWORD s_hight, float weight_scale, float hight_scale, int MODE);

#endif // !1
