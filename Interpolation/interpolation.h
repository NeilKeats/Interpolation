#pragma once
#ifndef INTERPOLATION
#define INTERPOLATION

#include<Windows.h>   
#include<inttypes.h>
#include<string>

#define MODE_NEAREST_NEIGHBOUR 0x00
#define MODE_BICUBIC 0x01
#define MODE_SPLINE 0x02
#define MODE_BC_KERNEL 0x03
#define MODE_LANCZOS_KERNEL 0x04
#define MODE_BILINEAR 0x05

static int A_MODE[] = { 
	MODE_NEAREST_NEIGHBOUR,
	MODE_BICUBIC,
	MODE_SPLINE,
	MODE_BC_KERNEL,
	MODE_LANCZOS_KERNEL,
	MODE_BILINEAR};


static std::string MODE_NAME[] = { 
	"NN",
	"BICUBIC",
	"BSPLINE",
	"BC_KERNEL",
	"LANCZOS_KERNEL",
	"BILINEAR"};


float bubic_conv_kernel(float a, float x);

float lanczos_conv_kernel(float a, float x);

void gradientX(const float *Imdata, float * dx, float *dbuffer, DWORD width, DWORD height);

void gradientY(const float *Imdata, float * dx, float *dbuffer, DWORD width, DWORD height);

void bicubic_coeff(const float *Imdata,float * coeff, DWORD width, DWORD height);

void bicubic_spline_coeff(const float *f_data, float * coeff, DWORD s_width, DWORD s_height);

float cal_bicubic(float *coeff, float s_x, float s_y, DWORD s_width, DWORD s_height);

float cal_bicubic_kernel(const float *f_data, float s_x, float s_y, DWORD s_width, DWORD s_height, float a, int MODE);

float nearest_neighbour(const float *f_data, float s_x, float s_y, DWORD s_width, DWORD s_height);

float bilinear(const float* f_data, float s_x, float s_y, DWORD s_width, DWORD s_height);

void fill_data(const uint8_t* s_data, float * f_data, DWORD s_width, DWORD s_height);

void interpolation(const uint8_t *s_data, uint8_t *d_data, DWORD s_width, DWORD s_height, float width_scale, float height_scale, int MODE);


#endif