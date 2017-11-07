#include "interpolation.h"
#include<iostream>  
#include<malloc.h>  
#include<stdlib.h>  
#include<stdio.h>  
#include<string.h> 
#include<math.h>

static float m_dBicubicMatrix[16][16] = { 
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ -3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0 },
	{ -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0 },
	{ 9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1 },
	{ -6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1 },
	{ 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0 },
	{ -6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1 },
	{ 4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1 } };

static double m_dBSplineControlMatrix[4][4] = {
	{ 71 / 56.0, -19 / 56.0, 5 / 56.0, -1 / 56.0 },
	{ -19 / 56.0, 95 / 56.0, -25 / 56.0, 5 / 56.0 },
	{ 5 / 56.0, -25 / 56.0, 95 / 56.0, -19 / 56.0 },
	{ -1 / 56.0, 5 / 56.0, -19 / 56.0, 71 / 56.0 } };

static double m_dBSplineFunctionMatrix[4][4] = { 
	{ -1 / 6.0, 3 / 6.0, -3 / 6.0, 1 / 6.0 },
	{ 3 / 6.0, -6 / 6.0, 3 / 6.0, 0 },
	{ -3 / 6.0, 0, 3 / 6.0, 0 },
	{ 1 / 6.0, 4 / 6.0, 1 / 6.0, 0 } };

void gradientX(const float *f_data, float * dx, float *dbuffer, DWORD weight, DWORD hight) {
	const DWORD local_weight = weight + 4;
	for(int i= 0; i<hight+4;++i)
		for (int j = 2; j < weight + 2; ++j) {
			dbuffer[local_weight * i + j] =	 
				2.0 / 3 * f_data[local_weight * i + j + 1] 
				- 1.0 / 12 * f_data[local_weight * i + j + 2]
				- 2.0 / 3 * f_data[local_weight * i + j -1]
				+ 1.0 / 12 * f_data[local_weight * i + j -2 ];
		}	

	for(int i=0; i<hight; ++i)
		for (int j = 0; j < weight; ++j) {
			int row = i + 2;
			int col = j + 2;
			dx[i*weight + j] = dbuffer[row*local_weight+col];
		}
}

void gradientY(const float *f_data, float * dy, float *dbuffer, DWORD weight, DWORD hight) {
	const DWORD local_weight = weight + 4;
	for( int i =2 ; i< hight+2 ; ++i)
		for (int j = 0; j < weight; ++j) {
			dbuffer[local_weight * i + j] =   
				2.0 / 3 * f_data[local_weight * (i - 1 ) + j ]
				- 1.0 / 12 * f_data[local_weight * (i - 2 ) + j ]
				- 2.0 / 3 * f_data[local_weight * (i + 1 ) + j ]
				+ 1.0 / 12 * f_data[local_weight * (i + 2 ) + j ];
		}

	for (int i = 0; i<hight; ++i)
		for (int j = 0; j < weight; ++j) {
			int row = i + 2;
			int col = j + 2;
			dy[i*weight + j] = dbuffer[row*local_weight + col];
		}
}

void bicubic_coeff(const float *f_data, float * coeff, DWORD s_weight, DWORD s_hight){
	
	float *dx, *dy, *dxy,*d_buffx, *d_buffy;
	dx = (float*)malloc(sizeof(float)*s_weight*s_hight);
	dy = (float*)malloc(sizeof(float)*s_weight*s_hight);
	dxy = (float*)malloc(sizeof(float)*s_weight*s_hight);
	d_buffx = (float*)malloc(sizeof(float)*(s_weight + 4)*(s_hight + 4));
	d_buffy = (float*)malloc(sizeof(float)*(s_weight + 4)*(s_hight + 4));

	gradientY(f_data, dy, d_buffy, s_weight, s_hight);
	gradientX(f_data, dx, d_buffx, s_weight, s_hight);
	gradientY(d_buffx, dxy, d_buffy, s_weight, s_hight);

	float dTao[16], dAlpha[16];

	for(int i = 1 ; i<s_hight ; ++i ){
		for (int j = 0; j < s_weight - 1; ++j) {
			dTao[0] = f_data[(i + 2)*(s_weight + 4) + j + 2];
			dTao[1] = f_data[(i + 2)*(s_weight + 4) + j + 3];
			dTao[2] = f_data[(i + 1)*(s_weight + 4) + j + 2];
			dTao[3] = f_data[(i + 1)*(s_weight + 4) + j + 3];

			
			dTao[4] = dx[(i)*s_weight  + j ];
			dTao[5] = dx[(i)*s_weight  + j + 1];
			dTao[6] = dx[(i - 1)*s_weight  + j ];
			dTao[7] = dx[(i - 1)*s_weight  + j + 1 ];

			dTao[8] = dy[(i)*s_weight + j ];
			dTao[9] = dy[(i)*s_weight + j + 1];
			dTao[10] = dy[(i - 1)*s_weight + j];
			dTao[11] = dy[(i - 1)*s_weight + j + 1];
			
			dTao[12] = dxy[(i)*s_weight + j ];
			dTao[13] = dxy[(i)*s_weight + j + 1];
			dTao[14] = dxy[(i - 1)*s_weight + j];
			dTao[15] = dxy[(i - 1)*s_weight + j + 1];
			
			/*
			dTao[4] = dTao[5] = dTao[6] = dTao[7] = 0;
			dTao[8] = dTao[9] = dTao[10] = dTao[11] = 0;
			dTao[12] = dTao[13] = dTao[14] = dTao[15] = 0;
			*/

			for (int k = 0; k < 16; k++)
			{
				dAlpha[k] = 0;
				for (int l = 0; l < 16; l++)
				{
					dAlpha[k] += (m_dBicubicMatrix[k][l] * dTao[l]);
				}
			}

			float *fBicubic = coeff + ((i*s_weight) + j) * 16;
			for (int i = 0; i < 16; i++)
				fBicubic[i] = dAlpha[i];

		}
	}

	free(dx);
	free(dy);
	free(dxy);
	free(d_buffx);
	free(d_buffy);

}

float cal_bicubic(float *local_coeff, float x, float y) {
	/*
	float x = s_x - floor(s_x);
	float y;
	if (s_y <= 0.0)
		y = 1.0;
	else
		y = ceil(s_y) - s_y;
	*/

	float p_x[4],temp[4];
	/*
	p_x[0] = 1;
	p_x[1] = x;
	p_x[2] = x*x;
	p_x[3] = x*x*x;
	*/
	p_x[0] = pow(x, 0);
	p_x[1] = pow(x, 1);
	p_x[2] = pow(x, 2);
	p_x[3] = pow(x, 3);

	//
	float result = 0;
	temp[0] = pow(y, 0);
	temp[1] = pow(y, 1);
	temp[2] = pow(y, 2);
	temp[3] = pow(y, 3);
	for(int k=0 ; k<4 ; ++k)
		for (int j = 0; j < 4; ++j) {
			result += local_coeff[k*4 + j] * temp[k] * p_x[j];
			//result += local_coeff[k * 4 + j] * temp[j] * p_x[k];
		}

	/*
	for(int i=0;i<4;++i)
		temp[i] = p_x[0] * local_coeff[0+i] 
				+ p_x[1] * local_coeff[4+i] 
				+ p_x[2] * local_coeff[8+i] 
				+ p_x[3] * local_coeff[12+i];
	*/

		/*
		temp[i] = p_x[0] * local_coeff[i * 4 + 1]
				+ p_x[1] * local_coeff[i * 4 + 2]
				+ p_x[2] * local_coeff[i * 4 + 3]
				+ p_x[3] * local_coeff[i * 4 + 4];
		*/
	//return (temp[0] + temp[1]*y + temp[2]*y*y + temp[3]*y*y*y);
	return result;
}

void bicubic_spline_coeff(const float *f_data, float * coeff, DWORD s_weight, DWORD s_hight) {
	
	int i, j, k, l, m, n;
	float  Omiga[4][4], Beta[4][4];

	for (i = 0; i < s_hight; ++i) {
		for (j = 0; j < s_weight; ++j) {
			float *local_coeff = coeff + (i*s_weight +j)*16;

			//store neighbour value into Omeiga
			for (k = 0; k < 4; ++k)
				for (l = 0; l < 4; ++l)
					//reflect to original data (i + 1 - k, j - 1 +l)
					//reflect to filled data (i + 1 -k +2 , j - 1 +l +2)
					Omiga[k][l] = f_data[(i+1-k+2)*(s_weight+4) + (j-1+l+2)] ;

			//Beta
			for (k = 0; k < 4; k++)
			{
				for (l = 0; l < 4; l++)
				{
					Beta[k][l] = 0;
					for (m = 0; m < 4; m++)
					{
						for (n = 0; n < 4; n++)
						{
							Beta[k][l] += 
								m_dBSplineControlMatrix[k][m] *
								m_dBSplineControlMatrix[l][n] * 
								Omiga[n][m];
						}
					}
				}
			}

			//calculate p array;
			for (k = 0; k < 4; k++)
			{
				for (l = 0; l < 4; l++)
				{
					//dTBspline[i][j][k][l] = 0;
					local_coeff[k*4 + l] = 0 ;
					for (m = 0; m < 4; m++)
					{
						for (n = 0; n < 4; n++)
						{
							//dTBspline[i][j][k][l] += m_dBSplineFunctionMatrix[k][m] * m_dBSplineFunctionMatrix[l][n] * dBeta[n][m];
							local_coeff[k * 4 + l] += m_dBSplineFunctionMatrix[k][m] * m_dBSplineFunctionMatrix[l][n] * Beta[n][m];
						}
					}
				}
			}

			//trans p to a;
			for (k = 0; k < 2; k++)
			{
				for (l = 0; l < 4; l++)
				{
					float m_dTemp = local_coeff[k * 4 + l];
					local_coeff[k * 4 + l] = local_coeff[(3-k) * 4 + (3-l)];
					local_coeff[(3 - k) * 4 + (3 - l)] = m_dTemp;
				}
			}

		}
	
	}

}

void fill_data(const char* s_data, float * f_data, DWORD s_weight, DWORD s_hight) {
//将原 w*h 大小的图像 填充到 (w+4) * (h+4)中，即上下左右增加了两行数据。
//原边界外的元素的值用近邻元素的像素值代替
	DWORD f_weight = s_weight + 4;
	DWORD f_hight = s_hight +4 ;

	int s_i, s_j;
	for(int i=0; i < f_hight ;++i)
		for (int j = 0; j < f_weight; ++j) {
			
			s_i = (i - 2) < 0 ? 0 : (i - 2);
			s_i = (i + 2) >= s_hight ? s_hight-1 : s_i;

			s_j = (j - 2) < 0 ? 0 : (j - 2);
			s_j = (j + 2) >= s_weight ? s_weight - 1 : s_j;
			f_data[i*f_weight+j] = (float)s_data[s_i * s_weight + s_j];
			
			/*
			if (i < 2 || (i + 2) >= s_hight || j < 2 || j + 2 >= s_weight)
				f_data[i*f_weight + j] = 0;
			else
				f_data[i*f_weight + j] = (float)s_data[(i-2) * s_weight + (j-2)];
			*/
		}
}

void interpolation(const char *s_data, char *d_data, DWORD s_weight, DWORD s_hight, float weight_scale, float hight_scale, int MODE) {
	//prefix: d for destination  image; s for source image
	//
	if (weight_scale <= 0 || hight_scale <= 0)
		return;
	DWORD d_weight = s_weight * weight_scale;
	DWORD d_hight = s_hight * hight_scale;

	//buffer with filling
	float *f_data = nullptr;
	f_data = (float *)malloc(sizeof(float)*(s_weight + 4)*(s_hight + 4));
	fill_data(s_data, f_data, s_weight, s_hight);

	//coefficient
	float *coeff;
	coeff = (float *)malloc(sizeof(float)*s_weight*s_hight * 16);
	//compute coeef
	if (MODE == MODE_BICUBIC)
		bicubic_coeff(f_data, coeff, s_weight, s_hight);
	else if (MODE == MODE_SPLINE)
		bicubic_spline_coeff(f_data, coeff, s_weight, s_hight);

	char *output = d_data;
	float s_x, s_y;
	//is_x  int s_x;  is_y  int s_y
	int is_x, is_y;
	float temp;
	for (int i = 0; i < d_hight; ++i) {
		//coefficient array obtained from the left corner element.
		//for y, using ceil(),causes our image data start from left bottom
		//for x, using floor()

		//destination y ->  source y
		s_y = (float)i / (float)(d_hight - 1) * (float)(s_hight - 1);
		is_y = ceil(s_y);
		is_y = is_y <= 0 ? 1 : is_y;
		for (int j = 0; j < d_weight; ++j) {
			if (j > 128) {
				int xx = 0;
				xx++;
			}

			//destination x -> source x
			s_x = (float)j / (float)(d_weight - 1) * (float)(s_weight - 1);
			
			is_x = floor(s_x);
			is_x = is_x >= s_weight-1? s_weight - 2:is_x;
			
			temp = cal_bicubic(coeff + (is_y*s_weight + is_x) * 16, s_x-(float)is_x, float(is_y) - s_y);
			output[i*d_weight + j] = (char)temp;
		}
	}

	free(f_data);
	free(coeff);

}




