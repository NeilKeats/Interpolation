#include "interpolation.h"
#include<iostream>  
#include<malloc.h>  
#include<stdlib.h>  
#include<stdio.h>  
#include<string.h> 
#include<math.h>
#include<time.h>
//#include<cmath>

#define TIME_COUNT

const float M_PI=3.14159265358979f;

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

inline float bubic_conv_kernel(float a, float x) {
	float tmp = x < 0 ? -x : x;
	if (tmp <= 1)
		return ((a + 2)*pow(tmp, 3) - (a + 3)*pow(tmp, 2) + 1);
	else if (1 < tmp < 2)
		return (a*pow(tmp, 3) - 5 * a*pow(tmp, 2) + 8 * a*tmp - 4 * a);
	else
		return 0;
};

inline float lanczos_conv_kernel(float a, float x) {

	if (x == 0)
		return 1;
	else if (x<a && x > -a) {
		float result = (a*sinf(M_PI*x)*sinf(M_PI*x/a))/(M_PI*M_PI*x*x);
		return result;
	}
	else
		return 0;
}

void gradientX(const float *f_data, float * dx, float *dbuffer, DWORD width, DWORD height) {
	const DWORD local_width = width + 4;
	for(int i= 0; i<height+4;++i)
		for (int j = 2; j < width + 2; ++j) {
			dbuffer[local_width * i + j] =	 
				2.0 / 3 * f_data[local_width * i + j + 1] 
				- 1.0 / 12 * f_data[local_width * i + j + 2]
				- 2.0 / 3 * f_data[local_width * i + j -1]
				+ 1.0 / 12 * f_data[local_width * i + j -2 ];
		}	

	for(int i=0; i<height; ++i)
		for (int j = 0; j < width; ++j) {
			int row = i + 2;
			int col = j + 2;
			dx[i*width + j] = dbuffer[row*local_width+col];
		}
}

void gradientY(const float *f_data, float * dy, float *dbuffer, DWORD width, DWORD height) {
	const DWORD local_width = width + 4;
	for( int i =2 ; i< height+2 ; ++i)
		for (int j = 0; j < width + 4 ; ++j) {
			dbuffer[local_width * i + j] =   
				2.0 / 3 * f_data[local_width * (i - 1 ) + j ]
				- 1.0 / 12 * f_data[local_width * (i - 2 ) + j ]
				- 2.0 / 3 * f_data[local_width * (i + 1 ) + j ]
				+ 1.0 / 12 * f_data[local_width * (i + 2 ) + j ];
		}

	for (int i = 0; i<height; ++i)
		for (int j = 0; j < width; ++j) {
			int row = i + 2;
			int col = j + 2;
			dy[i*width + j] = dbuffer[row*local_width + col];
		}
}

void bicubic_coeff(const float *f_data, float * coeff, DWORD s_width, DWORD s_height){
	
	float *dx, *dy, *dxy,*d_buffx, *d_buffy;
	dx = (float*)malloc(sizeof(float)*s_width*s_height);
	dy = (float*)malloc(sizeof(float)*s_width*s_height);
	dxy = (float*)malloc(sizeof(float)*s_width*s_height);
	d_buffx = (float*)malloc(sizeof(float)*(s_width + 4)*(s_height + 4));
	d_buffy = (float*)malloc(sizeof(float)*(s_width + 4)*(s_height + 4));

	gradientY(f_data, dy, d_buffy, s_width, s_height);
	gradientX(f_data, dx, d_buffx, s_width, s_height);
	gradientY(d_buffx, dxy, d_buffy, s_width, s_height);

	float dTao[16], dAlpha[16];

	for(int i = 1 ; i<s_height ; ++i ){
		for (int j = 0; j < s_width-1 ; ++j) {
			dTao[0] = f_data[(i + 2)*(s_width + 4) + j + 2];
			dTao[1] = f_data[(i + 2)*(s_width + 4) + j + 3];
			dTao[2] = f_data[(i + 1)*(s_width + 4) + j + 2];
			dTao[3] = f_data[(i + 1)*(s_width + 4) + j + 3];

			
			dTao[4] = dx[(i)*s_width  + j ];
			dTao[5] = dx[(i)*s_width  + j + 1];
			dTao[6] = dx[(i - 1)*s_width  + j ];
			dTao[7] = dx[(i - 1)*s_width  + j + 1 ];

			dTao[8] = dy[(i)*s_width + j ];
			dTao[9] = dy[(i)*s_width + j + 1];
			dTao[10] = dy[(i - 1)*s_width + j];
			dTao[11] = dy[(i - 1)*s_width + j + 1];
			
			dTao[12] = dxy[(i)*s_width + j ];
			dTao[13] = dxy[(i)*s_width + j + 1];
			dTao[14] = dxy[(i - 1)*s_width + j];
			dTao[15] = dxy[(i - 1)*s_width + j + 1];
			
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

			float *fBicubic = coeff + ((i*s_width) + j) * 16;
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

void bicubic_spline_coeff(const float *f_data, float * coeff, DWORD s_width, DWORD s_height) {

	int i, j, k, l, m, n;
	float  Omiga[4][4], Beta[4][4];

	for (i = 0; i < s_height; ++i) {
		for (j = 0; j < s_width; ++j) {
			float *local_coeff = coeff + (i*s_width + j) * 16;

			//store neighbour value into Omeiga
			for (k = 0; k < 4; ++k)
				for (l = 0; l < 4; ++l)
					//reflect to original data (i + 1 - k, j - 1 + l)
					//reflect to filled data (i + 1 -k +2 , j - 1 +l +2)
					Omiga[k][l] = f_data[(i + 1 - k + 2)*(s_width + 4) + (j - 1 + l + 2)];

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
					local_coeff[k * 4 + l] = 0;
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
					local_coeff[k * 4 + l] = local_coeff[(3 - k) * 4 + (3 - l)];
					local_coeff[(3 - k) * 4 + (3 - l)] = m_dTemp;
				}
			}

			//for (int i = 4; i < 16; ++i)
			//	local_coeff[i] = 0;

		}

	}

}

float cal_bicubic(float *coeff, float s_x, float s_y, DWORD s_width, DWORD s_height) {
	//coefficient array obtained from the left corner element.
	//for y, using ceil(),causes our image data start from left bottom
	//for x, using floor()

	//is_x  int s_x;  is_y  int s_y
	int is_y, is_x;
	float x, y;
	is_y = ceil(s_y);
	is_y = is_y <= 0 ? 1 : is_y;
	//is_y = is_y <= 0 ? 0 : is_y;
	y = (float)is_y - s_y;

	is_x = floor(s_x);
	if (s_x >= s_width -1 )
		is_x = is_x;
	is_x = is_x >= s_width - 1 ? s_width - 2 : is_x;
	//is_x = is_x >= s_width - 1 ? s_width - 1 : is_x;
	x = s_x - (float)is_x;

	float *local_coeff = coeff + (is_y*s_width + is_x) * 16;
	
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

float cal_bicubic_kernel(const float *f_data, float s_x, float s_y, DWORD s_width, DWORD s_height, float a, int MODE) {
	int is_y, is_x;
	float x, y;
	is_y = ceil(s_y);
	//is_y = is_y <= 0 ? 1 : is_y;
	is_y = is_y <= 0 ? 0 : is_y;
	y = (float)is_y - s_y;

	is_x = floor(s_x);
	//is_x = is_x >= s_width - 1 ? s_width - 2 : is_x;
	is_x = is_x >= s_width - 1 ? s_width - 1 : is_x;
	x = s_x - (float)is_x;

	float c_j_[4], c_i_[4];
	if (MODE == MODE_BC_KERNEL) {
		c_j_[0] = bubic_conv_kernel(a, (1 + x));
		c_j_[1] = bubic_conv_kernel(a, (x));
		c_j_[2] = bubic_conv_kernel(a, (1 - x));
		c_j_[3] = bubic_conv_kernel(a, (2 - x));

		c_i_[0] = bubic_conv_kernel(a, (1 + y));
		c_i_[1] = bubic_conv_kernel(a, (y));
		c_i_[2] = bubic_conv_kernel(a, (1 - y));
		c_i_[3] = bubic_conv_kernel(a, (2 - y));
	}
	else if (MODE == MODE_LANCZOS_KERNEL) {
		c_j_[0] = lanczos_conv_kernel(a, (1 + x));
		c_j_[1] = lanczos_conv_kernel(a, (x));
		c_j_[2] = lanczos_conv_kernel(a, (1 - x));
		c_j_[3] = lanczos_conv_kernel(a, (2 - x));

		c_i_[0] = lanczos_conv_kernel(a, (1 + y));
		c_i_[1] = lanczos_conv_kernel(a, (y));
		c_i_[2] = lanczos_conv_kernel(a, (1 - y));
		c_i_[3] = lanczos_conv_kernel(a, (2 - y));
	}
	float result = 0;

	float sum = 0;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			sum += c_i_[i]*c_j_[j];

	//store neighbour value into Omeiga
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j){
			//reflect to original data (is_y + 1 - i, is_x - 1 +l)
			//reflect to filled data (is_y + 1 - i + 2 , is_x - 1 + l +2)
			int yy = (is_y + 1 - i + 2);
			int xx = (is_x - 1 + j + 2);
			//result += f_data[(is_y + 1 - i + 2)*(s_width + 4) + (is_x - 1 + j + 2)]*c_i_[i]*c_j_[j];
			float coeff = c_i_[i] * c_j_[j];
			result += f_data[yy*(s_width + 4) + xx] * coeff;
		}
	return result;
}

float nearest_neighbour(const float *f_data, float s_x, float s_y, DWORD s_width, DWORD s_height) {
	int is_y, is_x;
	is_y = floor(s_y);
	is_y = s_y - (float)is_y >= 0.5 ? (is_y + 1) : is_y;
	//is_y = is_y <= 0 ? 1 : is_y;
	//is_y = is_y <= 0 ? 0 : is_y;

	is_x = floor(s_x);
	is_x = s_x - (float)is_x >= 0.5 ? (is_x + 1) : is_x;
	//is_x = is_x >= s_width - 1 ? s_width - 2 : is_x;
	//is_x = is_x >= s_width - 1 ? s_width - 1 : is_x;

	int location = (is_y+2)*(s_width+4) + is_x + 2;

	return f_data[location];
}

float bilinear(const float* f_data, float s_x, float s_y, DWORD s_width, DWORD s_height) {
	int y_0, x_0, y_1, x_1;
	float f_y0_x0, f_y0_x1, f_y1_x0, f_y1_x1;
	float f_y0, f_y1, f_inter;

	y_0 = floor(s_y);
	y_1 = y_0 + 1;

	x_0 = floor(s_x);
	x_1 = x_0 + 1;

	f_y0_x0 = f_data[(y_0+2)*(s_width + 4) + (x_0 + 2)];
	f_y0_x1 = f_data[(y_0+2)*(s_width + 4) + (x_1 + 2)];
	f_y1_x0 = f_data[(y_1+2)*(s_width + 4) + (x_0 + 2)];
	f_y1_x1 = f_data[(y_1+2)*(s_width + 4) + (x_1 + 2)];

	f_y0 = f_y0_x0 + (s_x - (float)x_0)*(f_y0_x1 - f_y0_x0);
	f_y1 = f_y1_x0 + (s_x - (float)x_0)*(f_y1_x1 - f_y1_x0);
	f_inter = f_y0 + (s_y - (float)y_0)*(f_y1 - f_y0);

	return f_inter;
}

void fill_data(const uint8_t* s_data, float * f_data, DWORD s_width, DWORD s_height) {
//将原 w*h 大小的图像 填充到 (w+4) * (h+4)中，即上下左右增加了两行数据。
//原边界外的元素的值用近邻元素的像素值代替
	DWORD f_width = s_width + 4;
	DWORD f_height = s_height +4 ;

	int s_i, s_j;
	for(int i=0; i < f_height ;++i)
		for (int j = 0; j < f_width; ++j) {
			
			s_i = (i - 2) < 0 ? 0 : (i - 2);
			s_i = (i - 2) >= (int)s_height ? s_height-1 : s_i;

			s_j = (j - 2) < 0 ? 0 : (j - 2);
			s_j = (j - 2) >= (int)s_width ? s_width - 1 : s_j;
			f_data[i*f_width+j] = (float)s_data[s_i * s_width + s_j];
			
			/*
			if (i < 2 || (i + 2) >= s_height || j < 2 || j + 2 >= s_width)
				f_data[i*f_width + j] = 0;
			else
				f_data[i*f_width + j] = (float)s_data[(i-2) * s_width + (j-2)];
			*/
		}
}

void interpolation(const uint8_t *s_data, uint8_t *d_data, DWORD s_width, DWORD s_height, float width_scale, float height_scale, int MODE) {
	//prefix: d for destination  image; s for source image
	//
	clock_t start, end;
	double duration;
	start = clock();

	if (width_scale <= 0 || height_scale <= 0)
		return;
	DWORD d_width = s_width * width_scale;
	DWORD d_height = s_height * height_scale;

	//buffer with filling
	float *f_data = nullptr;
	f_data = (float *)malloc(sizeof(float)*(s_width + 4)*(s_height + 4));
	fill_data(s_data, f_data, s_width, s_height);

	//coefficient
	float *coeff=nullptr;
	float a;

	//compute coeef
	if (MODE == MODE_BICUBIC) {
		coeff = (float *)malloc(sizeof(float)*s_width*s_height * 16);
		bicubic_coeff(f_data, coeff, s_width, s_height);
	}
	else if (MODE == MODE_SPLINE) {
		coeff = (float *)malloc(sizeof(float)*s_width*s_height * 16);
		bicubic_spline_coeff(f_data, coeff, s_width, s_height);
	}
	else if (MODE == MODE_BC_KERNEL)
		a = -0.75;
	else if (MODE == MODE_LANCZOS_KERNEL)
		//using a = 3, casues some kind of stripes?. unknown
		a = 2;
	else if (MODE == MODE_NEAREST_NEIGHBOUR)
		;
	else if (MODE == MODE_BILINEAR)
		;

	uint8_t *output = d_data;
	float s_x, s_y;
	float temp;
	for (int i = 0; i < d_height; ++i) {
		//destination y ->  source y
		s_y = (float)i / (float)(d_height - 1) * (float)(s_height - 1);
		for (int j = 0; j < d_width; ++j) {
			//destination x -> source x
			s_x = (float)j / (float)(d_width - 1) * (float)(s_width - 1);

			if (MODE == MODE_BICUBIC || MODE == MODE_SPLINE)
				temp = cal_bicubic(coeff, s_x, s_y, s_width, s_height);
			else if (MODE == MODE_BC_KERNEL || MODE == MODE_LANCZOS_KERNEL)
				temp = cal_bicubic_kernel(f_data, s_x, s_y, s_width, s_height, a, MODE);
			else if (MODE == MODE_NEAREST_NEIGHBOUR)
				temp = nearest_neighbour(f_data, s_x, s_y, s_width, s_height);
			else if (MODE == MODE_BILINEAR)
				temp = bilinear(f_data, s_x, s_y, s_width, s_height);

			if (temp >= 255)
				temp = 255.0;
			else if (temp <= 0)
				temp = 0;
			output[i*d_width + j] = (uint8_t)temp;
		}
	}

	free(f_data);
	if(coeff!=nullptr)
		free(coeff);
	end = clock();
	duration = end - start;
	duration = duration / CLOCKS_PER_SEC;
	std::cout	<< "Time of scaling picture: "
				<<s_width<<"x"<<s_height 
				<<" in rate:"<<width_scale<<"x"<<height_scale
				<<" using method:"<<MODE_NAME[MODE]
				<<"	finished in:"<<duration<<"s"<<std::endl;
}




