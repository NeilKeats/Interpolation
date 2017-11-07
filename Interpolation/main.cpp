#include<iostream>  
#include<Windows.h>  
#include<malloc.h>  
#include<stdlib.h>  
#include<stdio.h>  
#include<string.h>  
#include "BMP_Image.h"
#include "interpolation.h"
using namespace std;

#define SCALERATE 2
void main()
{
	//char fileName[48];
	//const char *fileName = "D:\\Projects\\test.bmp";
	//const char *fileName = "D:\\Projects\\data\\dic challenge12\\crop_oht_cfrp_0.bmp";
	const char *fileName = "D:\\Projects\\lena512.bmp";
	FILE *fp = nullptr;
	int r, g, b, pix;

	//cout << "请输入要打开文件的名字：";
	//cin >> fileName;
	if ((fp = fopen(fileName, "rb")) == NULL)
	{
		cout << "文件未找到！";
		exit(0);
	}

	bmp_i inImage(fp); 

	fclose(fp);

	//show image
	/*
	HWND wnd;                                 //窗口句柄  
	HDC dc;                                   //绘图设备环境句柄
	wnd = GetForegroundWindow();               //获取窗口句柄  
	dc = GetDC(wnd);                           //获取绘图设备  
	int x = 40;
	int y = 40;
	char *p = inImage.buf;
	for (int j = 0; j<inImage.h; j++)
	{
		for (int i = 0; i<inImage.w; i++)
		{
			//b = *p++; g = *p++; r = *p++;
			b = *p++;
			g = r = b;
			pix = RGB(r, g, b);
			SetPixel(dc, x + i, y + inImage.h - j, pix);
		}
	}
	*/

	bmp_i OutImage = inImage;
	OutImage.buf = nullptr;
	OutImage.resize(SCALERATE);
	char *tmp = OutImage.table;
	OutImage.table = (char* )malloc(OutImage.bf.bfOffBits-54);
	memcpy(OutImage.table,tmp, OutImage.bf.bfOffBits - 54);
	interpolation(inImage.buf, OutImage.buf, inImage.w, inImage.h, SCALERATE, SCALERATE, MODE_SPLINE);
	FILE *fpw = fopen("TestOutputx1_spline.bmp","wb");
	OutImage.write_image(fpw);
	//inImage.write_image(fpw);
	fclose(fpw);

	/*
	x = 40;
	y = 340;
	p = OutImage.buf;
	for (int j = 0; j<OutImage.h; j++)
	{
		for (int i = 0; i<OutImage.w; i++)
		{
			//b = *p++; g = *p++; r = *p++;  D:\Projects\test.bmp
			b = *p++;
			g = r = b;
			pix = RGB(r, g, b);
			SetPixel(dc, x + i, y + OutImage.h - j, pix);
		}
	}*/

	//system("pause");
		//write image


	//return fp;  
}