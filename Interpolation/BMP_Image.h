#pragma once
#ifndef BMP_Image
#define BMP_Image

#include<iostream>  
#include<Windows.h>  
#include<malloc.h>  
#include<stdlib.h>  
#include<stdio.h>  
#include<string.h>  

class bmp_i {
public:
	char *buf;                                //定义文件读取缓冲区  
	char *table;
	//FILE *fp;                                 //定义文件指针  
	//FILE *fpw;                                //定义保存文件指针  
	DWORD w, h;                                //定义读取图像的长和宽  
	DWORD bitCorlorUsed;                      //定义  
	DWORD bitSize;                            //定义图像的大小  
	BITMAPFILEHEADER bf;                      //图像文件头  
	BITMAPINFOHEADER bi;                      //图像文件头信息  

	bmp_i() { buf = nullptr; };
	bmp_i(FILE *fp);
	~bmp_i() {
		free(buf);
		free(table);
	}
	void write_image(FILE *fpw);
	void resize(int scale);
};

bmp_i::bmp_i(FILE *fp) {

	fread(&bf, sizeof(BITMAPFILEHEADER), 1, fp);//读取BMP文件头文件  
	fread(&bi, sizeof(BITMAPINFOHEADER), 1, fp);//读取BMP文件头文件信息  
	w = bi.biWidth;                            //获取图像的宽  
	h = bi.biHeight;                           //获取图像的高  
	bitSize = bi.biSizeImage;                  //获取图像的size  
	//buf = (char*)malloc(w*h * 3);                //分配缓冲区大小  
	buf = (char*)malloc(w*h);                //分配缓冲区大小  
	int xxx = bf.bfOffBits - (sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));
	int xxx_0 = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	table = (char*)malloc(xxx);
	fseek(fp, xxx_0, 0);
	fread(table, xxx, 1, fp);
	fseek(fp, bf.bfOffBits, 0);//定位到像素起始位置  
	//fread(buf, 1, w*h * 3, fp);                   //开始读取数据  
	fread(buf, 1, w*h, fp);                   //开始读取数据  
}

void bmp_i::write_image(FILE *fpw) {
	if (buf == nullptr || fpw == nullptr)
		return;
	FILE *fp = fpw;
	fwrite(&(this->bf), sizeof(BITMAPFILEHEADER), 1, fp);  //写入文件头  
	fwrite(&(this->bi), sizeof(BITMAPINFOHEADER), 1, fp);  //写入文件头信息  
	int xxx = bf.bfOffBits - (sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));
	fwrite(this->table, xxx,1,fp);
	//fseek(fp, long(sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)), 0);
	char *p = this->buf;
	for (int j = 0; j<h; j++)
	{
		for (int i = 0; i<w; i++)
		{
			fwrite(p++, 1, 1, fp);
		}
	}
};

void bmp_i::resize(int scale) {

	if (scale < 1)
		return;
	this->bi.biWidth = this->bi.biWidth*scale;
	this->w = this->bi.biWidth;

	this->bi.biHeight = this->bi.biHeight*scale;
	this->h = this->bi.biHeight;

	this->bi.biSizeImage = this->bi.biSizeImage * scale * scale;
	this->bitSize = this->bi.biSizeImage;

	this->bf.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + this->bitSize;

	this->buf = (char*)malloc(w*h);
}

#endif // !BMP_Image
