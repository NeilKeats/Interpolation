#pragma once
#ifndef BMP_Image
#define BMP_Image

#include<iostream>  
#include<fstream>
#include<Windows.h>  
#include<malloc.h>  
#include<stdlib.h>  
#include<stdio.h>  
#include<string.h>  
#include<inttypes.h>
#include<string>

class bmp_i {
public:
	uint8_t *buf;                                //定义文件读取缓冲区  
	uint8_t *table;
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
	void resize(float scale);
	void write_buffer(const char *filename);
};

bmp_i::bmp_i(FILE *fp) {

	fread(&bf, sizeof(BITMAPFILEHEADER), 1, fp);//读取BMP文件头文件  
	fread(&bi, sizeof(BITMAPINFOHEADER), 1, fp);//读取BMP文件头文件信息  
	w = bi.biWidth;                            //获取图像的宽  
	h = bi.biHeight;                           //获取图像的高  
	bitSize = bi.biSizeImage;                  //获取图像的size  
	//buf = (char*)malloc(w*h * 3);                //分配缓冲区大小  
	buf = (uint8_t*)malloc(w*h);                //分配缓冲区大小  
	int xxx = bf.bfOffBits - (sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));
	int xxx_0 = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	table = (uint8_t*)malloc(xxx);
	fseek(fp, xxx_0, 0);
	fread(table, xxx, 1, fp);
	fseek(fp, bf.bfOffBits, 0);//定位到像素起始位置  
	//fread(buf, 1, w*h * 3, fp);                   //开始读取数据  
	int aligned_width = (w+3)/4*4;
	for (int i = 0; i<h; ++i){
		fseek(fp, bf.bfOffBits+aligned_width*i, 0);
		fread(buf+w*i, 1, w, fp);                   //开始读取数据  
	}
}

void bmp_i::write_image(FILE *fpw) {

	if (buf == nullptr || fpw == nullptr)
		return;
	FILE *fp = fpw;
	fwrite(&(this->bf), sizeof(BITMAPFILEHEADER), 1, fp);  //写入文件头  
	fwrite(&(this->bi), sizeof(BITMAPINFOHEADER), 1, fp);  //写入文件头信息  
	int xxx = bf.bfOffBits - (sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));

	RGBQUAD *board = new RGBQUAD[256];
	for (int i = 0; i < 256; ++i) {
		board[i].rgbBlue = i;
		board[i].rgbGreen = i;
		board[i].rgbRed = i;
		board[i].rgbReserved = 0;
	}

	fwrite(board, sizeof(RGBQUAD),256,fp);

	//fseek(fp, long(sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)), 0);
	uint8_t *p = this->buf;
	int aligned_width = (w + 3) / 4 * 4;
	int i, j;
	uint8_t empty = 0;
	for ( i = 0; i<h; i++)
	{
		fseek(fp, bf.bfOffBits + aligned_width*i, 0);
		for ( j = 0; j<w; j++)
		{
			fwrite(p++, 1, 1, fp);
		}
		for (; j<aligned_width; j++)
			fwrite(&empty, 1, 1, fp);

	}
};

void bmp_i::resize(float scale) {

	if (scale <= 0)
		return;
	this->bi.biWidth = this->bi.biWidth*scale;
	this->w = this->bi.biWidth;

	this->bi.biHeight = this->bi.biHeight*scale;
	this->h = this->bi.biHeight;

	this->bi.biSizeImage = this->bi.biSizeImage * scale * scale;
	this->bitSize = this->bi.biSizeImage;

	this->bf.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + this->bitSize;

	this->buf = (uint8_t*)malloc(w*h);
}

void bmp_i::write_buffer(const char *filename) {

	if (buf == nullptr)
		return;

	std::ofstream outfile;
	outfile.open(filename);

	uint8_t *p = this->buf;
	int i, j;
	std::string str_tmp;
	uint8_t empty = 0;
	for (i = 0; i<h; i++)
	{
		for (j = 0; j<w; j++)
		{
			str_tmp = std::to_string(*p++);
			outfile << str_tmp <<" ";
		}
		outfile << std::endl;
	}
	outfile.close();
}

#endif // !BMP_Image
