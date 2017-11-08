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
	uint8_t *buf;                                //�����ļ���ȡ������  
	uint8_t *table;
	//FILE *fp;                                 //�����ļ�ָ��  
	//FILE *fpw;                                //���屣���ļ�ָ��  
	DWORD w, h;                                //�����ȡͼ��ĳ��Ϳ�  
	DWORD bitCorlorUsed;                      //����  
	DWORD bitSize;                            //����ͼ��Ĵ�С  
	BITMAPFILEHEADER bf;                      //ͼ���ļ�ͷ  
	BITMAPINFOHEADER bi;                      //ͼ���ļ�ͷ��Ϣ  

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

	fread(&bf, sizeof(BITMAPFILEHEADER), 1, fp);//��ȡBMP�ļ�ͷ�ļ�  
	fread(&bi, sizeof(BITMAPINFOHEADER), 1, fp);//��ȡBMP�ļ�ͷ�ļ���Ϣ  
	w = bi.biWidth;                            //��ȡͼ��Ŀ�  
	h = bi.biHeight;                           //��ȡͼ��ĸ�  
	bitSize = bi.biSizeImage;                  //��ȡͼ���size  
	//buf = (char*)malloc(w*h * 3);                //���仺������С  
	buf = (uint8_t*)malloc(w*h);                //���仺������С  
	int xxx = bf.bfOffBits - (sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));
	int xxx_0 = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	table = (uint8_t*)malloc(xxx);
	fseek(fp, xxx_0, 0);
	fread(table, xxx, 1, fp);
	fseek(fp, bf.bfOffBits, 0);//��λ��������ʼλ��  
	//fread(buf, 1, w*h * 3, fp);                   //��ʼ��ȡ����  
	int aligned_width = (w+3)/4*4;
	for (int i = 0; i<h; ++i){
		fseek(fp, bf.bfOffBits+aligned_width*i, 0);
		fread(buf+w*i, 1, w, fp);                   //��ʼ��ȡ����  
	}
}

void bmp_i::write_image(FILE *fpw) {

	if (buf == nullptr || fpw == nullptr)
		return;
	FILE *fp = fpw;
	fwrite(&(this->bf), sizeof(BITMAPFILEHEADER), 1, fp);  //д���ļ�ͷ  
	fwrite(&(this->bi), sizeof(BITMAPINFOHEADER), 1, fp);  //д���ļ�ͷ��Ϣ  
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
