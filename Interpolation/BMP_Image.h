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
	char *buf;                                //�����ļ���ȡ������  
	char *table;
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
	void resize(int scale);
};

bmp_i::bmp_i(FILE *fp) {

	fread(&bf, sizeof(BITMAPFILEHEADER), 1, fp);//��ȡBMP�ļ�ͷ�ļ�  
	fread(&bi, sizeof(BITMAPINFOHEADER), 1, fp);//��ȡBMP�ļ�ͷ�ļ���Ϣ  
	w = bi.biWidth;                            //��ȡͼ��Ŀ�  
	h = bi.biHeight;                           //��ȡͼ��ĸ�  
	bitSize = bi.biSizeImage;                  //��ȡͼ���size  
	//buf = (char*)malloc(w*h * 3);                //���仺������С  
	buf = (char*)malloc(w*h);                //���仺������С  
	int xxx = bf.bfOffBits - (sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER));
	int xxx_0 = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	table = (char*)malloc(xxx);
	fseek(fp, xxx_0, 0);
	fread(table, xxx, 1, fp);
	fseek(fp, bf.bfOffBits, 0);//��λ��������ʼλ��  
	//fread(buf, 1, w*h * 3, fp);                   //��ʼ��ȡ����  
	fread(buf, 1, w*h, fp);                   //��ʼ��ȡ����  
}

void bmp_i::write_image(FILE *fpw) {
	if (buf == nullptr || fpw == nullptr)
		return;
	FILE *fp = fpw;
	fwrite(&(this->bf), sizeof(BITMAPFILEHEADER), 1, fp);  //д���ļ�ͷ  
	fwrite(&(this->bi), sizeof(BITMAPINFOHEADER), 1, fp);  //д���ļ�ͷ��Ϣ  
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
