#include<iostream>  
#include<Windows.h>  
#include<malloc.h>  
#include<stdlib.h>  
#include<stdio.h>  
#include<string.h>  
#include "BMP_Image.h"
#include "interpolation.h"
using namespace std;

#define SCALERATE 2.0
#define WIDTH_SCALERATE 4.0
#define HEIGHT_SCALERATE 1.0


const string inputDir = 
	"D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\INPUT\\";
const string outputDir =
	"D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\OUTPUT1\\";
const string inFileName = 
	"mr";
const string postfix = ".bmp";



void test(const bmp_i *inImage) {
	bmp_i OutImage(inImage);
	OutImage.resize(WIDTH_SCALERATE, HEIGHT_SCALERATE);

	for (int i = 0; i<5; ++i){
		interpolation(
			inImage->buf, OutImage.buf,
			inImage->w, inImage->h,
			WIDTH_SCALERATE, HEIGHT_SCALERATE, A_MODE[i]);
		string outputfile = outputDir + inFileName + string("_") + MODE_NAME[i] + postfix;
		bmp_file_write(&OutImage, outputfile.c_str());
	}
};

void main()
{

	FILE *fp = nullptr, *fpw = nullptr;
	
	string inputfile = inputDir + inFileName + postfix;

	bmp_i *inImage = nullptr;
	inImage = bmp_file_read(inputfile.c_str());
	if (inImage == nullptr)
		exit(1);

	//converte to gray
	/*
	inImage.converte_to_gray();
	FILE *fpw = fopen("D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\OUTPUT\\test.bmp", "wb");
	inImage.write_image(fpw);
	fclose(fpw);
	return;
	*/


	/*
	bmp_i OutImage(inImage);
	OutImage.resize(SCALERATE, SCALERATE);
	for(int i=0;i<5;++i)
	interpolation(	inImage->buf, OutImage.buf,
					inImage->w, inImage->h, 
					SCALERATE, SCALERATE, A_MODE[i]);
	bmp_file_write(&OutImage,outputName);
	*/
	test(inImage);

	delete inImage;
	system("pause");

}