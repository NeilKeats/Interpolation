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
#define WIDTH_SCALERATE 19.98020833
#define HEIGHT_SCALERATE 1.0

const string postfix = ".bmp";

const string inputDir =
	"D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET\\";
	//"D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\INPUT\\";
const string outputDir =
	"D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_OUTPUT\\";
const string inFileName =
	"epoxy_tensile_g_0";
	//"lena512";
const string outFileName =	"epoxy_tensile_g_0_cut_scale_";

void converte_gray_batch() {
	const string i_dir = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_INPUT\\";
	const string o_dir = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_OUTPUT\\";
	const string i_fname = "epoxy_tensile_";
	const string o_fname = "epoxy_tensile_g_";
	const string order[20] = {"0","1","2","3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9",
							"10","11","12","13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19"};
	bmp_i *p_image = nullptr;
	for (int i = 0; i < 20; ++i) {
		string filename = i_dir + i_fname + order[i] + postfix;
		string o_filename = o_dir + o_fname + order[i] + postfix;
		p_image = bmp_file_read(filename.c_str());
		if (p_image == nullptr)
			exit(1);
		p_image->converte_to_grey();
		bmp_file_write(p_image,o_filename.c_str());

		if(p_image != nullptr)
			delete p_image;
		p_image = nullptr;
	}
}

void cross_width_merge(bmp_i *d_image, const bmp_i *s_image, uint8_t offset, uint8_t width_scale) {
	//check width * height
	if (d_image->h != s_image->h || d_image->w != s_image->w*width_scale)
		return;
	DWORD d_width = d_image->w;
	DWORD s_width = s_image->w;

	if (offset == 0)
		for (int i = 0; i < s_image->h; ++i) {
			for (int j = 0; j < s_width ; ++j) {
				d_image->buf[i*d_width + j*width_scale + offset]
					= s_image->buf[i*s_width + j ];
			}
		}

	if(offset!= 0 )
		for (int i = 0; i < s_image->h; ++i) {
			for (int j = 0; j < s_width -1 ; ++j) {
				d_image->buf[i*d_width + j*width_scale + (20 - offset)]
					= s_image->buf[i*s_width + j + 1];
			}
			//d_image->buf[i*d_width + (s_width - 1)*width_scale + offset] = 255;
		}
	/*
	if(offset == 20)
		for (int i = 0; i<s_image->h; ++i)
			for (int j = 0; j <s_width-1; ++j) {
				d_image->buf[i*d_width + j*width_scale + offset]
					= s_image->buf[i*s_width + j];
			}
	*/
}

void merge_image() {
	const string i_dir = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_INPUT\\";
	const string o_dir = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_OUTPUT\\";
	const string i_fname = "epoxy_tensile_g_";
	const string o_fname = "epoxy_tensile_merge";
	const string order[21] = { "0","1","2","3" ,"4" ,"5" ,"6" ,"7" ,"8" ,"9",
		"10","11","12","13" ,"14" ,"15" ,"16" ,"17" ,"18" ,"19", "20" };
	bmp_i *p_image = nullptr;
	bmp_i *m_image = nullptr;
	m_image = bmp_file_read("D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_INPUT\\epoxy_tensile_g_0.bmp");
	if (m_image == nullptr)
		return;
	m_image->resize(20,1);

	for (int i = 0; i < 20; ++i) {
		string filename = i_dir + i_fname + order[i] + postfix;
		p_image = bmp_file_read(filename.c_str());
		if (p_image == nullptr)
			exit(1);

		cross_width_merge(m_image,p_image,i,20);

		delete p_image;
		p_image = nullptr;
	}

	m_image->cutimage(m_image->w - 19 ,m_image->h,0,0);
	string o_filename = o_dir + o_fname + postfix;
	bmp_file_write(m_image, o_filename.c_str());
	delete m_image;
	return;
}

void test(const bmp_i *inImage) {

	const string ref_filename = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_OUTPUT\\epoxy_tensile_merge.bmp";
	bmp_i *RefImage = nullptr;
	RefImage = bmp_file_read(ref_filename.c_str());
	if (RefImage == nullptr)
		return;

	for (int i = 0; i<5; ++i){
		float t_aver, r_aver, bias;
		bmp_i OutImage(inImage);
		OutImage.resize(WIDTH_SCALERATE, HEIGHT_SCALERATE);

		interpolation(
			inImage->buf, OutImage.buf,
			inImage->w, inImage->h,
			WIDTH_SCALERATE, HEIGHT_SCALERATE, A_MODE[i]);

		bias = bmp_compare(&OutImage,RefImage,&t_aver,&r_aver);
		
		/*
		float bias_ = (r_aver - t_aver);
		bias_ = bias_ < 0 ? -bias_ : bias_;
		*/

		cout << "target_aver:" << t_aver
			<< "	ref_aver:" << r_aver
			<< "	bias:" << bias << endl;

		string outputfile = outputDir + outFileName + MODE_NAME[i] + postfix;
		//OutImage.cutimage(OutImage.w - 20, OutImage.h, 0, 0);
		
		//bmp_file_write(&OutImage, outputfile.c_str());
	}

	delete RefImage;

};

void main()
{
	string inputfile = inputDir + inFileName + postfix;

	//string inputfile = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_OUTPUT\\epoxy_tensile_merge.bmp";

	bmp_i *inImage = nullptr;
	inImage = bmp_file_read(inputfile.c_str());
	if (inImage == nullptr)
		exit(1);

	/*
	inImage->cutimage(inImage->w-1,inImage->h,0,0);
	string outf = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_INPUT\\epoxy_tensile_g_0_cut.bmp";
	bmp_file_write(inImage, outf.c_str());
	*/

	//converte to gray
	/*
	inImage->converte_to_grey();
	string outf = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_INPUT\\epoxy_tensile_g_20.bmp";
	bmp_file_write(inImage,outf.c_str());
	return;
	*/


	/*
	const string outputName = "D:\\Codes\\VS\\CXX\\Interpolation\\Interpolation\\DATASET_OUTPUT\\op.bmp";
	bmp_i OutImage(inImage);
	float downscale = 0.0500495281789;
	OutImage.resize(downscale, 1.0);
	interpolation(	inImage->buf, OutImage.buf,
					inImage->w, inImage->h, 
					downscale, 1.0, MODE_BICUBIC);
	bmp_file_write(&OutImage, outputName.c_str());
	*/
	
	test(inImage);

	//converte_gray_batch();

	//merge_image();

	if(inImage != nullptr)
		delete inImage;
	system("pause");

}