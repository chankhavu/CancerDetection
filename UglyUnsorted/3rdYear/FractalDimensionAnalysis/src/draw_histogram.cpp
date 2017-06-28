#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include "EasyBMP/EasyBMP.h"
#include "EasyBMP/EasyBMP_Font.h"

int count(const char* filename, int* hist, float a, int c){
	std::ifstream f;
	f.open(filename);
	int max = 0;

	for (int i=0; i<c; i++)
		hist[i] = 0;

	double x;
	while (f >> x){
		int t = hist[(int)trunc((x-a) * 100.0)]++;
		if (t > max)
			max = t;
	}
	return max;
}

void draw_col(BMP& result, int value, int max, int ind, int num, float a, int c, int cw){
	int y = (num+1) * (100 + 10) - 1;
	int x = ind * cw;

	int r, g, b;
	if (num == 0)
		r = 128, g = 128, b = 128;
	if (num == 1)
		r = 255, g = 0, b = 0;
	if (num == 2)
		r = 0, g = 255, b = 0;
	if (num == 3)
		r = 0, g = 0, b = 255;

	int h = (100*value)/max;
	for (int j = y; j > y-h; j--)
		for (int i = x; i < x+cw; i++){
			result(i,j)->Red = r;
			result(i,j)->Green = g;
			result(i,j)->Blue = b;
		}

	for (int j = y; j > y-100; j--){
		result(x, j)->Red = 150;
		result(x, j)->Green = 150;
		result(x, j)->Blue = 150;
	}

	if ((int)round((a + (float)ind*0.01)*100.0) % 5 == 0)
		for (int j = y; j > y-100; j--){
			result(x, j)->Red = 0;
			result(x, j)->Green = 0;
			result(x, j)->Blue = 0;
		}

	if ((int)round((a + (float)ind*0.01)*100.0) % 10 == 0){
		for (int j=y; j > y-110; j--){
			result(x,j)->Red = 0;
			result(x,j)->Green = 0;
			result(x,j)->Blue = 0;
		}

		char text_string[10];
		RGBApixel FontColor;
		FontColor.Red = 0;
		FontColor.Green = 0;
		FontColor.Blue = 0;
		sprintf(text_string, "%.1f", a+ind*0.01);
		PrintString(result, text_string, x+1, y-108, 6, FontColor);
	}

}

int main(int argc, char** argv){
	if (argc != 9){
		std::cerr << "Wrong number of arguments (need 8): " << std::endl <<
			"[red txt] [green txt] [blue txt] [gray txt] [result bitmap] [from float value] [to float value] [column width]" << std::endl;
		return 0;
	}

	float x, y;
	int cw;
	sscanf(argv[6], "%f", &x);
	sscanf(argv[7], "%f", &y);
	sscanf(argv[8], "%d", &cw);	
	x = trunc(x*100.0)/100.0;
	y = ceil(y*100.0)/100.0;

	BMP result;
	int c = (int)(ceil(y*100.0) - trunc(x*100.0));
	int w = (c+1)*cw;
	int h = 110 * 4;
	result.SetSize(w, h);
	result.SetBitDepth(24);	

	int* r = new int[c+1];
	int* g = new int[c+1];
	int* b = new int[c+1];
	int* s = new int[c+1];

	int mr = count(argv[1], r, x, c);
	int mg = count(argv[2], g, x, c);
	int mb = count(argv[3], b, x, c);
	int ms = count(argv[4], s, x, c);

	for (int i=0; i<c; i++){
		draw_col(result, s[i], ms, i, 0, x, c, cw);
		draw_col(result, r[i], mr, i, 1, x, c, cw);
		draw_col(result, g[i], mg, i, 2, x, c, cw);
		draw_col(result, b[i], mb, i, 3, x, c, cw);
	}

	result.WriteToFile(argv[5]);

	return 0;
}
