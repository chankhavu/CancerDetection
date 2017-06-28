#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include <cstring>
#include <queue>
#include <stack>
#include <algorithm>
#include <cstdlib>

#include "EasyBMP/EasyBMP.h"

const int MinWidth = 60;
const int MinHeight = 60;

const int MaxRemovalWidth = 40;
const int MaxRemovalHeight = 40;

typedef struct {
    int Top;
    int Left;
    int Width;
    int Height;

    int x;
    int y;

    int number;
} CellObject ;

void RemoveProtuberances(int** GrayImg, int Width, int Height){
	std::vector< std::pair<int,int> > tobe_removed;
	for (int i=0; i<Width; i++)
		for (int j=0; j<Height; j++){
			int count = 0;
			for (int di=-1; di<=1; di++)
				for (int dj=-1; dj<=1; dj++)
					if ((i+di) < 160 && (j+dj) < 160 && ((i+di) >= 0) && ((j+dj) >=0))
						if (GrayImg[i+di][j+dj] >= 0)
							count ++;
			if (count <= 5)
				tobe_removed.push_back(std::make_pair(i, j));
		}

	for (int i=0; i<tobe_removed.size(); i++)
		GrayImg[tobe_removed[i].first][tobe_removed[i].second] = -1;
}

void RemoveDarkHoles(int** GrayImg, int Width, int Height, int x, int y){
	if (GrayImg[x][y] != -1)
		return ;

	CellObject cell;
    cell.x = x;
    cell.y = y;

    cell.Top = y;
    cell.Left = x;
    cell.Width = 1;
    cell.Height = 1;
    cell.number = -3;

    int number = -3;

    std::vector< std::pair<int, int> > used(0);
    std::queue< std::pair<int, int> > que;

    used.push_back(std::make_pair(x, y));
    que.push(std::make_pair(x, y));
    GrayImg[x][y] = number;

    while (!que.empty()){
        std::pair<int ,int> p = que.front();
        que.pop();
        int i = p.first, j = p.second;

        if (j < cell.Top)
            cell.Top = j;
        if (j > cell.Top + cell.Height)
            cell.Height = j - cell.Top;

        if (i < cell.Left)
            cell.Left = i;
        if (i > cell.Left + cell.Width)
            cell.Width = i - cell.Left;


        for (int di = -1; di <= 1; di++)
            for (int dj = -1; dj <= 1; dj++)
				if (dj * di == 0)
                if (i + di < Width && i + di > 0 && j + dj < Height && j + dj > 0)
                    if (GrayImg[i + di][j + dj] == -1){
                        que.push(std::make_pair(i + di, j + dj));
                        used.push_back(std::make_pair(i + di, j + dj));
                        GrayImg[i + di][j + dj] = number;
                    }
    }

    if (cell.Width < MaxRemovalWidth || cell.Height < MaxRemovalHeight){
        for (unsigned int i=0; i<used.size(); i++)
            GrayImg[used[i].first][used[i].second] = 1;
        return ;
    }

	for (unsigned int i=0; i<used.size(); i++)
		GrayImg[used[i].first][used[i].second] = -2;

}

void AddCellObject(int** GrayImg, int Width, int Height, std::vector<CellObject>& cells, int x, int y, int radius = 1){

    // if it's background or is already known cell
    if (GrayImg[x][y] != 0)
        return ;

    CellObject cell;
    cell.x = x;
    cell.y = y;

    cell.Top = y;
    cell.Left = x;
    cell.Width = 1;
    cell.Height = 1;
    cell.number = cells.size()+1;

    int number = cells.size()+1;

    std::vector< std::pair<int, int> > used(0);
    std::queue< std::pair<int, int> > que;

    used.push_back(std::make_pair(x, y));
    que.push(std::make_pair(x, y));
    GrayImg[x][y] = number;

    while (!que.empty()){
        std::pair<int ,int> p = que.front();
        que.pop();
        int i = p.first, j = p.second;

        if (j < cell.Top)
            cell.Top = j;
        if (j > cell.Top + cell.Height)
            cell.Height = j - cell.Top;

        if (i < cell.Left)
            cell.Left = i;
        if (i > cell.Left + cell.Width)
            cell.Width = i - cell.Left;


        for (int di = -radius; di <= radius; di++)
            for (int dj = -radius; dj <= radius; dj++)
                if (i + di < Width && i + di > 0 && j + dj < Height && j + dj > 0)
                    if (GrayImg[i + di][j + dj] == 0){
                        que.push(std::make_pair(i + di, j + dj));
                        used.push_back(std::make_pair(i + di, j + dj));
                        GrayImg[i + di][j + dj] = number;
                    }
    }

    if (cell.Width < MinWidth || cell.Height < MinHeight){
        for (unsigned int i=0; i<used.size(); i++)
            GrayImg[used[i].first][used[i].second] = -1;
        return ;
    }

    cell.number = cells.size() + 1;
    cells.push_back(cell);
    return ;
}

void CreateHistogram(int** GrayImg, int* Hist, int Width, int Height, int* MinValue, int* MaxValue){
    *MinValue = GrayImg[0][0];
    *MaxValue = GrayImg[0][0];

    for (int i=0; i<256; i++)
        Hist[i] = 0;

    for (int i=0; i<Width; i++)
        for (int j=0; j<Height; j++){
            Hist[ GrayImg[i][j] ]++;
            if (GrayImg[i][j] > *MaxValue)
                *MaxValue = GrayImg[i][j];
            if (GrayImg[i][j] < *MinValue)
                *MinValue = GrayImg[i][j];
        }
}

int OtsuThreshold(int** GrayImg, int Width, int Height){
    // First, create histogram
    int Histogram[256], MinValue, MaxValue;
    CreateHistogram(GrayImg, Histogram, Width, Height, &MinValue, &MaxValue);

    int m = 0, n = 0;
    for (int t = MinValue; t <= MaxValue; t++)
        m += (t - MinValue) * Histogram[t], n += Histogram[t];

    float MaxSigma = -1.0;
    int Threshold = 0;
    int Alpha1 = 0, Beta1 = 0;

    for (int t = MinValue; t <= MaxValue; t++){
        Alpha1 += (t - MinValue) * Histogram[t];
        Beta1 += Histogram[t];

        // possibility of 1st class, w2 = 1-w1
        float w1 = (float)Beta1 / (1. * n);

        // a = a1 - a2, where a1 and a2 are arithmetic mean for classes 1 and 2
        float a = (float)Alpha1 / (1. * Beta1) - (float)(m - Alpha1) / (1. * (n - Beta1));

        // now calculate sigma
        float sigma = w1 * (1.0 - w1) * a * a;

        if (sigma > MaxSigma){
            MaxSigma = sigma;
            Threshold = t;
        }
    }

    return Threshold;
}

double LogarithmicRelation(CellObject cell, int** GrayImg, int eps){
    // eps size in pixels, assuming that the cell can be fit in a square 1x1
    double count = 0.0;

    for (int i=cell.Left; i < cell.Left + cell.Width; i += eps)
        for (int j=cell.Top; j < cell.Top + cell.Height; j += eps){
            bool inc = false;
            for (int di = 0; di < eps; di++)
                for (int dj = 0; dj < eps; dj++)
                     if (i + di <= cell.Left + cell.Width && j + dj <= cell.Top + cell.Height)
                         if (GrayImg[i + di][j + dj] == cell.number)
                             inc = true;
            if (inc)
                count += 1.0;
        }

    return std::log(count);
}

std::pair<double, double> LinearRegression(std::vector< std::pair<double, double> >& data)
{
    //x^t * x
    double xtx[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
	

    for (unsigned int i = 0; i < data.size(); i++)
    {
        xtx[0][1] += data[i].first;
        xtx[0][0] += data[i].first * data[i].first;
    }
    xtx[1][0] = xtx[0][1];
    xtx[1][1] = data.size();

    //inverse
    double xtxInv[2][2] = { {0.0, 0.0}, {0.0, 0.0} };
    double d = 1.0/(xtx[0][0]*xtx[1][1] - xtx[1][0]*xtx[0][1]);
    xtxInv[0][0] = xtx[1][1]*d;
    xtxInv[0][1] = -xtx[0][1]*d;
    xtxInv[1][0] = -xtx[1][0]*d;
    xtxInv[1][1] = xtx[0][0]*d;

    //times x^t
    double* xtxInvxt[2];
    xtxInvxt[0] = new double[data.size()];
    xtxInvxt[1] = new double[data.size()];

    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < data.size(); j++)
        {
            xtxInvxt[i][j] = xtxInv[i][0]*data[j].first + xtxInv[i][1];
        }
    }

    //times y
    double theta[2] = { 0.0, 0.0 };
    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < data.size(); j++)
        {
            theta[i] += xtxInvxt[i][j]*data[j].second;
        }
    }

    delete[] xtxInvxt[0];
    delete[] xtxInvxt[1];

    return std::make_pair(theta[0], theta[1]);
}

double FractalDimension(CellObject cell, int** GrayImg, int step = 1){
    int l = std::max(cell.Width, cell.Height);

    std::vector< std::pair<double,double> > boxCount;


	int bound = std::min(100, l/8);
	if (l < 400)
		bound = l/7;

	if (l < 200)
		bound = l/4.2;	
	   for (int s = 1; s <= bound; s += step){
        double d = LogarithmicRelation(cell, GrayImg, s);
        double e = (double)s / (double)l;
        double r = -std::log(e);

        boxCount.push_back(std::make_pair(r, d));        
    }

	std::pair<double,double> reg = LinearRegression(boxCount);
    return reg.first;
}


int main(int argc, char* argv[])
{

    if (argc != 2)
        return 0;

	std::string InputImage = argv[1];

#ifdef BITMAP
	BMP input_imgtext;
	input_imgtext.ReadFromFile(argv[1]);

	int height = input_imgtext.TellHeight();
	int width = input_imgtext.TellWidth();

	int** GrayImg = new int* [width];
	for (int i=0; i<width; i++)
		GrayImg[i] = new int [height];

	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++){
			int r = (int) input_imgtext(i, j)->Red;
			int g = (int) input_imgtext(i, j)->Green;
			int b = (int) input_imgtext(i, j)->Blue;

			GrayImg[i][j] = (int) (r+g+b)/3;

		}

	int Threshold = OtsuThreshold(GrayImg, width, height);

	std::vector<CellObject> cells;

	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++)
			if (GrayImg[i][j] > Threshold)
				GrayImg[i][j] = -1;
			else
				GrayImg[i][j] = 0;

	RemoveProtuberances(GrayImg, width, height);
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++)
			if (GrayImg[i][j] == 0)
				AddCellObject(GrayImg, width, height, cells, i, j, 2);

	for (int i=0; i<159; i++)
		for (int j=0; j<159; j++)
			if (GrayImg[i][j] == -1)
				RemoveDarkHoles(GrayImg, width, height, i, j);
	
	std::cerr << cells.size() << std::endl;
	CellObject c = *std::max_element(cells.begin(), cells.end(), 
			[](CellObject& a, CellObject& b) -> bool 
				{ return std::max(a.Width, a.Height) < std::max(b.Width, b.Height); });

//#ifdef DEBUG
	BMP cellobj_img;
	cellobj_img.SetSize(width, height);
	std::cerr << "c number: " << c.number << std::endl;
	for (int i=0; i<width; i++)
		for (int j=0; j<height; j++){
			RGBApixel temp;
			if (GrayImg[i][j] == c.number){
				temp.Red = 255;
				temp.Green = 255;
				temp.Blue = 255;
			}
			else {
				temp.Red = 0;
				temp.Green = 0;
				temp.Blue = 0;
			}
			cellobj_img.SetPixel(i, j, temp);
		}
	cellobj_img.WriteToFile("debug.bmp");
//#endif			

	std::cout << FractalDimension(c, GrayImg, 1) << std::endl;

	for (int i=0; i<160; i++)
		delete[] GrayImg[i];
    delete[] GrayImg;


#else
	std::ifstream input_imgtext;
	input_imgtext.open(InputImage.c_str());

	int** GrayImg = new int* [160];
	for (int i=0; i<160; i++)
		GrayImg[i] = new int [160];

	for (int i=0; i<160; i++)
		for (int j=0; j<160; j++){
			double a;			
			input_imgtext >> a;
			GrayImg[i][j] = std::round(a);			
		}
	std::cerr << "processing " << argv[1] << std::endl;

	input_imgtext.close();

	int Threshold = OtsuThreshold(GrayImg, 159, 159);

	std::vector<CellObject> cells;

	for (int i=0; i<159; i++)
		for (int j=0; j<159; j++)
			if (GrayImg[i][j] > Threshold)
				GrayImg[i][j] = -1;
			else
				GrayImg[i][j] = 0;

	RemoveProtuberances(GrayImg, 159, 159);
	for (int i=0; i<159; i++)
		for (int j=0; j<159; j++)
			if (GrayImg[i][j] == 0)
				AddCellObject(GrayImg, 159, 159, cells, i, j);

	for (int i=0; i<159; i++)
		for (int j=0; j<159; j++)
			if (GrayImg[i][j] == -1)
				RemoveDarkHoles(GrayImg, 159, 159, i, j);

	CellObject c = *std::max_element(cells.begin(), cells.end(), 
			[](CellObject& a, CellObject& b) -> bool 
				{ return std::max(a.Width, a.Height) < std::max(b.Width, b.Height); });

	std::cout << FractalDimension(c, GrayImg, 1) << std::endl;

	for (int i=0; i<160; i++)
		delete[] GrayImg[i];
    delete[] GrayImg;
#endif

    return 0;
}

