#include <Magick++.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <numeric>

typedef struct {
	double r, g, b;	
} point3d;

inline double dist2(point3d a, point3d b){
	return (a.r-b.r)*(a.r-b.r) + (a.g-b.g)*(a.g-b.g) + (a.b-b.b)*(a.b-b.b);
}

class mapped_image {
private:
	int height;
	int width;
	point3d **pixels;
public:
	mapped_image(const std::string& image_name){
		Magick::Image image(image_name);
		Magick::PixelPacket* pixel_cache = image.getPixels(0, 0, image.columns(), image.rows());
		this->width = image.columns(), this->height = image.rows();

		pixels = new point3d* [height];
		for (int i=0; i<height; i++)
			pixels[i] = new point3d [width];

		for (int i=0; i<height; i++)
			for (int j=0; j<width; j++){
				Magick::Color c = image.pixelColor(j, i); //*(pixel_cache + i*width + j);
				pixels[i][j] = {
					(double)c.redQuantum(), (double)c.greenQuantum(), (double)c.blueQuantum() };
			}
	}

	~mapped_image(){
		for (int i=0; i<height; i++)
			delete[] pixels[i];
		delete[] pixels;
	}

	int get_height(){ return height; }
	int get_width(){ return width; }
	point3d operator() (int column, int row){
		return pixels[row][column];
	}
};

Magick::Image marginal_filter(mapped_image& image, int radius){
	Magick::Image filtered_image(Magick::Geometry(image.get_width(), image.get_height()), "white");	

	int column = 0, row = 0;
	
	while (column < image.get_width()){
		while (row < image.get_height()){
			point3d d0, d1;
			double distance = -1., diameter = -1.;
			
			std::vector<point3d> v;
			for (int i=std::max(row-radius, 0); i<std::min(row+radius+1, image.get_height()); i++)
				for (int j=std::max(column-radius, 0); j<std::min(column+radius+1, image.get_width()); j++)
					v.push_back(image(j, i));
			
			point3d center = v[v.size()/2+1];
			point3d average = {0, 0, 0};
			average = std::accumulate(v.begin(), v.end(), average, [&](point3d a, point3d b){ 
					return point3d({ a.r + b.r, a.g + b.g, a.b + b.b }); });
			average = { average.r / v.size(), average.g / v.size(), average.b / v.size() };

			std::vector<double> tr, tg, tb;
			for (int i=0; i<v.size(); i++)
				tr.push_back(v[i].r), tg.push_back(v[i].g), tb.push_back(v[i].b);
			std::sort(tr.begin(), tr.end());
			std::sort(tg.begin(), tg.end());
			std::sort(tb.begin(), tb.end());
			point3d median = { tr[tr.size()/2+1], tg[tg.size()/2+1], tb[tb.size()/2+1] };

			auto itr = std::min_element(tr.begin(), tr.end(), [&](double a, double b){ return 
					(abs(a-median.r)+abs(a-average.r)+abs(a-center.r)) <
					(abs(b-median.r)+abs(b-average.r)+abs(b-center.r)); });
			auto itg = std::min_element(tg.begin(), tg.end(), [&](double a, double b){ return 
					(abs(a-median.g)+abs(a-average.g)+abs(a-center.g)) < 
					(abs(b-median.g)+abs(b-average.g)+abs(b-center.g)); });
			auto itb = std::min_element(tb.begin(), tb.end(), [&](double a, double b){ return 
					(abs(a-median.b)+abs(a-average.b)+abs(a-center.b)) < 
					(abs(b-median.b)+abs(b-average.b)+abs(b-center.b)); });

			//filtered_image.pixelColor(column, row, 
			//		Magick::Color(tr[tr.size()/2+1], tg[tg.size()/2+1], tb[tb.size()/2+1], 0));
			filtered_image.pixelColor(column, row, Magick::Color((int)*itr, (int)*itg, (int)*itb));
			row++;
		}

		row = 0;
		column++;
	}

	return filtered_image;
}

int main(int argc, char* argv[]){
	
#ifdef IMG_DEBUG
	mapped_image image("a.bmp");
	Magick::Image filtered_image = marginal_filter(image, 4);
	filtered_image.magick("PNG");
	filtered_image.write("filtered_marginal");
	return 0;
#endif
	
	if (argc != 4)
		return 0;

	int kernel_size = std::atoi(argv[3]);

	mapped_image image(argv[1]);
	Magick::Image filtered_image = marginal_filter(image, 4);
	filtered_image.magick("PNG");
	filtered_image.write(argv[2]);
	return 0;

}
