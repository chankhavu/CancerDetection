#include <Magick++.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <cstdlib>

typedef struct {
	double r, g, b;	
} point3d;

inline double dist2(point3d a, point3d b){
	return (a.r-b.r)*(a.r-b.r) + (a.g-b.g)*(a.g-b.g) + (a.b-b.b)*(a.b-b.b);
}
typedef struct {
	point3d f;
	point3d s;
} zip;

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


void affine_move(std::vector<zip>& v, point3d u){
	for (int i=0; i<v.size(); i++)
		v[i].f = { v[i].f.r + u.r, v[i].f.g + u.g, v[i].f.b + u.b };
}

void affine_rotationxy(std::vector<zip>& v, double xy_angle){
	for (int i=0; i<v.size(); i++){
		double r = v[i].f.r * cos(xy_angle) - v[i].f.g * sin(xy_angle);
		double g = v[i].f.r * sin(xy_angle) + v[i].f.g * cos(xy_angle);
		v[i].f = {r, g, v[i].f.b};
	}
}
void affine_rotationxz(std::vector<zip>& v, double xz_angle){
	for (int i=0; i<v.size(); i++){
		double r = v[i].f.r * cos(xz_angle) - v[i].f.b * sin(xz_angle);
		double b = v[i].f.r * sin(xz_angle) + v[i].f.b * cos(xz_angle);
		v[i].f = {r, v[i].f.g, b};
	}
}

void affine_scale(std::vector<zip>& v){
	double maxr = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;
	double minr = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;

	double maxg = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;
	double ming = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;

	double maxb = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;
	double minb = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;

	point3d u;
	if (minr < 0.)		
		affine_move(v, u = {-minr, 0, 0}), maxr -= minr;
	if (ming < 0.)
		affine_move(v, u = {0, -ming, 0}), maxg -= ming;
	if (minb < 0.)
		affine_move(v, u = {0, 0, -minb}), maxb -= minb;

	for (int i=0; i<v.size(); i++){
		v[i].f.r /= maxr;
		v[i].f.g /= maxg;
		v[i].f.b /= maxb;
	}
}

void petunin_order(std::vector<zip>& v, point3d d0, point3d d1){	
	std::cerr << std::fixed << std::setprecision(6) ;
	point3d u = { -d0.r, -d0.g, -d0.b };
	affine_move(v, u);

	point3d d = { d1.r - d0.r, d1.g - d0.g, d1.b - d0.b };
	double xy_angle = -atan2(d.g, d.r);
	affine_rotationxy(v, xy_angle);

	double tr = d.r * cos(xy_angle) - d.g * sin(xy_angle);
	double tg = d.r * sin(xy_angle) + d.g * cos(xy_angle);
	d = { tr, tg, d.b };

	double xz_angle = -atan2(d.b, d.r);
	affine_rotationxz(v, xz_angle);

	       tr = d.r * cos(xz_angle) - d.b * sin(xz_angle);
	double tb = d.r * sin(xz_angle) + d.b * cos(xz_angle);
	d = { tr, d.g, tb };

//	std::cerr << d1.r << " " << d1.g << " " << d1.b << " " << xy_angle << " " << xz_angle << std::endl;
	assert(abs(d.g) < 0.0001 && abs(d.b) < 0.0001);

//	affine_scale(v);
//	point3d m = {0.5, 0.5, 0.5};	

	
	double maxr = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;
	double minr = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;

	double maxg = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;
	double ming = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;

	double maxb = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;
	double minb = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;

	point3d m = { (maxr+minr)/2., (maxg+ming)/2., (maxb+minb)/2. };
	
	std::sort(v.begin(), v.end(), [&](zip a, zip b){
				//return dist2(a.f, m) < dist2(b.f, m);
				return (pow((a.f.r - m.r)/(maxr - minr), 2) + 
					   pow((a.f.g - m.g)/(maxg - ming), 2) + 
					   pow((a.f.b - m.b)/(maxb - minb), 2))  < 
					(pow((b.f.r - m.r)/(maxr - minr), 2) + 
					   pow((b.f.g - m.g)/(maxg - ming), 2) + 
					   pow((b.f.b - m.b)/(maxb - minb), 2));
			});
}



Magick::Image petunin_filter(mapped_image& image, int radius){
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

			for (int i=0; i<v.size(); i++)
				for (int j=i; j<v.size(); j++)
					if ((distance = dist2(v[i], v[i])) > diameter)
						d0 = v[i], d1 = v[j], diameter = distance;
				
			std::vector<zip> vv;
			for (int i=0; i<v.size(); i++)	
				vv.push_back({v[i], v[i]});
			petunin_order(vv, d0, d1);

			point3d mid = vv[0].s;
			int r = mid.r, g = mid.g, b = mid.b;
			filtered_image.pixelColor(column, row, Magick::Color(r, g, b, 0));

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
	Magick::Image filtered_image = petunin_filter(image, 1);
	filtered_image.magick("PNG");
	filtered_image.write("filtered_b");
#endif

	if (argc != 4)
		return 0;

	int kernel_size = std::atoi(argv[3]);

	mapped_image image(argv[1]);
	Magick::Image filtered_image = petunin_filter(image, kernel_size);
	filtered_image.magick("PNG");
	filtered_image.write(argv[2]);

	return 0;
}
