#include <iostream>
#include <iomanip>
#include <queue>
#include <string>
#include <algorithm>
#include <Magick++.h>
#include <tuple>
#include <assert.h>
#include <cstdlib>

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
				Magick::Color c = *(pixel_cache + i*width + j);
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

class screen_row {
protected:
	int size;
	int left, right;
	int column;
	int row;
	mapped_image& image;

	point3d d0, d1;
	double distance = -1.;
public:
	std::tuple<point3d, point3d, double> get_diameter(){
		return std::make_tuple(d0, d1, distance);
	}

	screen_row(mapped_image& map, int c, int r, int rad): image(map), column(c), row(r){
		right = std::min(rad+1, image.get_width()-c);
		left = std::max(-rad, 0-c);
		size = right - left;
	}

	int width(){ return size; }
	point3d operator[] (int index){
	//	std::cerr << column + left << " " << column + left + size<< std::endl;;
		return image(column + left + index, row);
	}

	friend void try_update_diameter(screen_row& row, screen_row& new_row){
		for (int i=0; i<row.width(); i++)
			for (int j=0; j<new_row.width(); j++){
				point3d a = row[i], b = new_row[j];
				if (a.r > b.r) 
					std::swap(a, b);
				double new_dist = dist2(a, b);
				if (new_dist > row.distance)
					row.d0 = a, row.d1 = b, row.distance = new_dist;	
			}
		return; 
	}
};

typedef struct {
	point3d f;
	point3d s;
} zip;

class kernel {
protected:
	mapped_image& image;
	std::list<screen_row> screen;
	int column, row, radius;

	point3d d0, d1;
	double diameter;
public:
	kernel(mapped_image& img, int rad, int c, int r): 
		image(img), radius(rad)
	{
		this->move_to(c, r);

		for (auto it0 = screen.begin(); it0 != screen.end(); it0++)
			for (auto it1 = it0; it1 != screen.end(); it1++)
				try_update_diameter(*it1, *it0);
		diameter = 1.0;
		for (auto it = screen.begin(); it != screen.end(); it++){
			auto t = (*it).get_diameter();
			if (std::get<2>(t) > this->diameter)
				this->d0 = std::get<0>(t), this->d1 = std::get<1>(t), diameter = std::get<2>(t);
		}
	}
	
	void move_to(int c, int r) {
		if (!screen.empty())
			screen.clear();
		this->column = c, this->row = r;
		for (int i=std::max(-radius, 0-row); i<std::min(radius+1, image.get_height()-row); i++)
			screen.push_front(screen_row(image, column, row+i, radius));
	}

	void move_down() {
		diameter = -1.;
		screen.pop_back();
		row++;
		
		if (row + radius >= image.get_height()) {
			for (auto it = screen.begin(); it != screen.end(); it++){
				auto t = (*it).get_diameter();
				if (std::get<2>(t) > this->diameter)
					this->d0 = std::get<0>(t), this->d1 = std::get<1>(t), diameter = std::get<2>(t);
			}
			return;
		}

		screen.push_front(screen_row(image, column, row+radius-1, radius)); 

		screen_row& new_row = screen.front();
		for (auto it = screen.begin(); it != screen.end(); it++){
			try_update_diameter(*it, new_row);
			auto t = (*it).get_diameter();
			if (std::get<2>(t) > this->diameter)
				this->d0 = std::get<0>(t), this->d1 = std::get<1>(t), diameter = std::get<2>(t);
		}

	}

	std::tuple<point3d, point3d, double> get_diameter(){
		return std::make_tuple(d0, d1, diameter);
	}

	std::vector<zip> to_zipped(){
		std::vector<zip> v;
		for (auto it = screen.begin(); it != screen.end(); it++)
			for (int i=0; i<(*it).width(); i++)
				v.push_back({ (*it)[i], (*it)[i] });
		return v;
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

	 maxr = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;
	 minr = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;

	 maxg = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;
	 ming = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;

	 maxb = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;
	 minb = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;

//	std::cerr << minr << " " << ming << " " << minb << std::endl;
	assert(minr >= -0.001 && ming >= -0.001 && minb >= -0.001);
	assert(maxr <= 1.001 && maxg <= 1.001 && maxb <= 1.001);
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

//	       tr = d.r * cos(xz_angle) - d.b * sin(xz_angle);
//	double tb = d.r * sin(xz_angle) + d.b * cos(xz_angle);
//	d = { tr, d.g, tb };

//	std::cerr << d1.r << " " << d1.g << " " << d1.b << " " << xy_angle << " " << xz_angle << std::endl;
//	assert(abs(d.g) < 0.0001 && abs(d.b) < 0.0001);

	affine_scale(v);
	point3d m = {0.5, 0.5, 0.5};	

	/*
	double maxr = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;
	double minr = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.r<b.f.r;})->f.r;

	double maxg = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;
	double ming = std::min_element(v.begin(), v.end(), [](zip a, zip b){return a.f.g<b.f.g;})->f.g;

	double maxb = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;
	double minb = std::max_element(v.begin(), v.end(), [](zip a, zip b){return a.f.b<b.f.b;})->f.b;

	point3d m = { (maxr+minr)/2., (maxg+ming)/2., (maxb+minb)/2. };
	*/
	std::sort(v.begin(), v.end(), [&](zip a, zip b){return dist2(a.f, m) < dist2(b.f, m);});
}

Magick::Image petunin_filter(mapped_image& image, int radius){
	Magick::Image filtered_image(Magick::Geometry(image.get_width(), image.get_height()), "white");	

	int column = 0, row = 0;
	kernel petunin_kernel(image, radius, column, 0);

	while (column < image.get_width()){
		while (row < image.get_height()){
	//		std::cerr << row << std::endl;
			point3d d0, d1;
			auto t = petunin_kernel.get_diameter();
			d0 = std::get<0>(t), d1 = std::get<1>(t);
	//		std::cerr << std::setprecision(6) << std::fixed << std::get<2>(t) << std::endl;
			
			std::vector<zip> v = petunin_kernel.to_zipped();
			petunin_order(v, d0, d1);

			point3d mid = v[0].s;
			int r = mid.r, g = mid.g, b = mid.b;
			filtered_image.pixelColor(column, row, Magick::Color(r, g, b));

			petunin_kernel.move_down();
			row++;
		}

		row = 0;
		column++;
		petunin_kernel.move_to(column, row);
	}

	return filtered_image;
}

int main(int argc, char* argv[]){
#ifdef IMG_DEBUG
	mapped_image image("a.bmp");
	Magick::Image filtered_image = petunin_filter(image, 7);
	filtered_image.magick("PNG");
	filtered_image.write("filtered_a");
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
