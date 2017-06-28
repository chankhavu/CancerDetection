#include <iostream>
#include <vector>
#include <algorithm>
#include <Magick++.h>

class custom_order {
public:
	virtual void operator() (std::vector<std::tuple<int, int, int>>& pixels) = 0;	
};

static class filter {
private:
	void apply_kernel(unsigned int column, unsigned int row, unsigned int radius,
			Magick::Image& image, Magick::Image filtered_image, custom_order& order){
		unsigned int top = std::max(row - radius, 0);
		unsigned int bottom = std::min(row + radius, image.rows());
		unsigned int left = std::max(column - radius, 0);
		unsigned int right = std::min(column + radius, image.columns());

		std::vector<std::tuple<int, int, int>>	
	}
public:
	void operator ()(Magick::Image& image, custom_order& order){
		Magick::Image filtered_image(Magick::Geometry(image.columns(), image.rows()), "black");
	}
};

int main()
{
#ifdef BMP_DEBUG
	Magick::Image sample_bmp("a.bmp");
	std::cout << sample_bmp.columns() << " " << sample_bmp.rows() << std::endl;
#endif
	return 0;
}
