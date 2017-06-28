// Автор: Чан Ха Ву (tran.thecoder@gmail.com)
//
//
// Внимание: 
// Данный код тесно связан с публикацией "Непараметрический
// критерий эквивалентности генеральных совокупностей, основанных
// на мере близости между выборками" Д.А.Клюшин, Ю.И.Петунин. 
// 
// Переменные и функции в коде названы в соответствии с параграфом 3
// "Мера близости между выборками (Р-статистика)" вышеупомянутой
// публикации.

#ifndef __METRIC_H__
#define __METRIC_H__

#include <vector>
#include <algorithm>
#include <iostream>
#include <array>
#include <tuple>

using std::vector;
using std::cerr;
using std::tuple;
using std::make_tuple; // c++11

namespace P {

	const double G = 3.0;

	// возвращает тройку значении (h, p1, p2), где h - мера близости (или
	// Р-статистика) выборок x и z, а (p1, p2) - доверительный интервал для
	// h, построенный по правилу три сигма
	tuple<double, double, double> metric(vector<double>& x, vector<double>& z){

		vector<double> ordered_x(x);
		sort(ordered_x.begin(), ordered_x.end());

		// some lambdas with capturing (c++11) 
		auto calc_h = [&](size_t i, size_t j) -> double{
			if (i >= j)
				return 0.0;

			size_t m = z.size();
			size_t h = 0;
			for (size_t k=0; k<m; k++)
				if (ordered_x[i] < z[k] && z[k] < ordered_x[j])
					h++;
			return double(h) / double(m);
		};

		auto calc_p = [&](size_t i, size_t j) -> double {
			if (i >= j)
				return 0.0;

			size_t q = j - i;
			size_t n = x.size();
			
			return double(1.0 * q) / (n+1.0);			
		};

		// some lambdas without capturing (c++11)
		auto p1 = [](double h, unsigned int m) -> double {
			double _x = (double(h*m) + double(G*G)/2.0)/double(m + G*G);
			double si = G * std::sqrt(h*(1.0-h)*double(m) + double(G*G)/4.0)/double(m + G*G);
			return _x - si;
		};

		auto p2 = [](double h, unsigned int m) -> double {
			double _x = (double(h*m) + double(G*G)/2.0)/double(m + G*G);
			double si = G * std::sqrt(h*(1.0-h)*double(m) + double(G*G)/4.0)/double(m + G*G);
			return _x + si;
		};

		size_t n = x.size();
		size_t N = n*(n-1)/2;
		size_t L = 0;

		// parameters
		unsigned int m = z.size();

		for (size_t i=0; i<n-1; i++)
			for (size_t j=i+1; j<n; j++){
				double h = calc_h(i, j);
				double p = calc_p(i, j);

				if (p1(h, m) < p && p < p2(h, m))
					L++;
			}

		double h = double(L) / double(N);
		double p = p1(h, N);
		double P = p2(h, N);

		return make_tuple(h, p, P);
	}

}; // end of namespace P

#endif
