#include <iostream>
#include<fstream>
#include "Header.h"
//#include "Alt.h"
#include<math.h>
using namespace std;
void input_grid(double X[], double Y[], int& n_x, int& n_y);
int main()
{
	int n_skv;
	double skvminX1, skvminY1, skvmaxX1, skvmaxY1, skvminX2, skvminY2, skvmaxX2, skvmaxY2, r;
	double xmin, xmax, ymin, ymax, k;
	double ki_1, ki_2, phi, s2, tetta1, tetta2;
	int n_x, n_y;
	input_xy(&xmin, &xmax, &ymin, &ymax, &ki_1, &ki_2, &phi, &s2);
	input_wells(&n_skv, &skvminX1, &skvminY1, &skvmaxX1, &skvmaxY1, &skvminX2, &skvminY2, &skvmaxX2, &skvmaxY2, &tetta1, &tetta2);
	input_mesh(&n_x, &n_y, &k);
	double* x;
	x = new double[n_x];

	double* y;
	y = new double[n_y];


	double* resX1 = new double[2];
	double* resY1 = new double[2];
	double* resX2 = new double[2];
	double* resY2 = new double[2];

	double* skvgridX1 = new double[5000];
	double* skvgridY1 = new double[5000];
	double* skvgridX2 = new double[5000];
	double* skvgridY2 = new double[5000];
	double* mid = new double[5000];

	//создание равномерной сетки по х и по у
	createevengridXY(x, n_x, xmin, xmax);
	createevengridXY(y, n_y, ymin, ymax);
	outgrid(x, n_x, y, n_y);
	FindSkvXY(x, n_x, skvminX1, skvmaxX1, resX1);
	FindSkvXY(x, n_x, skvminX2, skvmaxX2, resX2);
	FindSkvXY(y, n_y, skvminY2, skvmaxY2, resY2);
	FindSkvXY(y, n_y, skvminY1, skvmaxY1, resY1);

	double* resX = new double[5000];
	double* resY = new double[5000];

	LocalElemTriangle(x, n_x, y, n_y);
	//LocalElem(x, n_x, y, n_y);
	BC1(x, n_x, y, n_y);
	
	BC2_Triangle(resX, n_x, resY, n_y, skvminX1, skvminY1, skvminX2, skvminY2, tetta1, tetta2);
	outgrid(resX, n_x, resY, n_y);
	out_mat(ki_1, ki_2, phi, s2, n_x, x);
}

void input_grid(double X[], double Y[], int& n_x, int& n_y)
{
	ifstream in("gridX.txt");
	in >> n_x ;

	//in >> n_x >> 0;
	for (int i = 0; i < n_x; i++)
	{
		double x;
		in >> x;
		X[i] = x;
	}
	in.close();
	in.open("gridY.txt");
	in >> n_y;
	for (int i = 0; i < n_y; i++)
	{
		double y;
		in >> y;
		Y[i] = y;
	}
	in.close();
}

void input_mesh(int* n_x, int* n_y, double* k)
{
	ifstream in("mesh_xy.txt");
	in >> *n_x >> *n_y >> *k;
	in.close();
}

void input_xy(double* xmin, double* xmax, double* ymin, double* ymax, double* ki_1, double* ki_2, double* phi, double* s2)
{
	ifstream in("domain_xy.txt");
	in >> *xmin >> *xmax >> *ymin >> *ymax >> *ki_1 >> *ki_2 >> *phi >> *s2;
	in.close();
}

void input_wells(int* n, double* skvminX1, double* skvminY1, double* skvmaxX1, double* skvmaxY1, double* skvminX2, double* skvminY2, double* skvmaxX2, double* skvmaxY2, double* tetta1, double* tetta2)
{
	ifstream in("wells.txt");
	in >> *n;
	in >> *skvminX1 >> *skvmaxX1 >> *skvminY1 >> *skvmaxY1 >> *skvminX2 >> *skvmaxX2 >> *skvminY2 >> *skvmaxY2 >> *tetta1 >> *tetta2;
	in.close();
}

void createevengridXY(double grid[], int n, double xstart, double xend) // равномерная сетка: n - количество узлов
{
	grid[0] = xstart;
	grid[n - 1] = xend;
	double h = (xend - xstart) / (n - 1);
	for (int i = 1; i < n - 1; i++)
	{
		grid[i] = grid[i - 1] + h;
	}
}
void FindSkvXY(double grid[], int n, double skvmin, double skvmax, double res[]) // узлы равномерной сетки res0 - левая граница res1 - правая граница + 
{
	for (int i = 0; i < n; i++)
	{
		if (grid[i] <= skvmin)
		{
			res[0] = grid[i];
		}
	}

	for (int i = n - 1; i >= 0; i--)
	{

		if (grid[i] >= skvmax)
		{
			res[1] = grid[i];
		}
	}

}
int LeftRightGridXY(double skvgrid[], double skvmin, double skvmax, double left, double right, double k) // разбиение сетки со сгущением слева и справа
{
	double r = (skvmax - skvmin) / 2.0;

	int flag = 1;
	double step = r;
	int n = 1;
	double* skvleftgrid;
	skvleftgrid = new double[1000];

	skvleftgrid[0] = skvmin;
	if ((skvmin - step) > left)
	{
		skvleftgrid[1] = skvmin - step;
		n++;
		for (int i = 2; flag; i++)
		{
			step = step * k;

			if (skvleftgrid[i - 1] - step <= left)
			{
				skvleftgrid[i] = left;
				flag = 0;
				n++;
			}
			else {
				skvleftgrid[i] = skvleftgrid[i - 1] - step;
				n++;
			}
		}
	}
	else
	{
		if (skvleftgrid[0] != left) {
			skvleftgrid[1] = left; n++;
		}
	}

	double* skvrightgrid;
	skvrightgrid = new double[1000];
	flag = 1;
	step = r;
	int m = 1;

	skvrightgrid[0] = skvmax;
	if (skvmax + step <= right)
	{
		skvrightgrid[1] = skvmax + step;
		m++;
		for (int i = 2; flag; i++)
		{
			step = step * k;

			if (skvrightgrid[i - 1] + step >= right)
			{
				if (skvrightgrid[i - 1] == right)
				{
					flag = 0;
				}
				else {
					skvrightgrid[i] = right;
					m++;
					flag = 0;
				}
			}
			else {
				skvrightgrid[i] = skvrightgrid[i - 1] + step; m++;
			}

		}
	}

	int i = 0;

	for (int i = 0; i < n; i++)
	{
		skvgrid[i] = skvleftgrid[n - i - 1];
	}
	i = n;

	for (int j = 0; j < m; j++, i++)
	{
		skvgrid[i] = skvrightgrid[j];
	}
	return i;
}

int cregridXY(double grid[], int n, double resultgrid[], double grid1[], double grid2[], int n1, int n2)// формирование общей сетки
{
	int flag = 1;
	int m = 0;

	for (int i = 0; i < n; i++)
	{

		if (grid[i] == grid1[0] || (resultgrid[m - 1] == grid1[0]))
		{
			for (int i = 0; i < n1; i++)
			{
				resultgrid[m + i] = grid1[i];
			}

			m = m + n1;
			i++;

		}

		if (grid[i] == grid2[0] || (resultgrid[m - 1] == grid2[0]))
		{
			if (resultgrid[m - 1] == grid2[0]) m--;
			for (int i = 0; i < n2; i++)
			{
				resultgrid[m + i] = grid2[i];
			}
			m = m + n2;
			i++;
		}
		while (resultgrid[m - 1] > grid[i]) i++;

		if (resultgrid[m] == grid1[0] || resultgrid[m] == grid2[0] || resultgrid[m - 1] == grid[i]) m--;

		resultgrid[m] = grid[i];
		m++;

	}
	return m;
}
void LocalElem(double gridX[], int nX, double gridY[], int nY)
{
	int m = 1;
	FILE* out;
	fopen_s(&out, "out.txt", "w");
	fprintf_s(out, "%d %d\n", (nX - 1) * (nY - 1), 4);

	for (int i = 0; i < nY - 1; i++)
	{
		for (int j = 0; j < nX - 1; j++)
		{
			fprintf_s(out, "%d %d %d %d \n", m - 1, m, m + nX - 1, m + nX);
			m++;
		}
		m++;
	}

	fclose(out);

}

void LocalElemTriangle(double gridX[], int nX, double gridY[], int nY)
{
	int m = 0;
	FILE* out;
	fopen_s(&out, "outTriangle.txt", "w");
	fprintf_s(out, "%d %d\n", 2 * (nX - 1) * (nY - 1), 3);

	for (int i = 0; i < nY - 1; i++)
	{
		for (int j = 0; j < nX - 1; j++)
		{
			fprintf_s(out, "%d %d %d \n%d %d %d \n", m, m + 1, m + nX, m + 1, m + nX, m + nX + 1);
			m++;
		}
		m++;
	}

	fclose(out);
}

void BC2_Triangle(double gridX[], int nX, double gridY[], int nY, double skvminX1, double skvminY1, double skvminX2, double skvminY2, double tetta1, double tetta2)
{
	FILE* out;
	fopen_s(&out, "BC2_Tri.txt", "w");
	int x = 0;
	int y = 0;
	for (int i = 0; (i < nX && gridX[i] <= skvminX1); i++)
	{
		x = i + 1;
	}
	for (int j = 0; (j < nY && gridY[j] <= skvminY1); j++)
	{
		y = j;
	}

	int x1y1 = x + nX * y;
	int x2y1 = x + 1 + nX * y;
	int x1y2 = x1y1 + nX;
	int x2y2 = x1y2 + 1;
	int elem;
	fprintf_s(out, "%d\n", 8);
	if (y != 0)
	{
		elem = (x2y1 - y) * 2;
	}
	else elem = x2y1 - 1;
	fprintf_s(out, "%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n", elem - 1, 0, 1, tetta1, elem - 1, 2, 0, tetta1, elem, 1, 2, tetta1, elem, 0, 2, tetta1);


	x = 0;
	y = 0;
	for (int i = 0; (i < nX && gridX[i] <= skvminX2); i++)
	{
		x = i + 1;
	}
	for (int j = 0; (j < nY && gridY[j] <= skvminY2); j++)
	{
		y = j;
	}

	x1y1 = x + nX * y;
	x2y1 = x + 1 + nX * y;
	x1y2 = x1y1 + nX;
	x2y2 = x1y2 + 1;


	if (y != 0)
	{
		elem = (x2y1 - y) * 2;
	}
	else elem = x2y1 - 1;
	fprintf_s(out, "%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n", elem - 1, 0, 1, tetta2, elem - 1, 2, 0, tetta2, elem, 1, 2, tetta2, elem, 0, 2, tetta2);

	fclose(out);
}

void BC2(double gridX[], int nX, double gridY[], int nY, double skvminX1, double skvminY1, double skvminX2, double skvminY2, double tetta1, double tetta2)
{


	FILE* out;
	fopen_s(&out, "BC2.txt", "w");
	int x = 0;
	int y = 0;
	for (int i = 0; (i < nX && gridX[i] <= skvminX1); i++)
	{
		x = i + 1;
	}
	for (int j = 0; (j < nY && gridY[j] <= skvminY1); j++)
	{
		y = j;
	}

	int x1y1 = x + nX * y;
	int x2y1 = x + 1 + nX * y;
	int x1y2 = x1y1 + nX;
	int x2y2 = x1y2 + 1;

	fprintf_s(out, "%d\n", 8);
	fprintf_s(out, "%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n", x1y1 - y - 1, 0, 1, tetta1, x1y1 - y - 1, 1, 3, tetta1, x1y1 - y - 1, 3, 2, tetta1, x1y1 - y - 1, 2, 0, tetta1);


	x = 0;
	y = 0;
	for (int i = 0; (i < nX && gridX[i] <= skvminX2); i++)
	{
		x = i + 1;
	}
	for (int j = 0; (j < nY && gridY[j] <= skvminY2); j++)
	{
		y = j;
	}

	x1y1 = x + nX * y;
	x2y1 = x + 1 + nX * y;
	x1y2 = x1y1 + nX;
	x2y2 = x1y2 + 1;

	fprintf_s(out, "%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n%d %d %d %.15e\n\n", x1y1 - y - 1, 0, 1, tetta2, x1y1 - y - 1, 1, 3, tetta2, x1y1 - y - 1, 3, 2, tetta2, x1y1 - y - 1, 2, 0, tetta2);

	fclose(out);

}

void BC1(double gridX[], int nX, double gridY[], int nY)
{

	FILE* out;
	fopen_s(&out, "BC1.txt", "w");
	fprintf_s(out, "%d\n", 2 * nX + 2 * (nY - 2));

	for (int i = 0; i < nX; i++) //нижняя граница
	{
		fprintf_s(out, "%d\n", i);
	}

	for (int i = 1; i < nY - 1; i++) //левая и правая границы
	{
		fprintf_s(out, "%d\n", nX * i);
		fprintf_s(out, "%d\n", nX * (i + 1) - 1);
	}

	for (int i = nX - 1; i >= 0; i--) //верхняя граница
	{
		fprintf_s(out, "%d\n", nX * nY - i - 1);
	}

	fclose(out);
}
void outgrid(double gridX[], int nX, double gridY[], int nY)
{
	FILE* out;
	fopen_s(&out, "grid.txt", "w");
	fprintf_s(out, "%d\n", nX * nY);

	for (int i = 0; i < nY; i++)
	{
		for (int j = 0; j < nX; j++)
		{
			fprintf_s(out, "%.15e %.15e\n", gridX[j], gridY[i]);
		}
	}


	fclose(out);
}
void out_mat(double ki_1, double ki_2, double  phi, double s2, int nX, double x[])
{
	FILE* out;
	fopen_s(&out, "mat.txt", "w");
	int n = nX / 2;
	fprintf_s(out, "%.15e %.15e\n%d\n%.15e %.15e\n%.15e %.15e\n", s2, phi, 2, ki_1, x[n], ki_2, x[nX - 1]);
	fclose(out);
}