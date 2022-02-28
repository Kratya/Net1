#pragma once

int LeftGridAAAAAAAA(double skvgrid[], double skvmin, double skvmax, double left, double right, double g[], int n)
{
	auto fake = (right - left) / 2;
	auto leftX = skvmin - g[n] + fake;
	auto h = g[1] - g[0];


	auto k = 1 - h / leftX;
	auto x0 = g[0];
	skvgrid[0] = x0;
	skvgrid[1] = x0 = h;

	int i = 1;
	while (x0 + h * k < leftX - fake)
	{
		h = h * k;
		x0 += h;
		skvgrid[i++] = x0;
	}
	skvgrid[i++] = leftX - fake;
	return i;
}

int RightGridAAAAAAAA(double skvgrid[], double skvmin, double skvmax, double left, double right, double g[], int n)
{
	auto fake = (right - left) / 2;
	auto leftX = g[n - 1] - skvmax + fake;
	auto h = g[1] - g[0];

	double* temp = new double[10000];


	auto k = 1 - h / leftX;
	auto x0 = g[0];
	temp[0] = x0;
	temp[1] = x0 = h;

	int i = 1;
	while (x0 + h * k < leftX - fake)
	{
		h = h * k;
		x0 += h;
		temp[i++] = x0;
	}
	temp[i++] = leftX - fake;
	i--;
	int j = 0;
	h = g[1] - g[0];
	for (i; i >= 0; i--)
	{
		skvgrid[j++] = leftX - temp[i] + skvmax - fake;
	}
	return j;
}

int MidGridAAAAAAAA(double skvgrid[], double skvmin, double skvmax, double left, double right, double g[], int n)
{
	double* temp1 = new double[10000];
	double* temp2 = new double[10000];
	auto nx1 = LeftGridAAAAAAAA(temp1, skvmin, skvmax, left, right, g, n);
	auto nx2 = RightGridAAAAAAAA(temp2, skvmin, skvmax, left, right, g, n + 1);
	for (int i = 0; i < nx1; i++)
	{
		temp1[i] += temp2[nx2 - 1];
	}
	int i = 0;
	for (i; i < nx2; i++)
	{
		skvgrid[i] = temp2[i];
	}
	int j = 1;
	for (j; j < nx2; j++)
	{
		skvgrid[i++] = temp1[j];
	}
	return i;
}

int createGridX(double res[], double left[], double mid[], double right[], int nl, int nr, int nm)
{
	int j = 0;
	for (int i = 0; i < nl; i++)
	{
		res[j++] = left[i];
	}
	for (int i = 0; i < nm; i++)
	{
		res[j++] = mid[i];
	}
	for (int i = 0; i < nr; i++)
	{
		res[j++] = right[i];
	}
	return nl + nm + nr;
}

int createGridY(double res[], double left[], double right[], int nl, int nr)
{
	int j = 0;
	for (int i = 0; i < nl; i++)
	{
		res[j++] = left[i];
	}
	for (int i = 0; i < nr; i++)
	{
		res[j++] = right[i];
	}
	return nl + nr;
}