void input_mesh(int* n_x, int* n_y, double* k);

//void input_xy(double* xmin, double* xmax, double* ymin, double* ymax, double* ki, double* phi, double* s2);
void input_xy(double* xmin, double* xmax, double* ymin, double* ymax, double* ki_1, double* ki_2, double* phi, double* s2);

void input_wells(int* n, double* skvminX1, double* skvminY1, double* skvmaxX1, double* skvmaxY1, double* skvminX2, double* skvminY2, double* skvmaxX2, double* skvmaxY2, double* tetta1, double* tetta2);

void createevengridXY(double grid[], int n, double xstart, double xend);


void FindSkvXY(double grid[], int n, double skvmin, double skvmax, double res[]);


int LeftRightGridXY(double skvgrid[], double skvmin, double skvmax, double left, double right, double k);

int cregridXY(double grid[], int n, double resultgrid[], double grid1[], double grid2[], int n1, int n2);

//void LocalElem(double gridX[], int nX, double gridY[], int nY);

void LocalElemTriangle(double gridX[], int nX, double gridY[], int nY);

void BC2_Triangle(double gridX[], int nX, double gridY[], int nY, double skvminX1, double skvminY1, double skvminX2, double skvminY2, double tetta1, double tetta2);
//void BC2(double gridX[], int nX, double gridY[], int nY, double skvminX1, double skvminY1, double skvminX2, double skvminY2, double tetta1, double tetta2);

void BC1(double gridX[], int nX, double gridY[], int nY);

void outgrid(double gridX[], int nX, double gridY[], int nY);

//void out_mat(double ki, double  phi, double s2);
void out_mat(double ki_1, double ki_2, double  phi, double s2, int nX, double x[]);