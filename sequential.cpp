#include<bits/stdc++.h>

using namespace std;

#define NX 201
#define NY 201

double ExactSol(double x, double y)
{
    return (double)(sin(2.0*x*y));
}

double ExactDelSq(double x, double y)
{
    return (double)(-4.0*sin(2.0*x*y)*(x*x+y*y));
}

double RMSfind(int nx, int ny, double a[NX][NY])
{
    double v = 0.0;

    for(int j = 0; j < ny; j++)
        for(int i = 0; i < nx; i++)
            v += a[i][j] * a[i][j];

    v = sqrt ( v / ( double ) ( nx * ny )  );

    return v;
}


void RightSide(int nx, int ny, double a[NX][NY])
{
    double x,y;

    for(int j = 0; j < ny; j++)
    {
        y = (double)(j)/(double)(ny-1);
        for(int i = 0; i < nx; i++)
        {
            x = (double)(i)/(double)(nx-1);
            if(i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                a[i][j] = ExactSol(x, y);
            else
                a[i][j] = -ExactDelSq(x, y);
        }
    }

    return;
}

void OSiter(int nx, int ny, double dx, double dy, double rhsinit[NX][NY], double mat[NX][NY], double mat_new[NX][NY]) 
{
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            if(i == 0 || j == 0 || i == nx - 1 || j == ny - 1)
                mat_new[i][j] = mat[i][j];
            else
                mat_new[i][j] = 0.25 * (mat[i-1][j] + mat[i][j+1] + mat[i][j-1] + mat[i+1][j] - rhsinit[i][j] * dx * dy );
        }
    }
    return;
}




int main(int argc, char *argv[])
{
    cout << "Solving del^2 P = 4*(x^2 + y^2)*sin(2*x*y)"<<endl;
    cout << "Region 0 <= X <= 1, 0 <= Y <= 1"<<endl;
    cout << "Grid : 200x200 "<<endl;
    
    bool converged = false;
    long nx = NX, ny = NY;
    double dx = 1.0 / ( double ) ( nx - 1 ), dy = 1.0 / ( double ) ( ny - 1 ), x, y;
    double toleranceLevel = 0.000001;

    double RHSinit[NX][NY];
    RightSide(nx, ny, RHSinit); //RHSinit is passed by reference. Arrays are always passed by ref. We are initializing RHSinit here with actual values

    int i, j;
    double mat[NX][NY], mat_err, mat_new[NX][NY], mat_new_err;
    double mat_diff[NX][NY], mat_exact[NX][NY];
    
    for(j = 0; j < ny; j++)
    {
        y = (double)(j)/(double)(ny-1);
        for(i = 0; i < nx; i++)
        {
            x = (double)(i)/(double)(nx-1);
            if(i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                mat_new[i][j] = RHSinit[i][j] ; //new 2D array we are going to work with, where inside points are same as RHSinit and boundary points are set to 0.
            else
                mat_new[i][j] = 0.0; //setting bpoundary points to 0
        }
    }

    mat_new_err = RMSfind(nx, ny, mat_new);

    for(j = 0; j < ny; j++)
    {
        y = (double)(j)/(double)(ny - 1 );
        for(i = 0; i < nx; i++)
        {
            x = (double)(i)/(double)(nx - 1);
            mat_exact[i][j] = ExactSol(x, y);
        }
    }
    

    mat_err = RMSfind(nx, ny, mat_exact);

    mat_err ++; mat_err--; //getting rid of warning

    auto start = chrono::high_resolution_clock::now();
    ios_base::sync_with_stdio(false);

    long iter;
    long iter_max = 1000000; //in case it doesn't converge 
    
    double diff;
    
    for(iter = 1; iter <= iter_max; iter++)
    {
        for ( j = 0; j < ny; j++ )
            for ( i = 0; i < nx; i++ )
                mat[i][j] = mat_new[i][j];

        OSiter(nx, ny, dx, dy, RHSinit, mat, mat_new);

        mat_err = mat_new_err;
        mat_new_err = RMSfind(nx, ny, mat_new);

        for(j = 0; j < ny; j++)
            for(i = 0; i < nx; i++)
                mat_diff[i][j] = mat_new[i][j] - mat[i][j];

        diff = RMSfind(nx, ny, mat_diff);

        if(diff <= toleranceLevel)
        {
            converged = true;
            break;
        }

    }

    if(converged)
        cout<<"The iteration has converged"<<endl<<"Num of steps : "<<iter<<endl;
    else
        cout<<"The iteration has NOT converged."<<endl;
    
    auto end = chrono::high_resolution_clock::now();

    // total time taken
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;

    cout<<"Time taken : "<<fixed<<time_taken<<setprecision(9)<<" sec"<<endl;

    return 0;
}

# undef NX
# undef NY
