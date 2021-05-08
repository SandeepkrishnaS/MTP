#include<bits/stdc++.h>
#include<omp.h>

using namespace std;

#define NX 201
#define NY 201


double RMSfind(int nx, int ny, double a[NX][NY])
{
    double v = 0.0;

    for(int j = 0; j < ny; j++)
        for(int i = 0; i < nx; i++)
            v += a[i][j] * a[i][j];
        
    v = sqrt(v/(double)(nx * ny));
    return v;
}

double ExactSol(double x, double y)
{
    return (double)(sin(2.0*x*y));
}

double ExactDelSq(double x, double y)
{
    return (double)(-4.0*sin(2.0*x*y)*(x*x+y*y));
}

void RightSide(int nx, int ny, double rhsinit[NX][NY])
{
    double x,y;

    for(int j = 0; j < ny; j++)
    {
        y = (double)(j)/(double)(ny-1);
        for(int i = 0; i < nx; i++)
        {
            x = (double)(i)/(double)(nx-1);
            if(i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                rhsinit[i][j] = ExactSol(x, y);
            else
                rhsinit[i][j] = -ExactDelSq(x, y);
        }
    }
    return;
}
 
void OSiter(int nx, int ny, double dx, double dy, double rhsinit[NX][NY], int iterold, int iternew, double mat[NX][NY], double mat_new[NX][NY]) //OneStepIteration()
{
    int i,j,iter;

    #pragma omp parallel shared ( dx, dy, rhsinit, iternew, iterold, nx, ny, mat, mat_new ) private ( i, iter, j )
    for(iter = iterold + 1; iter <= iternew; iter++)
    {
        #pragma omp for          //loop construct for splitting the pathway of for loop to the threads
        for ( j = 0; j < ny; j++ )
            for ( i = 0; i < nx; i++ )
                 mat[i][j] = mat_new[i][j];
        
        //cout<<"-->"<<omp_get_thread_num()<<"<--";
        #pragma omp for
        for(j = 0; j < ny; j++)
            for(i = 0; i < nx; i++)
                if(i == nx - 1 || j == ny - 1 || i == 0 || j == 0)
                    mat_new[i][j] = rhsinit[i][j];
                else
                    mat_new[i][j] = 0.25 * (mat[i-1][j] + mat[i][j+1] + mat[i][j-1] + mat[i+1][j] - rhsinit[i][j] * dx * dy);
    }
    return;
}

int main(int argc, char *argv[])
{
    bool converged = false;
    long nx = NX, ny = NY;
    double tolerance = 0.000001;
    double dx = 1.0 / ( double ) ( nx - 1 ), dy = 1.0 / ( double ) ( ny - 1 );

    cout<<"Number of processors : "<<omp_get_num_procs()<<endl;
    
    int id; 
    # pragma omp parallel
    {
        id = omp_get_thread_num();
        if(id == 0)
            cout<<"Maximum threads : "<<omp_get_num_threads()<<endl; 
    }
    
    cout<<"Solving del^2 P = 4*(x^2 + y^2)*sin(2*x*y)"<<endl;
    cout<<"Region 0 <= X <= 1, 0 <= Y <= 1"<<endl;
    cout<<"Grid : 200x200 "<<endl;
    
    double RHSinit[NX][NY];
    RightSide(nx, ny, RHSinit);

    int i, j;
    double mat[NX][NY], mat_err, mat_new[NX][NY], mat_new_err ;
    double mat_exact[NX][NY], mat_diff[NX][NY];
    
    for(j = 0; j < ny; j++)
    {
        for(i = 0; i < nx; i++)
        {
        if(i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
            mat_new[i][j] = RHSinit[i][j];
        else
            mat_new[i][j] = 0.0;
        }
    }
    
    mat_new_err = RMSfind(nx, ny, mat_new);
    
    double x, y;
    for(j = 0; j < ny; j++)
    {
        y = (double)(j)/(double)(ny-1);
        for(i = 0; i < nx; i++)
        {
            x = (double)(i)/(double)(nx-1);
            mat_exact[i][j] = ExactSol(x, y);
        }
    }
    
    mat_err = RMSfind(nx, ny, mat_exact);
    mat_err++; mat_err--;

    for(j = 0; j < ny; j++)
    {
        for(i = 0; i < nx; i++)
            mat_diff[i][j] = mat_new[i][j] - mat_exact[i][j];
    }

    auto start = chrono::high_resolution_clock::now();   //time start
    ios_base::sync_with_stdio(false);

    long itnew = 0, itold, itmax = 1000000;
    double diff;

    while(true)
    {
        itold = itnew;
        itnew = itold + 100; //Change 100 to 1 if required to check no. of steps 

        OSiter(nx, ny, dx, dy, RHSinit, itold, itnew, mat, mat_new);

        mat_err = mat_new_err;
        mat_new_err = RMSfind( nx, ny, mat_new);

        for ( j = 0; j < ny; j++ )
            for ( i = 0; i < nx; i++ )
                mat_diff[i][j] = mat_new[i][j] - mat[i][j];
        
        diff = RMSfind( nx, ny, mat_diff );

        for ( j = 0; j < ny; j++ )
            for ( i = 0; i < nx; i++ )
                mat_diff[i][j] = mat_new[i][j] - mat_exact[i][j];

        if(itold > itmax)
            break;

        if(diff <= tolerance)
        {
            converged = true;
            break;
        }

    }

    if(converged)
        cout<<"The iteration has converged "<<endl<<"Steps : "<<itold<<" - "<<itnew<<endl;
    else
        cout<<"The iteration has NOT converged.\n";

    auto end = chrono::high_resolution_clock::now();
    double time_taken = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    cout<<"Time taken : "<<fixed<<time_taken<<setprecision(9)<< " sec"<<endl;
    return 0;
}

# undef NX
# undef NY
