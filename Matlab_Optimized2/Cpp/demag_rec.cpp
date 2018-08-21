//
//
//  THIS PROGRAM COMPUTES THE DEMAGNETIZING TENSOR OF A SINGLE-DOMAIN SQUARE PRISM
//  Multiple integrals are evaluated through the Monte Carlo method
//
//
//  O arquivo de entrada (IN.dat) deve conter tres linhas:
//
//  l (comprimento, em nm)
//  w (largura, nm)
//  t (espessura, nm)
//
//
//  O arquivo de saida (OUT.dat) contem uma linha:
//
//	d1, d2, d3
//
//	Os elementos do tensor DESMAGNETIZACAO sao obtidos fazendo
//
//	D11=d1/(4*Pi*V), D12=0,           D13=0
//	D21=0,           D22=d2/(4*Pi*V), D23=0
//	D31=0,           D32=0,           D33=d3/(4*Pi*V)
//
//	onde V e dado em nm^3
//
//
//  O (vetor) campo de desmagnetizacao normalizado é dado por
//
//  Hd/Ms = -(D11*mx+D12*my+D13*mz, D21*mx+D22*my+D23*mz, D31*mx+D32*my+D33*mz)
//        = -(D11*mx, D22*my, D33*mz)
//
//  Onde
//  mx = SIN(theta)*COS(phi)
//  my = SIN(theta)*SIN(phi)
//  mz = COS(theta)
//
//
//	A energia de desmagnetizacao (eV) e calculada como
//
//	Ud = K/2*((d1*COS(phi)^2 + d2*SIN(phi)^2)*SIN(theta)^2 + d3*COS(theta)^2)
//     = K/2*(d1*mx^2 + d2*my^2 + d3*mz^2)
//	onde
//
//	K     = mu0*Ms^2*JtoeV*10e-27/(4*Pi)
//	Ms    = 800 kA/m (valor do Permalloy)
//	mu0   = 4*Pi*10e-7 H/m
//	JtoeV = 6.242*10e18
//
//
//
#include <iostream>
#include <iomanip>
#include <cstdlib>  //Required for srand(), rand(),
#include <ctime>    //Required for time
#include <cmath>
#include <fstream>  //Required for ifstream, ofstream

using namespace std;

// Function prototypes
double frand(double, double);
double fu1(const double []);
double fu2(const double []);
double fu3(const double []);
double *demag(const double []);
const double NMC(1e6);



int main()
{
    // Declare and initialize global variables
    double wlt[3], *ptr_demag;
    int i, cont;

    // Inicia seed de acordo com o clock: gera # pseudoaleatorios cada vez que o programa rodar
    srand((int)time(NULL));

    //Define file streams for input and output
//    ifstream fin("IN.dat");
//    ofstream fout("OUT.dat");
    //Check for possible errors.
//    if(fin. fail())
//    {
//        cerr << "could not open input file IN.dat\n";
//        exit(1);
//    }
//    if(fout. fail())
//    {
//        cerr << "could not open output file OUT.dat\n";
//        exit(1);
//    }

    // Read particle's length, width and thickness and DIVIDE THEM BY 2
    //for(i=0; i<3; ++i) fin >> wlt[i];
    for(i=0; i<3; ++i) std::cin >> wlt[i];
    for(i=0; i<3; ++i) wlt[i] *= 0.5;
//    fin.close();

    ptr_demag = demag(wlt);

//    fout << setprecision(3) << std::scientific;
    std::cout << endl;;

    for(i=0; i<3; ++i)
    {
//        fout << ptr_demag[i] << setw(14);
        std::cout << ptr_demag[i] << setw(14);
    }
    std::cout << endl;

//    fout.close();
    return 0;
}

// o vetor de entrada de 5 elementos a[]={a0,a1,a2,a3,a4}
// frand gera um vetor com 5 numeros aleatorios, cada um variando de -a(i) < out(i) < a(i)
double *frand(const double a[])
{
    static double out[5];
    for(int i=0; i<5; ++i) out[i]=((double)rand()/RAND_MAX)*(a[i]*2) - a[i];
    return out;
}


inline double fu1(const double wlt[])
{
    double *r, den1, den2;
    // os 5 numeros aleatorios sao x,y,z,yp,zp; Fixo: xp = 0.5*w = wlt[0]
    double in[5]={wlt[0],wlt[1],wlt[2],wlt[1],wlt[2]};

    r = frand(in);

    den1 = pow((r[0]-wlt[0])*(r[0]-wlt[0]) + (r[1]-r[3])*(r[1]-r[3]) + (r[2]-r[4])*(r[2]-r[4]), -1.5);
    den2 = pow((r[0]+wlt[0])*(r[0]+wlt[0]) + (r[1]-r[3])*(r[1]-r[3]) + (r[2]-r[4])*(r[2]-r[4]), -1.5);

    return (wlt[0]-r[0])*den1 + (wlt[0]+r[0])*den2;
}


inline double fu2(const double wlt[])
{
    double *r, den1, den2;
    // os 5 numeros aleatorios sao x,y,z,xp,zp; Fixo: yp = 0.5*l = wlt[1]
    double in[5]={wlt[0],wlt[1],wlt[2],wlt[0],wlt[2]};

    r = frand(in);

    den1 = pow((r[0]-r[3])*(r[0]-r[3]) + (r[1]-wlt[1])*(r[1]-wlt[1]) + (r[2]-r[4])*(r[2]-r[4]), -1.5);
    den2 = pow((r[0]-r[3])*(r[0]-r[3]) + (r[1]+wlt[1])*(r[1]+wlt[1]) + (r[2]-r[4])*(r[2]-r[4]), -1.5);

    return (wlt[1]-r[1])*den1 + (wlt[1]+r[1])*den2;
}


inline double fu3(const double wlt[])
{
    double *r, den1, den2;
    // os 5 numeros aleatorios sao x,y,z,xp,yp; Fixo: zp = 0.5*t = wlt[2]
    double in[5]={wlt[0],wlt[1],wlt[2],wlt[0],wlt[1]};

    r = frand(in);

    den1 = pow((r[0]-r[3])*(r[0]-r[3]) + (r[1]-r[4])*(r[1]-r[4]) + (r[2]-wlt[2])*(r[2]-wlt[2]), -1.5);
    den2 = pow((r[0]-r[3])*(r[0]-r[3]) + (r[1]-r[4])*(r[1]-r[4]) + (r[2]+wlt[2])*(r[2]+wlt[2]), -1.5);

    return (wlt[2]-r[2])*den1 + (wlt[2]+r[2])*den2;
}


double *demag(const double wlt[])
{
    static double d[3]={0};
    double w = wlt[0], l = wlt[1], t = wlt[2];
    double s[] = {32*w*l*l*t*t/double(NMC), 32*w*w*l*t*t/double(NMC), 32*w*w*l*l*t/double(NMC)};
    int i;

    // este � o loop critico do codigo
    for(i=0; i<NMC; ++i)
    {
        d[0] += fu1(wlt);
        d[1] += fu2(wlt);
        d[2] += fu3(wlt);
    }

    for(i=0; i<3; ++i) d[i] *= (s[i]);

    return d;

}








