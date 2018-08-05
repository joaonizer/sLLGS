//
//
//  THIS PROGRAM COMPUTES THE DEMAGNETIZING TENSOR OF A SINGLE-DOMAIN SLANTED PRISM
//  Multiple integrals are evaluated through the Monte Carlo method
//
//
//  O arquivo de entrada (IN.dat) deve conter:
//	px(1),py(1): coordenadas do canto superior esquerdo (nm)
//	px(2),py(2): coordenadas do canto superior direito  (nm)
//	px(3),py(3): coordenadas do canto inferior direito  (nm)
//	px(4),py(4): coordenadas do canto inferior esquerdo (nm)
//	t	   : espessura (nm)
//
//	*************** RESTRICAO ****************
//	-------> px(1) e px(4) DEVEM ser iguais
//	-------> px(2) e px(3) DEVEM ser iguais
//	******************************************
//
//	O arquivo de saida (OUT_demag3D.dat) contem:
//	Linha 1: d1, d2, d3
//	Linha 2: d4, d5, d6
//	Linha 3: d7, d8, d9
//
//	Os elementos do tensor DESMAGNETIZACAO sao obtidos fazendo
//
//	D11=d1/(4*Pi*V), D12=d2/(4*Pi*V), D13=d3/(4*Pi*V)
//	D21=d4/(4*Pi*V), D22=d5/(4*Pi*V), D23=d6/(4*Pi*V)
//	D31=d7/(4*Pi*V), D32=d8/(4*Pi*V), D33=d9/(4*Pi*V)
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
//	Ud = K/2*(
//	          (d1*COS(phi)^2 + d5*SIN(phi)^2 + (d2+d4)*SIN(phi)*COS(phi))*SIN(theta)^2 +
//		      ((d3+d7)*COS(phi) + (d6+d8)*SIN(phi))*SIN(theta)*COS(theta) +
//		       d9*COS(theta)^2
//		     )
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
#include <iomanip>  //Required for setprecision(), setw()
#include <cstdlib>  //Required for srand(), rand()
//#include <ctime>    //Required for time
#include <cmath>    //Required for pow()
#include <fstream>  //Required for ifstream, ofstream

using namespace std;

// Global variables
const int NMC(1e6);
double w, t;             // w: particle's width, t: particle's thickness
double abu[2], abd[2];   //global vectors: yu(x)= abu[0]x +abu[1] and yd(x)= abd[0]x +abd[1]
// Function prototypes
double  frand(double, double);
double *frand_xz();
double *fu(double []);
double  fy(double, const double []);
double *demag(const double [], const double []);


int main()
{

    double px[4], py[4], *ptr, vol;
    double alfa2, beta2, alfa5, beta5, d[18];

    // Inicia seed de acordo com o clock: gera # pseudoaleatorios cada vez que o programa rodar
//    srand((int)time(NULL));
    srand(1);

    //Define file streams for input and output
//    ifstream fin("IN.dat");
//    ofstream fout("OUT.dat");

//    //Check for possible errors.
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

//  P1={x1,y1},...,P4={x4,y4}; xy={P1,P2,P3,P4,{-t/2,t/2}}={{xy[0][0],xy[0][1]},{xy[1][0],xy[1][1]},...,{xy[4][0],xy[4][1]}}
    //for(int i=0; i<4; ++i) fin >> px[i] >> py[i];
	for(int i=0; i<4; ++i) std::cin >> px[i] >> py[i];
    std::cin >> t;
    //fin.close();

//  definition of functions yu(x) and yd(x) coefficients: yu =au*x+bu=ab[0]x+ab[1] and yd =ad*x+bd=ab[2]x+ab[3]
//  the vectors abu and abd are global
         w =  px[1] - px[0];
    abu[0] = (py[1] - py[0])/w;                       // au
	abu[1] =  py[0] - abu[0]*px[0];                   // bu
	abd[0] = (py[2] - py[3])/w;                       // ad
	abd[1] =  py[3] - abd[0]*px[3];                   // bd

	alfa2 = -abu[0];
	beta2 =  1;
	alfa5 =  abd[0];
	beta5 = -1;

    ptr = demag(px,py);

    d[0] = ptr[0]+alfa2*ptr[1]-ptr[3]+alfa5*ptr[4];
    d[1] = beta2*ptr[1]+beta5*ptr[4];
//    d[2] = ptr[2]-ptr[5];
    d[2] = 0;

    d[3] = ptr[6]+alfa2*ptr[7]-ptr[9]+alfa5*ptr[10];
    d[4] = beta2*ptr[7]+beta5*ptr[10];
//    d[5] = ptr[8]-ptr[11];
    d[5] = 0;

//    d[6] = ptr[12]+alfa2*ptr[13]-ptr[15]+alfa5*ptr[16];
    d[6] = 0;
//    d[7] = beta2*ptr[13]+beta5*ptr[16];
    d[7] = 0;
    d[8] = ptr[14]-ptr[17];

//    fout << setprecision(3) << std::scientific;
    std::cout << endl;//setprecision(3) << std::scientific;

    for(int i=0; i<3; ++i)
    {
        //fout << setw(14) << *(d+3*i) <<setw(14)<< *(d+3*i+1) <<setw(14)<< *(d+3*i+2) << endl;
        std::cout << setw(14) << *(d+3*i) <<setw(14)<< *(d+3*i+1) <<setw(14)<< *(d+3*i+2) << endl;
    }
    std::cout << endl;
    vol = 0.5*(py[0]-py[3]+py[1]-py[2])*w*t;
    //cout << setw(14) << "Volume = " << vol << endl;
    //fout << setw(14) << vol << endl;



    //fout.close();
    return 0;
}



//******************************************************************************************/
double frand(double xmin, double xmax)
{
    return ((double)rand()/RAND_MAX)*(xmax-xmin)+xmin;
}


//******************************************************************************************/
double *frand_xz()
{
    // input vector is in={xmin,xmax,ymin,ymax,zmin,zmax,xpmin,xpmax,ypmin,ypmax,zmin,zmax}
    static double out[4];
    double in[8]={-0.5*w,0.5*w,-0.5*t,0.5*t,-0.5*w,0.5*w,-0.5*t,0.5*t};
    for(int i=0; i<4; ++i) out[i]=((double)rand()/RAND_MAX)*(in[2*i+1]-in[2*i])+in[2*i];
    return out;
}


//******************************************************************************************/
double *fu(double in[6])
{
     // input vector is  in={x(0),y(1),z(2),xp(3),yp(4),zp(5)}
    static double out[3];
    double den;
    den = pow((in[0]-in[3])*(in[0]-in[3]) + (in[1]-in[4])*(in[1]-in[4]) + (in[2]-in[5])*(in[2]-in[5]), -1.5);
    for(int i=0; i<3; ++i) *(out+i) = (*(in+i+3) - *(in+i))*den;
    return out;
}


//******************************************************************************************/
double fy(double x, const double ab[2])
{
    return ab[0]*x + ab[1];
}



//******************************************************************************************/
double *demag(const double px[3], const double py[3])
{
    static double abc_out[18];
    double *abc, a[6]={0}, b[6]={0}, c[6]={0};
    double *xz, x[6], yi[4];
    double yd_x,  yu_x,  Dy;
    double yd_xp, yu_xp, Dyp;

    double s[6];
    s[0] = w*t*t*(py[1]-py[2]);
    s[1] = w*w*t*t;
    s[2] = w*w*t;
    s[3] = w*t*t*(py[0]-py[3]);
    s[4] = s[1];
    s[5] = s[2];

    // main loop of the program
    for(int i=0; i<NMC; ++i)
    {

        xz    = frand_xz();           //obten #aleatorios x(0),z(1),x'(2),z'(3)
        yd_x  = fy(xz[0],abd);           // yd(x)
        yu_x  = fy(xz[0],abu);           // yu(x)
        Dy    = yu_x - yd_x;
        yd_xp = fy(xz[2],abd);           // yd(x')
        yu_xp = fy(xz[2],abu);
        Dyp   = yu_xp - yd_xp;           // yu(x')
        yi[0] = frand(yd_x, yu_x);        // yd(x)  <  y(x)  < yu(x)
        yi[1] = frand(yd_xp,yu_xp);      // yd(x') <  y(x') < yu(x')
        yi[2] = frand(py[2],py[1]);      // P3y    <  y     < P2y
        yi[3] = frand(py[3],py[0]);      // P4y    <  y     < P1y


        x[0] = xz[0];
        x[1] = yi[0];
        x[2] = xz[1];

        //********* S1 ************OK
        x[3] = 0.5*w;  // x'
        x[4] = yi[2];  // y'
        x[5] = xz[3];  // z'
        abc  = fu(x);
        *a  += *abc     * Dy;
        *b  += *(abc+1) * Dy;
        *c  += *(abc+2) * Dy;

        //********* S2 ************OK
        x[3] = xz[2];  // x'
        x[4] = yu_xp;  // y'
        x[5] = xz[3];  // z'
        abc  = fu(x);
        *(a+1) += *abc     * Dy;
        *(b+1) += *(abc+1) * Dy;
        *(c+1) += *(abc+2) * Dy;

        //********* S3 ************OK
        x[3] = xz[2];
        x[4] = yi[1];
        x[5] = 0.5*t;
        abc  = fu(x);
        *(a+2) += *abc     * Dy*Dyp;
        *(b+2) += *(abc+1) * Dy*Dyp;
        *(c+2) += *(abc+2) * Dy*Dyp;

        //********* S4 ************OK
        x[3] = -0.5*w;
        x[4] = yi[3];
        x[5] = xz[3];
        abc  = fu(x);
        *(a+3) += *abc     * Dy;
        *(b+3) += *(abc+1) * Dy;
        *(c+3) += *(abc+2) * Dy;

        //********* S5 ************OK
        x[3] = xz[2];  // x'
        x[4] = yd_xp;  // y'
        x[5] = xz[3];  // z'
        abc  = fu(x);
        *(a+4) += *abc     * Dy;
        *(b+4) += *(abc+1) * Dy;
        *(c+4) += *(abc+2) * Dy;

        //********* S6 ********************************************
        x[3] = xz[2];
        x[4] = yi[1];
        x[5] = -0.5*t;
        abc  = fu(x);
        *(a+5) += *abc     * Dy*Dyp;
        *(b+5) += *(abc+1) * Dy*Dyp;
        *(c+5) += *(abc+2) * Dy*Dyp;

    }


    for(int i=0; i<6; ++i)
    {
        abc_out[i]    = s[i]*a[i]/double(NMC);
        abc_out[i+6]  = s[i]*b[i]/double(NMC);
        abc_out[i+12] = s[i]*c[i]/double(NMC);
    }

    return abc_out;

}








