#include <iostream>
#include <armadillo>
#include <time.h>
#include <math.h>
#include <sstream>

using namespace std;
using namespace arma;


class Cholesky
{
    public:
        mat results;
        double error;
        void get_cholesky(mat A, mat b)
        {
            mat L, inv_l, inv_l_trans, Y;
            //Se calcula el lower de A
            chol(L, A, "lower");
            //La inversa del lower
            inv_l = inv(L);
            //La inversa de la transpuesta
            inv_l_trans = inv(L.t());
            //Se calcula Y con la inversa de L y b
            Y = inv_l*b;
            //Se obtiene el vector resultando al multiplicar L' con Y.
            mat result = inv_l_trans*Y;
            double error = norm(eye(size(A)) - inv(L * L.t())*A);
            this->results = result;
            this->error = error;

        }
};

class Qr
{
    public:
        mat results;
        double error;
        void get_QR(mat A, mat b)
        {
            mat Q, R, Y;
            //Se calcula Q y R
            qr(Q,R,A);

            //Se obtiene Y con la inversa de Q * b
            Y = inv(Q)*b;

            //Se obtiene el vector resultante con la inversa de R * Y
            vec result = inv(R)*Y;

            double error = norm(eye(size(A)) - inv(Q*R)*A);
            this->error = error;
            this->results = result;
        }
};


class Seidel
{
    public:
        mat results;
        mat errors;
        void get_seidel(mat A, mat b)
        {
            int iterations = 0;
            int max_iterations = 1000;
            int n = size(b)[0];
            mat sol = zeros(1,n).t();
            mat errors;
            mat error = zeros(1,1);
            while (iterations < max_iterations)
            {
                mat sol_old = sol;
                for (int i = 0; i < n; i++)
                {
                    double S = 0;
                    for (int j = 0; j < i-1 ; j++)
                    {
                        S = S + A(i,j)*sol(j);  
                    }

                    for (int j = i+1; j < n; j++)
                    {
                        S = S + A(i,j)*sol_old(j);
                    }

                    sol.row(i) = (1/A(i,i))*(b(i)-S);

                }
            mat diference = sol_old - sol;
            error(0,0) = norm(diference);
            errors = join_horiz(errors,error);
            iterations = iterations + 1;
            }
            this->results = sol;
            this->errors = errors;
        }
};

class Givens
{
    public:
        double error;
        mat results;
        void get_givens(mat A, mat b)
        {
            int m = size(A)[1];
            int n = size(A)[0];
            sp_mat Q,R;
            Q = Q.eye(m,n);
            R = A;
            for (int i = 0; i < m; i++)
            {
                for (int k = i+1; k < n ; k++)
                {
                    if (R(k,i) != 0)
                    {
                        double root = sqrt(R(k,i)*R(k,i) + R(i,i)*R(i,i));
                        double s = -R(k,i)/root;
                        double c = R(i,i)/root;
                        sp_mat g;
                        g = g.eye(m,n);
                        g(k,k) = c;
                        g(i,i) = c; 

                        g(i,k) = s;
                        g(k,i) = -s;

                        Q = Q*g;
                        R = g.t()*R;
                    }   
                }     
            }
            mat mat_Q(Q);
            mat mat_R(R);
            mat Y = inv(mat_Q)*b;
            mat sol = inv(mat_R)*Y;
            double error = norm(eye(size(A)) - inv(mat_Q*mat_R)*A);
            this->results = sol;
            this->error = error;
        }
};

class LeastSquares
{
    public:
        mat results;
        double error;
        void get_leastSquares (mat A, mat b)
        {
            mat A_trans = A.t();
            mat result = inv(A_trans*A)*A_trans*b;
            double error = 	norm(eye(size(A))-inv(inv(A_trans*A))*A_trans);
            this->results = result;
            this->error = error;
        }
};





timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
	    temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

int main()
    { 
   	timespec time1_chol, time2_chol;
   	timespec time1_qr, time2_qr;
   	timespec time1_givens, time2_givens;
   	timespec time1_seidel, time2_seidel;
   	timespec time1_leastSquares, time2_leastSquares;

    mat A;
    mat b;

    //A.load("./Sistemas/289x289/A289.dat", arma::raw_ascii);
    //b.load("./Sistemas/289x289/B289.dat", arma::raw_ascii);

    /*  A << 4 << 10 << 8 << endr 
    << 10 << 26 << 26 << endr
    << 8 << 26 << 91 << endr;

    b << 44 << endr
    << 128 << endr
    << 214 << endr;
*/
    A.load("./Sistemas/1089x1089/A1089.dat", arma::raw_ascii);
    b.load("./Sistemas/1089x1089/B1089.dat", arma::raw_ascii);

    //A.load("./Sistemas/4225x4225/A4225.dat", arma::raw_ascii);
    //b.load("./Sistemas/4225x4225/B4225.dat", arma::raw_ascii);
    

   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1_chol);
    Cholesky cholesky;
    cholesky.get_cholesky(A,b);
    double chol_error = cholesky.error;
    mat chol_result = cholesky.results;
   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2_chol);
	//cout<<diff(time1_chol,time2_chol).tv_sec<<":"<<diff(time1_chol,time2_chol).tv_nsec<<endl;

   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1_qr);
    Qr qr_method;
    qr_method.get_QR(A,b);
    double QR_error = qr_method.error;
    mat QR_result = qr_method.results;
   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2_qr);
	//cout<<diff(time1_qr,time2_qr).tv_sec<<":"<<diff(time1_qr,time2_qr).tv_nsec<<endl;

   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1_seidel);
    Seidel seidel;
    seidel.get_seidel(A,b);
    mat seidel_result = seidel.results;
    mat seidel_errors = seidel.errors;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2_seidel);
	//cout<<diff(time1_seidel,time2_seidel).tv_sec<<":"<<diff(time1_seidel,time2_seidel).tv_nsec<<endl;

   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1_givens);
    Givens givens;
    //givens.get_givens(A,b);
    //mat givens_result = givens.results;
    //double givens_error = givens.error;
   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2_givens);
	//cout<<diff(time1_givens,time2_givens).tv_sec<<":"<<diff(time1_givens,time2_givens).tv_nsec<<endl;

   	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1_leastSquares);
    LeastSquares leastSquares;
    leastSquares.get_leastSquares(A,b);
    mat leastSquares_result =  leastSquares.results;
    double leastSquares_error = leastSquares.error;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2_leastSquares);
	//cout<<diff(time1_leastSquares,time2_leastSquares).tv_sec<<":"<<diff(time1_leastSquares,time2_leastSquares).tv_nsec<<endl;

    ofstream file_1089 ("1089x1089results.txt");

    file_1089 << "---------------- Errores ------------------ \n\n";

    file_1089 << "Seidel: \n";
    file_1089 << seidel_errors;
    file_1089 << "\n";

    file_1089 << "QR: \n";
    file_1089 << QR_error;
    file_1089 << "\n";

    file_1089 << "Cholesky: \n";
    file_1089 << chol_error;
    file_1089 << "\n";

    file_1089 << "Givens: \n";
    //file_1089 << givens_error;
    file_1089 << "\n";

    file_1089 << "Mínimos cuadrados: \n";
    file_1089 << leastSquares_error;
    file_1089 << "\n";

    file_1089 << "--------------- Tiempos --------------- \n\n";

    std::ostringstream secs_seidel;
    std::ostringstream nano_seidel;
    std::string secs_seidelstr;
    std::string nano_seidelstr;

    secs_seidel << diff(time1_seidel,time2_seidel).tv_sec;
    secs_seidelstr = secs_seidel.str();

    nano_seidel << diff(time1_seidel,time2_seidel).tv_nsec;
    nano_seidelstr = nano_seidel.str();

    file_1089 << "Seidel: \n";
    file_1089 << secs_seidelstr + " : " + nano_seidelstr;
    file_1089 << "\n";


    std::ostringstream secs_qr;
    std::ostringstream nano_qr;
    std::string secs_qrstr;
    std::string nano_qrstr;

    secs_qr << diff(time1_qr,time2_qr).tv_sec;
    secs_qrstr = secs_qr.str();

    nano_qr << diff(time1_qr,time2_qr).tv_nsec;
    nano_qrstr = nano_qr.str();

  
    file_1089 << "QR: \n";
    file_1089 << secs_qrstr + " : " + nano_qrstr;
    file_1089 << "\n\n";

    std::ostringstream secs_chol;
    std::ostringstream nano_chol;
    std::string secs_cholstr;
    std::string nano_cholstr;

    secs_chol << diff(time1_chol,time2_chol).tv_sec;
    secs_cholstr = secs_chol.str();

    nano_chol << diff(time1_chol,time2_chol).tv_nsec;
    nano_cholstr = nano_chol.str();

    file_1089 << "Cholesky: \n";
    file_1089 << secs_cholstr + " : " + nano_cholstr;
    file_1089 << "\n\n";

    std::ostringstream secs_givens;
    std::ostringstream nano_givens;
    std::string secs_givensstr;
    std::string nano_givensstr;
    
    secs_givens << diff(time1_givens,time2_givens).tv_sec;
    secs_givensstr = secs_givens.str();

    nano_givens << diff(time1_givens,time2_givens).tv_nsec;
    nano_givensstr = nano_givens.str();

    file_1089 << "Givens: \n";
    file_1089 << secs_givensstr + " : " + nano_givensstr;
    file_1089 << "\n\n";

    std::ostringstream secs_ls;
    std::ostringstream nano_ls;
    std::string secs_lsstr;
    std::string nano_lsstr;

    secs_ls << diff(time1_leastSquares,time2_leastSquares).tv_sec;
    secs_lsstr = secs_ls.str();

    nano_ls << diff(time1_leastSquares,time2_leastSquares).tv_nsec;
    nano_lsstr = nano_ls.str();

    file_1089 << "Mínimos cuadrados: \n";
    file_1089 << secs_lsstr + " : " + nano_lsstr;
    file_1089 << "\n\n";

    file_1089.close(); 

    /*
    ofstream file_4225("4225x4225results.txt");

    file_4225 << "---------------- Errores ------------------ \n\n";

    file_4225 << "Seidel: \n";
    file_4225 << seidel_errors;
    file_4225 << "\n";

    file_4225 << "QR: \n";
    file_4225 << QR_error;
    file_4225 << "\n";

    file_4225 << "Cholesky: \n";
    file_4225 << chol_error;
    file_4225 << "\n";

    file_4225 << "Givens: \n";
    file_4225 << '0';
    file_4225 << "\n";

    file_4225 << "Mínimos cuadrados: \n";
    file_4225 << leastSquares_error;
    file_4225 << "\n";

    file_4225 << "--------------- Tiempos --------------- \n\n";

    std::ostringstream secs_seidel;
    std::ostringstream nano_seidel;
    std::string secs_seidelstr;
    std::string nano_seidelstr;

    secs_seidel << diff(time1_seidel,time2_seidel).tv_sec;
    secs_seidelstr = secs_seidel.str();

    nano_seidel << diff(time1_seidel,time2_seidel).tv_nsec;
    nano_seidelstr = nano_seidel.str();

    file_4225 << "Seidel: \n";
    file_4225 << secs_seidelstr + " : " + nano_seidelstr;
    file_4225 << "\n";


    std::ostringstream secs_qr;
    std::ostringstream nano_qr;
    std::string secs_qrstr;
    std::string nano_qrstr;

    secs_qr << diff(time1_qr,time2_qr).tv_sec;
    secs_qrstr = secs_qr.str();

    nano_qr << diff(time1_qr,time2_qr).tv_nsec;
    nano_qrstr = nano_qr.str();

  
    file_4225 << "QR: \n";
    file_4225 << secs_qrstr + " : " + nano_qrstr;
    file_4225 << "\n\n";

    std::ostringstream secs_chol;
    std::ostringstream nano_chol;
    std::string secs_cholstr;
    std::string nano_cholstr;

    secs_chol << diff(time1_chol,time2_chol).tv_sec;
    secs_cholstr = secs_chol.str();

    nano_chol << diff(time1_chol,time2_chol).tv_nsec;
    nano_cholstr = nano_chol.str();

    file_4225 << "Cholesky: \n";
    file_4225 << secs_cholstr + " : " + nano_cholstr;
    file_4225 << "\n\n";

    std::ostringstream secs_givens;
    std::ostringstream nano_givens;
    std::string secs_givensstr;
    std::string nano_givensstr;
    
    secs_givens << diff(time1_givens,time2_givens).tv_sec;
    secs_givensstr = secs_givens.str();

    nano_givens << diff(time1_givens,time2_givens).tv_nsec;
    nano_givensstr = nano_givens.str();

    file_4225 << "Givens: \n";
    file_4225 << secs_givensstr + " : " + nano_givensstr;
    file_4225 << "\n\n";

    std::ostringstream secs_ls;
    std::ostringstream nano_ls;
    std::string secs_lsstr;
    std::string nano_lsstr;

    secs_ls << diff(time1_leastSquares,time2_leastSquares).tv_sec;
    secs_lsstr = secs_ls.str();

    nano_ls << diff(time1_leastSquares,time2_leastSquares).tv_nsec;
    nano_lsstr = nano_ls.str();

    file_4225 << "Mínimos cuadrados: \n";
    file_4225 << secs_lsstr + " : " + nano_lsstr;
    file_4225 << "\n\n";


    file_4225.close();
    */
   
  return 0;
  }