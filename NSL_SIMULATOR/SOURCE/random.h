/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

class Random
{

private:
    const int m1{502};
    const int m2{1521};
    const int m3{4071};
    const int m4{2107};
    int l1, l2, l3, l4, n1, n2, n3, n4;

protected:
public:
    // constructors
    Random();
    // destructor
    ~Random();
    // methods
    void SetRandom(int *, int, int);
    void SaveSeed();
    double Rannyu(void);
    double Rannyu(double min, double max);
    double Gauss(double mean, double sigma);
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
