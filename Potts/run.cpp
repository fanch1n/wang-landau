#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>



/* Macro for transformation energy level number -> energy */
#define energy(b) (4.0*(b)-2.0*L*L)

/* Macro for reading parameters */
#define MAX_LEN 256
#define read_variable(string,type,variable)   \
 {                       \
  char input[MAX_LEN];            \
  fputs(string,stdout);            \
  fgets(input,MAX_LEN,stdin);         \
  sscanf(input,type,&variable);        \
 }

using namespace std; 
int main(int argc,char ** argv) 
{
  int L;   /* System size (L by L sites) */
  int top_b; /* Number of energy levels */
  double flat_thres;  /* Threshold for flat histogram */
  int b;   /* Current energy level number */
  int b_new; /* Energy Level number of proposed move */
  int seed;  /* Random number seed */
  int skip;  /* Number of moves to skip between checks of histogram */
  int min_steps;  /* Minimum number of moves for each f */
  double min_f;  /* f to start with */
  double dec_pow; /* Exponent of power law to decrease f after each run */
  double prob;/* Probability for acceptance of move */
  int mc_steps; /* Counters */
  char filename[MAX_LEN];
  FILE * out_file;

  L = 8; 
  min_f = 1.0 + 1E-4;
  dec_pow = 1.0/2.0;  
  flat_thres = 0.8; 
  min_steps = 1E5; 
  skip = 100; 

  seed = 1; 
  srand(seed);

  vector<double> g(L*L+1); 
  for(int i = 0; i < L*L + 1; i++)
  {
    g[i] = 0; 
  }
  vector<int> hist(L*L+1, 0); 
  vector<vector<int> > s(L, vector<int>(L)); 

  /* Form initial state */
  for(int i=0;i<L;i++)
  {
    for(int j=0;j<L;j++)
    {
      s[i][j] = 1;
    }
  }

  /* set upper limit for energy level number; depends on if L is odd or even */
  if(L%2==1)
   top_b=L*(L-1)+1; //odd
  else
   top_b=L*L+1; //even

  /* Set energy level number of initial state */
  b=0;


  double f = 2.71;  /* Factor to multiply g(E) after move is accepted */
  /* Repeat loop until f reaches minimum value */
  while(f>min_f) 
  {
   /* Only use logarithms of f and g for the calculation */
   double lnf=log(f);
   int c, cont=1;

   /* Start with flat histogram */
   for(int i = 0; i <= L*L; i++) {
    hist[i]=0;
   }

   /* Reset counter for number of steps */
   int n=0;
   c=skip+1;
   mc_steps=0;
   do {
    
    /* Perform one Monte Carlo sweep: randomly pick N=L*L sites */
    for(int i=0;i<L*L;i++) {
     
     int si,sj,neighsum;

     /* Choose random site for flipping */
     int site=(int)(L*L*(rand()/(RAND_MAX+1.0)));
     si=site/L;
     sj=site%L;
  
     /* Calculate sum over neighbor spins */
     neighsum= s[(si!=L-1)?si+1:0][sj];
     neighsum+=s[(si!=0)?si-1:L-1][sj];
     neighsum+=s[si][(sj!=L-1)?sj+1:0];
     neighsum+=s[si][(sj!=0)?sj-1:L-1];

     /* Number of new energy level */
     b_new=b+s[si][sj]*neighsum/2;

     /* Check if energy level is inside specified range */
     if(b_new<top_b) {

      /* Calculate probability for acceptance */
      prob=exp(g[b]-g[b_new]);

      /* Accept if proposed move has lower g (prob>=1.0) */
      /* or if random number [0.0,1.0] is less than prob */
      if((prob>=1.0) || ((double)rand()/RAND_MAX < prob)) {
       b=b_new;
       s[si][sj]*=-1;
      }

      /* Multiply g[b] by f: add logarithms */
      g[b]+=lnf;

      /* Update counters */
      ++hist[b];
      n++;
     }
    }

    /* Update counters for Monte Carlo steps */
    mc_steps++;
    c++;


    /* Check for flat histogram after skipping skip steps */
    if((mc_steps>=min_steps) && (c>=skip)) {

     int div;

     /* Reset skip counter */
     c=0;

     /* Check for flat histogram */
     /* Exclude "empty" energy levels, energy jumps by 8J from groundstate to first excited state */
     div=top_b>L*L-1?top_b-2:top_b-1;

     /* Look into each energy level and continue if fraction is larger than given threshold */
     int accum = 0; 
     for(int i=0; (i<top_b) && ((i==1) || (i==L*L-1) || ((double)hist[i]/n*div > flat_thres)); i++)
     {
       accum++; 
     }

     /* If loop made it to the end, histogram is flat */
     if(accum==top_b) cont=0;
     else cont=1;
    }
    else cont=1;
   } while(cont == 1);

   /* Normalize (logarithmic) g values, so g(0)=1 */
   double g0 = g[0]; 
   for(int i=0; i<top_b;i++) {
    g[i] -= g0;
   }

   /* Give status message with current f */
   printf("f-1=%.3e\tlog(f)=%g\tmc_steps=%d\n",f-1.,lnf,mc_steps);

   /* Decrease factor f by a power law */
   f=pow(f,dec_pow);
  }

  ofstream wl_fs("wl.out"); 

  wl_fs << "b\tE(b)\tlog(g(E))\tH(E)\n" << endl;
  wl_fs << "-----------------------------------------\n" << endl;
  cout << "top_b: " << top_b << endl; 
  cout << "g.size: " << g.size() << endl; 
  cout << "hist.size: " << hist.size() << endl; 
  for(int i=0;i<top_b;i++) 
  {
   if((i!=1)&&(i!=L*L-1))
   {
     wl_fs << i << "\t" << energy(i) << "\t" << g[i]+log(2.0) << "\t" << hist[i] << endl; 
   }
  }

  wl_fs.close(); 

  return 0;
}

