#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

string Rev = "$Rev: 128 $";
string Revision = Rev.substr(6,Rev.size()-8).c_str();

void ReadDatFile(const string& InFileName,
                 string& Header,
                 vector<double>& t,
                 vector<double>& amp,
                 vector<double>& phi);

void WriteDatFile(const string& OutFileName,
                  const string& Header,
                  const vector<double>& t,
                  const vector<double>& amp,
                  const vector<double>& phi);

vector<double> Unwrap(const vector<double>& In);

void MinimalGrid(const double ampTol,
                 const double phiTol,
                 vector<double>& t,
                 vector<double>& amp,
                 vector<double>& phi);


int main(int argc, const char* argv[]) {
  // Get input options
  if(argc<2 || strcmp(argv[1], "-h")==0 || strcmp(argv[1], "--help")==0) {
    cout << "usage:\n"
         << "  MinimizeGrid <FileName> [AmpTol [PhiTol]]\n"
         << "  MinimizeGrid --help\n\n"
         << "The output file name is 'FileName.minimal' and has data\n"
         << "such that the original data can be recovered to within\n"
         << "the given tolerances by linear interpolation (using, \n"
         << "e.g., the 'ReconstituteGrid' routine).  If no tolerance\n"
         << "is given, 1e-5 is assumed for both amplitude and phase;\n"
         << "if only one tolerance is given, it is used for both;\n"
         << "if two tolerances are given, the first is used for\n"
         << "amplitude and the second for phase.\n\n"
         << "Note that the tolerances are recorded as comments in\n"
         << "the output header."
         << endl;
    return 0;
  }
  string CommandLine = argv[0];
  for(int i=1; i<argc; ++i) {
    CommandLine = CommandLine + " " + argv[i];
  }
  string InFileName=argv[1], OutFileName=InFileName+".minimal";
  double ampTol=1.0e-5, phiTol=ampTol;
  if(argc>2) {
    ampTol = atof(argv[2]);
    cout << "ampTol = " << ampTol << endl;
    if(argc>3) {
      phiTol = atof(argv[3]);
    } else {
      phiTol = ampTol;
    }
    cout << "phiTol = " << phiTol << endl;
  }
  
  // Initialize data to be read from input file
  string Header("");
  vector<double> t(0), amp(0), phi(0);
  
  // Read data file into a vectors with time, amp, and phi
  cout << "Reading dat file..." << endl;
  ReadDatFile(InFileName, Header, t, amp, phi);
  int InitialLength = t.size();
  
  // Unwrap phi
  cout << "Unwrapping phi..." << endl;
  phi = Unwrap(phi);
  
  // Remove unecessary points
  cout << "Minimizing grid..." << endl;
  MinimalGrid(ampTol, phiTol, t, amp, phi);
  
  // Output data file
  cout << "Writing dat file..." << endl;
  string Comment = Header.substr(0,1);
  if(Comment.empty()) { Comment = "#"; }
  Header = Comment + Comment + " Minimized with `"
    + CommandLine + "` under svn Rev " + Revision + ".\n" + Header
    + Comment + " time\t\t   amp \t\t\tphi\n";
  WriteDatFile(OutFileName, Header, t, amp, phi);
  
  cout << "Finished.  Reduced waveform to " << t.size() << " points from "
       << InitialLength << " points -- a " << setprecision(3)
       << 100.0*(1.0 - double(t.size())/double(InitialLength))
       << "% reduction." << endl;
  
  return 0;
}




void ReadDatFile(const string& InFileName,
                 string& Header,
                 vector<double>& t,
                 vector<double>& amp,
                 vector<double>& phi)
{
  // Make sure the string and vectors are empty
  Header = "";
  t   = vector<double>(0);
  amp = vector<double>(0);
  phi = vector<double>(0);
  
  // Initialize variable
  string Temp="";
  double tVal=0.0, reVal=0.0, imVal=0.0;
  ifstream InFile(InFileName.c_str());
  
  // Read the file
  if(InFile.is_open()) {
    // Get the header
    while(InFile.peek()=='#' || InFile.peek()=='%') {
      getline(InFile, Temp);
      if(Temp.compare(1, 5, " time") != 0 || InFile.peek()=='#' || InFile.peek()=='%') {
        Header += Temp + "\n";
      }
    }
    // Get the data
    InFile >> tVal;
    while (!InFile.eof()) {
      t.push_back(tVal);
      InFile >> reVal;
      InFile >> imVal;
      InFile >> tVal;
      amp.push_back(sqrt(reVal*reVal + imVal*imVal));
      phi.push_back(atan2(imVal, reVal));
    }
    InFile.close();
  } else {
    cerr << "Couldn't open input file " << InFileName
         << " for reading."  << endl;
    exit(1);
  }
  return;
}


void WriteDatFile(const string& OutFileName,
                  const string& Header,
                  const vector<double>& t,
                  const vector<double>& amp,
                  const vector<double>& phi)
{
  ofstream OutFile(OutFileName.c_str());
  if(OutFile.is_open()) {
    OutFile << Header << flush;
    OutFile << setprecision(15);
    for(unsigned int i=0; i<t.size(); ++i) {
      OutFile << t[i] << " " << amp[i] << " " << phi[i] << endl;
    }
    OutFile.close();
  } else {
    cerr << "Couldn't open output file " << OutFileName
         << " for writing." << endl;
    exit(1);
  }
  return;
}


// Unwrap the phase data so that it's as smooth as possible.
// Compare Matlab's unwrap.m file.
vector<double> Unwrap(const vector<double>& In) {
  vector<double> Out = In;
  double Dp = 0.0;
  double Dps = 0.0;
  double CumCorr = 0.0;
  
  // Dp will contain the incremental phase variations;
  // Dps will contain the equivalents, confined to [-pi,pi)
  // CumCorr will contain the incremental phase corrections
  for(unsigned int i=1; i<In.size(); ++i) {
    Dp = In[i]-In[i-1];
    // C++'s fmod is unlike Matlab's negative 'a' values, so:
    if(Dp+M_PI<0) {
      Dps = M_PI - fmod(-Dp-M_PI, 2.0*M_PI);
    } else {
      Dps = fmod(Dp+M_PI, 2.0*M_PI) - M_PI;
    }
    if(Dps==-M_PI && Dp>0) { Dps = M_PI; }
    CumCorr += Dps - Dp;
    Out[i] += CumCorr;
  }
  
  return Out;
}


// Just a warning message if the minimum time step is being used
bool Cautioned = false;
void Caution(const double ampTol, const double phiTol=0.0) {
  if(!Cautioned) {
    Cautioned = true;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
         << "                 CAUTION!                 \n"
         << " The minimum timestep in the input data   \n"
         << " is being used.  This suggests that there \n"
         << " may be features so sharp that            \n"
         << " interpolation will not be capable of     \n"
         << " achieving the requested tolerances of    \n"
         << " AmpTol=" << ampTol << " and PhiTol="
         << (phiTol==0.0 ? ampTol : phiTol) << ", \n"
         << " except on precisely the original grid.   \n"
         << " Of course, this problem would probably   \n"
         << " be present even without trying to use a  \n"
         << " minimal grid, so you probably don't need \n"
         << " to worry at this step.  However, you     \n"
         << " might consider filtering the input data  \n"
         << " as this is typically a noise issue.      \n"
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
         << endl;
  }
  return;
}


// Returns true if the input (x,y) data can be interpolated
// between indices I0 and I1 to within a tolerance of Tol
// at the midpoint between I0 and I1.
bool Check(const vector<double>& x, const vector<double>& y,
           const int I0, const int I1, const double Tol)
{
  if(I1==I0+1) { return true; }
  int IMid = ((I0+I1) >> 1); // = floor(avg(I0+I1))
  if(abs(y[I0] + (x[IMid]-x[I0])*(y[I1]-y[I0])/(x[I1]-x[I0]) - y[IMid]) > Tol) {
    return false; // Interval is not fine enough
  }
  return true; // Interval is fine enough
}


// Given data (t,phi), some initial index I0, and a guess
// for I1, this function outputs the optimal index I1 so
// that (t,phi) can be interpolated between I0 and I1 to
// just within phiTol of the full input data set.
// Compare Numerical Recipes's 'hunt' function; this is
// basically a hunt for that optimal index.
int Hunt(const vector<double>& t, const vector<double>& phi,
         const double phiTol, const int I0, const int I1)
{
  int Inc=1, I1lo = I1, I1hi=I1+1;
  
  // Bracket the optimal I1 between I1lo and I1hi
  if( Check(t, phi, I0, I1lo, phiTol) ) {
    while( Check(t, phi, I0, I1hi, phiTol) ) {
      I1lo = I1hi;
      Inc *= 2;
      I1hi += Inc;
      if(I1hi>int(t.size())) {
        I1hi = t.size();
        break;
      }
    }
  } else {
    if(I1lo<=I0+2) {
      return I0+2;
    }
    I1hi = I1lo;
    I1lo -= 1;
    while( ! Check(t, phi, I0, I1lo, phiTol) ) {
      I1hi = I1lo;
      Inc *= 2;
      if(Inc > I1hi) {
        I1lo = I0+2;
        break;
      } else {
        I1lo = I1hi-Inc;
      }
    }
  }
  
  // Now use bisection between I1lo and I1hi
  while(I1hi-I1lo != 1) {
    int I1m=((I1hi+I1lo)>>1);
    if( Check(t, phi, I0, I1m, phiTol) ) {
      I1lo = I1m;
    } else {
      I1hi = I1m;
    }
  }
  
  return I1lo;
}


// This function does the actual work of selecting a minimal
// grid from the input data.
void MinimalGrid(const double ampTol, const double phiTol,
                 vector<double>& t,
                 vector<double>& amp,
                 vector<double>& phi)
{
  // The objective here will be to create a vector of bool's, the same
  // length as t.  The truth value will correspond to whether or not
  // that time step should be kept in the final data.  We begin by
  // assuming that the very first and last steps should obviously be
  // kept.  Then, there are two stages.  First is a coarse stage,
  // which steps through the data making intervals small enough to
  // reproduce the phi data at the interval's midpoint to within
  // phiTol, but no smaller.  Second is the finer stage, which goes
  // through each interval, checking that every single point in the
  // input data can be reproduced to within phiTol and ampTol.  If
  // that's not true, the interval is split evenly into two, and the
  // algorithm proceeds with the earlier interval.  Finally, the input
  // t, amp, and phi vectors are replaced by the smaller vectors given
  // by our vector of bool's.
  
  unsigned int I0 = 0;
  unsigned int I1 = ((t.size()-1) >> 1); // = midpoint of the input data set
  unsigned int NumPoints = 2;
  vector<bool> T(t.size(), false);
  T[0] = true;
  T[t.size()-1] = true;
  
  // Coarse -- check only phi at midpoints of each interval
  //   This loop starts from the beginning of the data set, and
  //   forms the smallest interval such that the phi Tolerance
  //   is achieved by linear interpolation.  Then, it moves to
  //   the end of that interval to find the next interval, etc.
  while(((I0+I1)>>1) < T.size()-1) {
    // hunt for optimal I1
    I1 = Hunt(t, phi, phiTol, I0, I1);
    
    if(!T[I1]) {
      T[I1] = true;
      ++NumPoints;
    }
    I0 = I1;
    I1 = 2*I1 - I0;
    if(I1<I0+2) { I1 = I0+2; }
  }
  
  // Fine -- check amp and phi at every point
  //   This loop goes through each of the intervals found above,
  //   and makes sure that every data point in both phi and amp
  //   can be reconstructed to within the given tolerances.  If
  //   not, it just adds the midpoint of the interval, and
  //   checks the new interval again.
  I0 = 0;
  I1 = 1;
  unsigned int i=1;
  while(i<t.size()) {
    if(i>I1) { // This could happen below
      while(i>I1) { ++I1; }
      I0 = i;
      while(!T[I0]) { --I0; }
    }
    while(!T[I1]) { ++I1; }
    if(i != I0 && i != I1) {
      if(abs(1-(amp[I0]+(t[i]-t[I0])*(amp[I1]-amp[I0])/(t[I1]-t[I0]))/amp[i]) > ampTol
         || abs(phi[I0]+(t[i]-t[I0])*(phi[I1]-phi[I0])/(t[I1]-t[I0])-phi[i]) > phiTol) {
        I1 = ((I0+I1)>>1);
        if(!T[I1]) {
          if(I1==I0+1) { Caution(ampTol, phiTol); }
          T[I1] = true;
          ++NumPoints;
        }
        continue;
      }
    }
    if(i==I1) {
      I0 = I1;
      ++I1;
    }
    ++i;
  }
  
  // Take only the smaller grid
  vector<double> tOut(NumPoints, 0.0);
  vector<double> ampOut(NumPoints, 0.0);
  vector<double> phiOut(NumPoints, 0.0);
  int Point = 0;
  for(unsigned int i=0; i<T.size(); ++i) {
    if(T[i]) {
      tOut[Point]   = t[i];
      ampOut[Point] = amp[i];
      phiOut[Point] = phi[i];
      ++Point;
    }
  }
  t   = tOut;
  amp = ampOut;
  phi = phiOut;
  
  return;
}
