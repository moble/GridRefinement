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
//int Revision = atoi(Rev.substr(5,Rev.size()-6).c_str());

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

void Interpolate(vector<double>& t,
                 vector<double>& amp,
                 vector<double>& phi,
                 const double dt);

double mindiff(const vector<double>& t) {
  if(t.size()<2) { return 1.0; }
  double dt = fabs(t[1]-t[0]);
  double dt2 = dt;
  for(unsigned int i=2; i<t.size(); ++i) {
    dt2 = fabs(t[i]-t[i-1]);
    dt = ((dt2<dt) ? dt2 : dt);
  }
  return dt;
}

int main(int argc, const char* argv[]) {
  // Get input options
  if(strcmp(argv[1], "-h")==0 || strcmp(argv[1], "--help")==0) {
    cerr << "usage:\n"
         << "  Reconstitute <FileName.minimal> [dt]\n"
         << "  Reconstitute --help\n\n"
         << "The output file is 'FileName' -- that is, '.minimal' is removed.\n"
         << "If the input FileName does not end in '.minimal', the output\n"
         << "file name is 'FileName.reconstituted'.\n\n"
         << "If no time step is given, dt is taken to be the smallest time\n"
         << "step in the input data." << endl;
    return 0;
  }
  string CommandLine = argv[0];
  for(int i=1; i<argc; ++i) {
    CommandLine = CommandLine + " " + argv[i];
  }
  string InFileName=argv[1], OutFileName=InFileName;
  size_t MinPos = InFileName.find(".minimal");
  if(MinPos != string::npos) {
    OutFileName = InFileName.substr(0, MinPos);
  } else {
    OutFileName = InFileName + ".reconstituted";
  }
  double dt = 0.0;
  if(argc>2) {
    dt = atof(argv[2]);
  }
  
  // Initialize data to be read from input file
  string Header("");
  vector<double> t(0), amp(0), phi(0);
  
  // Read data file into a vectors with time, amp, and phi
  cout << "Reading dat file..." << endl;
  ReadDatFile(InFileName, Header, t, amp, phi);
  if(dt==0.0) { dt = mindiff(t); }
  
  // Interpolate to the finer grid
  cout << "Interpolating to finer grid with dt=" << dt << " ..." << endl;
  Interpolate(t, amp, phi, dt);
  
  // Output data file
  cout << "Writing dat file..." << endl;
  string Comment = Header.substr(0,1);
  if(Comment.empty()) { Comment = "#"; }
  Header = Comment + Comment + " Reconstituted with `"
    + CommandLine + "` under svn Rev " + Revision + ".\n" + Header
    + Comment + " time\t\thplus\t\thcross\n";
  WriteDatFile(OutFileName, Header, t, amp, phi);
  
  cout << "Finished." << endl;
  
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
  
  // Initialize variables
  string Temp="";
  double tVal=0.0, ampVal=0.0, phiVal=0.0;
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
      InFile >> ampVal;
      InFile >> phiVal;
      InFile >> tVal;
      amp.push_back(ampVal);
      phi.push_back(phiVal);
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
      OutFile << t[i] << " " << amp[i]*cos(phi[i]) << " " << amp[i]*sin(phi[i])
              << endl;
    }
    OutFile.close();
  } else {
    cerr << "Couldn't open output file " << OutFileName
         << " for writing." << endl;
    exit(1);
  }
  return;
}


// Simple linear interpolation
void Interpolate(vector<double>& t,
                 vector<double>& amp,
                 vector<double>& phi,
                 const double dt)
{
  int N = int(floor((t[t.size()-1]-t[0])/dt));
  vector<double> T(N, 0.0);
  vector<double> Amp(N, 0.0);
  vector<double> Phi(N, 0.0);
  int i=0;
  for(int I=0; I<N; ++I) {
    T[I] = t[0] + I*dt;
    while(t[i+1]<T[I]) { ++i; }
    Amp[I] = amp[i] + (T[I]-t[i]) * (amp[i+1]-amp[i]) / (t[i+1]-t[i]);
    Phi[I] = phi[i] + (T[I]-t[i]) * (phi[i+1]-phi[i]) / (t[i+1]-t[i]);
  }
  t   = T;
  amp = Amp;
  phi = Phi;
  return;
}
