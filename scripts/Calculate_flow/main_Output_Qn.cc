#include<iostream>
#include <math.h>
#include<fstream>
#include<string>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include <vector>
#include <algorithm>

// *** The code to calculate observables output from jat+parton cascade framework *** //
// Copyright @  Wenbin Zhao, 2023 //

using namespace std;
int main(int argc, char* argv[] )
{
    int pid, total_number_of_particles;
    // Output file
    string output_filename2;
    output_filename2 = "QAQB.dat";// output files of final hadrons
    string output_filenamench;
    output_filenamench = "Nch.dat";// output files of final hadrons

    // Read the binary output
    stringstream strparticle_name;
    strparticle_name << "final_state_hard_hadrons.bin";
    std::ifstream InStream;
    InStream.precision(15);
    string particle_name = strparticle_name.str();
    InStream.open(particle_name.c_str(), std::ios::in | std::ios::binary);
    double mass, pT, phi, eta;
    float temp;
    double QnAB_0p0_3[11][4] = {0.0}; double QnAB_0p3_3[11][4] = {0.0}; double QnAB_0p5_3[11][4] = {0.0};
    //int MAMB_0_3[10] = {0}; int MAMB_0p3_3[10] = {0}; int MAMB_0p5_3[10] = {0}; 
    std::vector<int> NchVector; 
    double Nchbins[12] = {-1., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 500.,}; 
    while( InStream.read(reinterpret_cast<char*>(&total_number_of_particles), sizeof(int))) {
        int Nchtemp = 0;
        double QnA_0p0_3_real[4] = {0.}; double QnA_0p3_3_real[4] = {0.}; double QnA_0p5_3_real[4] = {0.}; 
        double QnB_0p0_3_real[4] = {0.}; double QnB_0p3_3_real[4] = {0.}; double QnB_0p5_3_real[4] = {0.}; 
        double QnA_0p0_3_imag[4] = {0.}; double QnA_0p3_3_imag[4] = {0.}; double QnA_0p5_3_imag[4] = {0.}; 
        double QnB_0p0_3_imag[4] = {0.}; double QnB_0p3_3_imag[4] = {0.}; double QnB_0p5_3_imag[4] = {0.}; 
        // Start one event
        for (auto i=0; i<total_number_of_particles; i++) {
            InStream.read(reinterpret_cast<char*>(&pid), sizeof(int));
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            mass = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            pT = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            phi = temp;
            InStream.read(reinterpret_cast<char*>(&temp), sizeof(float));
            eta = temp;
            //cout <<total_number_of_particles << " " <<  pid << "  " << mass << "  " << pT << "  " << phi << "  " << eta << endl;
            if (abs(pid)==211 || abs(pid) == 321 || abs(pid) == 2212) {
                Nchtemp ++;
                if (pT < 3.0) {
                    // pT > 0.5
                    if (pT > 0.5) {
                        for (auto iorder = 0; iorder < 4; iorder++) {
                            if (eta > 1.0) {
                                QnA_0p5_3_real[iorder] = QnA_0p5_3_real[iorder] + cos(iorder * 1. * phi);
                                QnA_0p5_3_imag[iorder] = QnA_0p5_3_imag[iorder] + sin(iorder * 1. * phi);
                            }
                            if (eta < -1.0) {
                                QnB_0p5_3_real[iorder] = QnB_0p5_3_real[iorder] + cos(iorder * 1. * phi);
                                QnB_0p5_3_imag[iorder] = QnB_0p5_3_imag[iorder] + sin(iorder * 1. * phi);
                            }
                        }
                    }
                    // pT > 0.3
                    if (pT > 0.3) {
                        for (auto iorder = 0; iorder < 4; iorder++) {
                            if (eta > 1.0) {
                                QnA_0p3_3_real[iorder] = QnA_0p3_3_real[iorder] + cos(iorder * 1. * phi);
                                QnA_0p3_3_imag[iorder] = QnA_0p3_3_imag[iorder] + sin(iorder * 1. * phi);
                            }
                            if (eta < -1.0) {
                                QnB_0p3_3_real[iorder] = QnB_0p3_3_real[iorder] + cos(iorder * 1. * phi);
                                QnB_0p3_3_imag[iorder] = QnB_0p3_3_imag[iorder] + sin(iorder * 1. * phi);
                            }
                        }
                    }
                    // pT > 0.0
                        for (auto iorder = 0; iorder < 4; iorder++) {
                            if (eta > 1.0) {
                                QnA_0p0_3_real[iorder] = QnA_0p0_3_real[iorder] + cos(iorder * 1. * phi);
                                QnA_0p0_3_imag[iorder] = QnA_0p0_3_imag[iorder] + sin(iorder * 1. * phi);
                            }
                            if (eta < -1.0) {
                                QnB_0p0_3_real[iorder] = QnB_0p0_3_real[iorder] + cos(iorder * 1. * phi);
                                QnB_0p0_3_imag[iorder] = QnB_0p0_3_imag[iorder] + sin(iorder * 1. * phi);
                            }
                        }
                }
            }
        }
        NchVector.push_back(Nchtemp);
        // Then Calculate the QAQB
        for (auto inch=0; inch<11; inch++) {
            if (Nchtemp > Nchbins[inch] && Nchtemp <= Nchbins[inch+1]) {
                for (int iorder=0; iorder<4; iorder++) {
                    QnAB_0p0_3[inch][iorder] = QnA_0p0_3_real[iorder] * QnA_0p0_3_real[iorder] 
                                             + QnA_0p0_3_imag[iorder] * QnA_0p0_3_imag[iorder];
                    QnAB_0p3_3[inch][iorder] = QnA_0p3_3_real[iorder] * QnA_0p3_3_real[iorder] 
                                             + QnA_0p3_3_imag[iorder] * QnA_0p3_3_imag[iorder];
                    QnAB_0p5_3[inch][iorder] = QnA_0p5_3_real[iorder] * QnA_0p5_3_real[iorder] 
                                             + QnA_0p5_3_imag[iorder] * QnA_0p5_3_imag[iorder];
                }
            }
        }
    }
    // Then output the Qn and NchVector

    ofstream output2(output_filename2.c_str()); 
    if (!output2.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filename2 << endl;
        return -1;
    }
    output2 << " # Qn ( n = 0,1,2,3) of 0.0<pT<3.0, 0.3<pT<3.0, 0.5<pT<3.0 " << endl;
    for (auto inch=0; inch<11; inch++) {
        for (int iorder=0; iorder<4; iorder++) {
            output2 << QnAB_0p0_3[inch][iorder] << "  "
                    << QnAB_0p3_3[inch][iorder] << "  "
                    << QnAB_0p5_3[inch][iorder] << "  ";
        }
        output2 << endl;
    }
        
    ofstream outputnch(output_filenamench.c_str()); 
    if (!outputnch.is_open() ) {
        cout << "cannot open output file:"<< endl
         << output_filenamench << endl;
        return -1;
    }
    for (int ii=0; ii<NchVector.size(); ii++) {
        outputnch << NchVector[ii] << endl;
    }
    
    output2.close();
    outputnch.close();
    InStream.close();
}

