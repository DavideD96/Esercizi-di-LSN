/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"
#include <iomanip>

using namespace std;

//ATTENZIONE CONTROLLARE SE IL DATA BLOCKING è IMPLEMENTATO CORRETTAMENTE.

int main(){

  Input();             //Inizialization
  int nconf = 1;
  //double dist_prova;
  //int nbin_prova = 100;
  //int bin_prova;
  //double size_prova = 0.025;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
       if(nmeasure%nInBlocks == 0 && nmeasure != 0){
	  Blocking();
       }
       //dist_prova = 1.091;
       //bin_prova = binSearcher(size_prova, nbin_prova, dist_prova);
       //cout << "bin prova = " << bin_prova <<endl;
     } 
  }
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> oldOrNot; //*

  ReadInput >> temp;

  ReadInput >> npart;

  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nInBlocks;

  nBlocks = (nstep/10)/nInBlocks;

  cout << "The program integrates Newton equations with the Verlet method " << endl;

  if(oldOrNot == 1){
  	cout << "The program starts the simulation using the configurations of the paricles in two consecutive instants." << endl; 
  }else{
	cout << "The program starts the simulation using an initial configuration and sorting the initial velocities." << endl;
  }

  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

  for(int i=0; i<nbin; i++){
    gr[i] = 0;
    gr_sqr[i] = 0;
    gr_mean[i] = 0;
  }

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4+nbin; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    //cout<<"x[i] = "<<x[i]<<"   "<<"y[i] = "<<y[i]<<"   "<<"z[i] = "<<z[i]<<endl;	
    //cout<<endl;   
  }
  ReadConf.close();

  if(oldOrNot == 1){

//Read old configuration
  	cout << "Read old configuration from file old.0 " << endl << endl;
  	ReadConf.open("old.0");
  	for (int i=0; i<npart; ++i){
    	  ReadConf >> xold[i] >> yold[i] >> zold[i];
    	  xold[i] = xold[i] * box;
    	  yold[i] = yold[i] * box;
    	  zold[i] = zold[i] * box;

	  //cout<<"xold[i] = "<<xold[i]<<"   "<<"yold[i] = "<<yold[i]<<"   "<<"zold[i] = "<<zold[i]<<endl;	
  
  	}
  	ReadConf.close();

  	double xnew_[m_part], ynew_[m_part], znew_[m_part], fx_[m_part], fy_[m_part], fz_[m_part], stima_kin_, stima_temp_;

  	for(int i=0; i<npart; ++i){ //Force acting on particle i
  	  fx_[i] = Force(i,0);
  	  fy_[i] = Force(i,1);
 	  fz_[i] = Force(i,2);
  	}

  	for(int i=0; i<npart; ++i){ //Verlet integration scheme

    	  xnew_[i] = Pbc( 2.0 * x[i] - xold[i] + fx_[i] * pow(delta,2) );
    	  ynew_[i] = Pbc( 2.0 * y[i] - yold[i] + fy_[i] * pow(delta,2) );
    	  znew_[i] = Pbc( 2.0 * z[i] - zold[i] + fz_[i] * pow(delta,2) );
 
    	  vx[i] = Pbc(xnew_[i] - x[i])/delta; //compute v(t+dt/2)
    	  vy[i] = Pbc(ynew_[i] - y[i])/delta;
    	  vz[i] = Pbc(znew_[i] - z[i])/delta;
	}

  	double t = 0.0;

  	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    	stima_temp_ = (2.0 / 3.0) * t/(double)npart; //Temperature
	double fattore_correttivo = sqrt(temp/stima_temp_);

	for (int i=0; i<npart; ++i){

    	  vx[i] *= fattore_correttivo; 
    	  vy[i] *= fattore_correttivo;  
    	  vz[i] *= fattore_correttivo;

	  //cout<<"vx[i] = "<<vx[i]<<"   "<<"vy[i] = "<<vy[i]<<"   "<<"vz[i] = "<<vz[i]<<endl;

    	  x[i] = Pbc( xnew_[i] - vx[i]*delta );
    	  y[i] = Pbc( ynew_[i] - vy[i]*delta );
    	  z[i] = Pbc( znew_[i] - vz[i]*delta );

	  xold[i] = x[i];
    	  yold[i] = y[i];
    	  zold[i] = z[i];

    	  x[i] = xnew_[i];
    	  y[i] = ynew_[i];
    	  z[i] = znew_[i];

	  //cout<<"xold[i] = "<<xold[i]<<"   "<<"yold[i] = "<<yold[i]<<"   "<<"zold[i] = "<<zold[i]<<endl;	
	  //cout<<"x[i] = "<<x[i]<<"   "<<"y[i] = "<<y[i]<<"   "<<"z[i] = "<<z[i]<<endl;	
	  //cout<<endl;     
 	}

  }else if(oldOrNot == 0){


//Prepare initial velocities
   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   	double sumv[3] = {0.0, 0.0, 0.0};
   	for (int i=0; i<npart; ++i){
     	  vx[i] = rand()/double(RAND_MAX) - 0.5;
     	  vy[i] = rand()/double(RAND_MAX) - 0.5;
     	  vz[i] = rand()/double(RAND_MAX) - 0.5;

     	  sumv[0] += vx[i];
     	  sumv[1] += vy[i];
     	  sumv[2] += vz[i];
   	}
   	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   	  double sumv2 = 0.0, fs;
   	for (int i=0; i<npart; ++i){
     	  vx[i] = vx[i] - sumv[0];
     	  vy[i] = vy[i] - sumv[1];
     	  vz[i] = vz[i] - sumv[2];

     	  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   	}
   	sumv2 /= (double)npart;

   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   	for (int i=0; i<npart; ++i){
     	  vx[i] *= fs;
     	  vy[i] *= fs;
     	  vz[i] *= fs;

     	  xold[i] = Pbc(x[i] - vx[i] * delta);
     	  yold[i] = Pbc(y[i] - vy[i] * delta);
     	  zold[i] = Pbc(z[i] - vz[i] * delta);
   	  }
   }
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    //cout<<"xold[i] = "<<x[i]<<"   "<<"yold[i] = "<<ynew<<"   "<<"zold[i] = "<<znew<<endl;

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update g(r)
     bin_size = box/(2*100);
     bin_number = binSearcher(bin_size, nbin, dr);
     //cout << "dist = " << dr << " bin = " << bin_number << endl;
     if(bin_number != -1){
         gr[bin_number] = gr[bin_number] + 2;
     }
     

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i){ t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    	//cout<<"vx[i] = "<<vx[i]<<"   "<<"vy[i] = "<<vy[i]<<"   "<<"vz[i] = "<<vz[i]<<endl;   
  }

    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    //cout<<"t = "<<t<<endl;
    //cout<<"temp = "<<stima_temp<<endl;

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    sum_epot += stima_pot;
    sum_ekin += stima_kin;
    sum_etot += stima_etot;
    sum_temp += stima_temp;

    nmeasure++;    

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  cout << "Print the old configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void Blocking(void){

  ofstream Epot, Ekin, Etot, Temp, Gr;

  //aggiorno al blocco attuale
  ave_ekin += sum_ekin/nInBlocks;
  ave_epot += sum_epot/nInBlocks;
  ave_etot += sum_etot/nInBlocks;
  ave_temp += sum_temp/nInBlocks;

  epot_sqr += pow(sum_epot/nInBlocks,2);
  ekin_sqr += pow(sum_ekin/nInBlocks,2);
  etot_sqr += pow(sum_etot/nInBlocks,2);
  temp_sqr += pow(sum_temp/nInBlocks,2);

  double r, normalization;

  for(int i=0; i<nbin; i++){

    r = i*bin_size;
    normalization = rho*(double)npart*(4.*M_PI/3.)*(pow(r+bin_size,3)-pow(r,3));

    gr_mean[i] += (double)gr[i]/((double)nInBlocks*normalization);
    gr_sqr[i] += ((double)gr[i]/((double)nInBlocks*normalization))*((double)gr[i]/((double)nInBlocks*normalization));
	
    
    Gr.open("output.gofr.0",ios::app);
    Gr << setw(20) << gr_mean[i]/(blockNumber+1);
    Gr.close();
    Gr.open("output.gofr_err.0",ios::app);

    if(blockNumber == 0){
        Gr << setw(20) << 0;
    }else{
        Gr << setw(20) << sqrt((1./blockNumber)*((1./(blockNumber+1))*(gr_sqr[i])-pow(gr_mean[i]/(blockNumber+1),2)));
    }

    Gr.close();
  }

  Gr.open("output.gofr.0",ios::app);
  Gr << endl;
  Gr.close();

  Gr.open("output.gofr_err.0",ios::app);
  Gr << endl;
  Gr.close();
  //cout << nmeasure << endl;
  //cout << etot_sqr << " " << ave_etot << endl;
  inc_epot = sqrt((1./blockNumber)*((1./(blockNumber+1))*(epot_sqr)-pow(ave_epot/(blockNumber+1),2)));
  inc_ekin = sqrt((1./blockNumber)*((1./(blockNumber+1))*(ekin_sqr)-pow(ave_ekin/(blockNumber+1),2)));
  inc_temp = sqrt((1./blockNumber)*((1./(blockNumber+1))*(temp_sqr)-pow(ave_temp/(blockNumber+1),2)));
  inc_etot = sqrt((1./blockNumber)*((1./(blockNumber+1))*(etot_sqr)-pow(ave_etot/(blockNumber+1),2)));

  //cout << ave_ekin << "  " << (1./nmeasure)*(ekin_sqr)- pow(ave_ekin,2)<< endl;

  sum_ekin = 0;
  sum_epot = 0;
  sum_etot = 0;
  sum_temp = 0;

  Epot.open("epot_err.dat", ios::app);
  Ekin.open("ekin_err.dat", ios::app);
  Etot.open("etot_err.dat", ios::app);
  Temp.open("temp_err.dat", ios::app);

  if(blockNumber == 0){

    Epot << ave_epot/(blockNumber+1) << " " << 0 << endl;
    Ekin << ave_ekin/(blockNumber+1) << " " << 0 << endl;
    Etot << ave_etot/(blockNumber+1) << " " << 0 << endl;
    Temp << ave_temp/(blockNumber+1) << " " << 0 << endl;

  }else{

    Epot << ave_epot/(blockNumber+1) << " " << inc_epot << endl;
    Ekin << ave_ekin/(blockNumber+1) << " " << inc_ekin << endl;
    Etot << ave_etot/(blockNumber+1) << " " << inc_etot << endl;
    Temp << ave_temp/(blockNumber+1) << " " << inc_temp << endl;
  }

  Epot.close();
  Ekin.close();
  Etot.close();
  Temp.close();

  blockNumber++;

  for(int i=0; i<nbin; i++){
    gr[i] = 0;
  }
}

int binSearcher(double bin_size, int nbin, double distance){

	int _K = nbin;
	double domain_selected[2]={0,bin_size*nbin};
	double interval_lenght = (domain_selected[1]-domain_selected[0])/_K;

	int range_int[2] = {1, _K}; //intervalli numerati da 1 a _K; man mano che vado avanti con la ricerca riduco questo range.
	int INT_NUMBER = 0; //numero dell'intervallo

	while (INT_NUMBER == 0){
		if (distance > domain_selected[1]){
		INT_NUMBER = -1;
		}else{
		
			if (_K % 2 == 0){									//PARI
				if (distance < (domain_selected[1]+domain_selected[0])/2.){		//se vale questa...
					domain_selected[1] = (domain_selected[1]+domain_selected[0])/2.;	//allora devo cercare nella metà di sinistra degli intervalli!
					range_int[1] = range_int[1] - _K/2;
				} else if (distance == (domain_selected[1]+domain_selected[0])/2.){
					INT_NUMBER = range_int[1] - _K/2 + 1; 					//se il numero è sul bordo appartiene all'intervallo a destra.
				} else {
					domain_selected[0] = (domain_selected[1]+domain_selected[0])/2.;
					range_int[0] = range_int[1] - _K/2 + 1;
				}
				_K = _K/2;
			} else {										//DISPARI
				if (distance < (domain_selected[1]+domain_selected[0])/2.){
					domain_selected[1] = (domain_selected[1]+domain_selected[0])/2.-interval_lenght/2.;
					range_int[1] = range_int[1] - _K/2 - 1;				//ricorda che k_/2 se k_ è dispari...	
					if (domain_selected[1] <= distance){
						INT_NUMBER = range_int[1] + 1;
					}
				} else {
					domain_selected[0] = (domain_selected[1]+domain_selected[0])/2.+interval_lenght/2.;
					range_int[0] = range_int[1] - _K/2 + 1;	
					if (distance < domain_selected[0]){
						INT_NUMBER = range_int[0] - 1;
					}
				}
				_K = _K/2;	//ricorda che k_/2 se k_ è dispari...
			}
		}
		if (range_int[0]==range_int[1] && INT_NUMBER==0){
			INT_NUMBER = range_int[0];
		}

	}
	return INT_NUMBER-1;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
