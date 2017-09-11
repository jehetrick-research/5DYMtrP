////////////////////////////////////
// Metropolis code for SU(3)
//
// This stage starts from su3dev.cpp, and adds ploop_less_slice() and tr_A(P)
// to the action.
//
// That's the plan...


#include "fermiqcd.h"
#include "trPheaders.h"

// code
#include "readmilcascii2.cpp"
#include "ploop3.cpp"


// counter for drand48() calls
//unsigned int rcount;


// This version of load_staples multiplies the staple by its appropriate
// beta: beta4 (s-s staples) or beta5 (s-t staples)

double load_staples(gauge_field &U, mdp_matrix_field &S, int dir, int parity,
		    float beta4, float beta5) {
   // Returns the average plaquette
  site x(U.lattice());
  int mu;
  mdp_matrix staple(U.nc, U.nc);
  double staple_sum=0;
  float b;
  int mu5;

  mu5 = U.ndim-1;

  forallsitesofparity(x, parity) {
     staple = 0;
     for(mu=0; mu<U.ndim; mu++) {
	if(mu != dir) {
	   b = (mu+dir <= mu5) ? beta4 : beta5;  
	   staple = staple + b*(
	   // upper staple contribution
	   U(x+dir, mu) * hermitian(U(x+mu, dir)) * hermitian(U(x, mu)) +
	   // lower staple contribution
	   hermitian(U(x+dir-mu, mu)) * hermitian(U(x-mu, dir)) * U(x-mu, mu) );
        }
     }
     S(x) = staple;
     //     staple_sum += real(trace(U(x,dir) * staple));
  }
  S.update();

  // dump staples
  /**
  forallsitesofparity(x, parity) {
     cout<<"STAPLE "<< x <<" "<< parity <<" "<< dir << endl;
     cout<< hermitian(S(x)) <<endl;
  }
  **/


  // Each of the 3 plaquettes connected to DIR are counted 2x -> 6 duplicity
  // Parity means sum on x =  nvol_gl/2. 
  // Re(Tr S) = 3 for Identity.
  // So normalization factor is:
  // 1/6 * 1/(nvol_gl/2) * 1/3 = 1/nvol_gl/9
  //  return staple_sum/U.lattice().nvol_gl/9;

  // don't compute the ave. plaquette from the staples.
  return 1.0;
}



mdp_matrix rand_su2subgrp(mdp_prng &r, int gen, float scale) {
   mdp_matrix R(3,3);
   float theta, c, s;
   float rr;

   do { theta = scale * (r.plain() - 0.5); }
   while( fabs(theta) < 0.25*scale );
   //   do { rr = drand48();  
   //        theta = scale * (rr - 0.5);
   //      rcount++; 
   //      cout <<"theta: rcount "<< rcount <<" "<< rr <<endl; }
   //   do { theta = scale * (drand48() - 0.5); rcount++; }
   //   while( fabs(theta) < 0.25*scale );

   c = cos(theta); s = sin(theta);

   R = mdp_identity(3);

   switch(gen){
   case 1:
      R(0,0) = R(1,1) = c;
      R(0,1) = R(1,0) = s*I;
      break;
   case 2:
      R(0,0) = R(1,1) = c;
      R(0,1) = s;
      R(1,0) = -s;
      break;
   case 3:
      R(0,0) = c + (s*I);
      R(1,1) = c + (-s*I);
      break;
   case 4:
      R(0,0) = R(2,2) = c;
      R(0,2) = R(2,0) = s*I;
      break;
   case 5:
      R(0,0) = R(2,2) = c;
      R(0,2) = s;
      R(2,0) = -s;
      break;
   case 6:
      R(1,1) = R(2,2) = c;
      R(1,2) = R(2,1) = s*I;
      break;
   case 7:
      R(1,1) = R(2,2) = c;
      R(1,2) = s;
      R(2,1) = -s;
      break;
   case 8:
      R(0,0) = R(1,1) = cos(theta/sqrt(3.)) + (I*sin(theta/sqrt(3.)));
      R(2,2) = cos(-2*theta/sqrt(3.)) + (I*sin(-2*theta/sqrt(3.)));
      break;
   }

   return R;
}



mdp_matrix make_change(gauge_field &U, mdp_site &x, float scale) {
   mdp_matrix R(3,3);
   float theta, c, s, r;
   //   float mpd_prng *r;
   //   float randno;
   int gen;
   

   //   r = U.lattice().random(x);

   
   do { theta = scale * (U.lattice().random(x).plain() - 0.5); }
   while( fabs(theta) < 0.25*scale );

   c = cos(theta); s = sin(theta);
   R = mdp_identity(3);

   gen = (int)(8.0*(U.lattice().random(x).plain())) + 1;

   switch(gen){
   case 1:
      R(0,0) = R(1,1) = c;
      R(0,1) = R(1,0) = s*I;
      break;
   case 2:
      R(0,0) = R(1,1) = c;
      R(0,1) = s;
      R(1,0) = -s;
      break;
   case 3:
      R(0,0) = c + (s*I);
      R(1,1) = c + (-s*I);
      break;
   case 4:
      R(0,0) = R(2,2) = c;
      R(0,2) = R(2,0) = s*I;
      break;
   case 5:
      R(0,0) = R(2,2) = c;
      R(0,2) = s;
      R(2,0) = -s;
      break;
   case 6:
      R(1,1) = R(2,2) = c;
      R(1,2) = R(2,1) = s*I;
      break;
   case 7:
      R(1,1) = R(2,2) = c;
      R(1,2) = s;
      R(2,1) = -s;
      break;
   case 8:
      R(0,0) = R(1,1) = cos(theta/sqrt(3.)) + (I*sin(theta/sqrt(3.)));
      R(2,2) = cos(-2*theta/sqrt(3.)) + (I*sin(-2*theta/sqrt(3.)));
      break;
   }

   return R;
}



int metropolis_trP(gauge_field &U, mdp_matrix_field &S, mdp_matrix_field &Pls,
		   mdp_matrix_field &T1, mdp_matrix_field &T2, 
		   float beta4, float beta5, float h) {

  site x(U.lattice());
  int mu, t, nt;
  float scale=1.1;  int gen;
  mdp_matrix deltaSU3(U.nc,U.nc);
  mdp_matrix newU(U.nc,U.nc);
  mdp_matrix tmpmat(U.nc, U.nc);
  float oldaction, newaction;
  int accept=0, reject=0;
  float r;
  double staple_plaq;
  mdp_complex zP;
  mdp_complex trPf;
  double trPA;

  zP = mdp_complex(0,0);

  for(int parity=1; parity>=0; parity--)
     for(mu=0; mu<U.ndim-1; mu++) {  // do space-like updates mu=0,1,2,...ndim-2
	staple_plaq = load_staples(U, S, mu, parity, beta4, beta5);
     //     cout << "partial plaq from staple: " << mu << " " << staple_plaq << endl;

     forallsitesofparity(x, parity) {

	deltaSU3 = make_change(U, x, scale);
	newU = deltaSU3 * U(x,mu);
	
	// note: beta4 and beta5 are already in the staples
	oldaction=0.333333*real(trace( U(x,mu) * S(x) ));
	newaction=0.333333*real(trace( newU * S(x) ));
	
	/**
	cout << " x: " << x << " mu: " << mu << endl;
	cout << "gen: " << gen << endl;
	cout << "deltaSU3: " << endl << deltaSU3 << endl;
	cout << "det(deltaSU3) = " << det(deltaSU3) << endl;
	cout << "S: " << endl << S(x) << endl;
	cout << "U: " << endl << U(x,mu) << endl;
	cout << "newU: " endl << newU << endl;
	cout << "oldaction: " << oldaction << " newaction: " << newaction << endl;
	//	cout << "r: " << U.lattice().random(x).plain() << endl;
	**/

	if( newaction > oldaction ){
	   U(x,mu)=newU;
	   //	   cout << "new action higher: accepted" << endl;
	   accept++;
	}
	else {
	   r = U.lattice().random(x).plain();
	   //	   cout << "exp(new-old): " << exp(newaction - oldaction) << endl;
	   if( r < exp(newaction - oldaction) ) {
	      U(x,mu) = newU;
	      accept++;
	      //	      cout << "r < exp: accepted" << endl;
	   }
	   //	} else{
	   //	   reject++;
	}
     }
  }
  U.update();

  // "time-like" update
  mu = U.ndim-1;
  nt = U.lattice().size(mu);

  for(int parity=1; parity>=0; parity--) {
     staple_plaq = load_staples(U, S, mu, parity, beta4, beta5);

     int opp_parity;
     opp_parity = parity==0 ? 1 : 1;

     for(t=0; t<nt; t++) {
	zP += ploop_less_slice_t(U, Pls, mu, t, T1, T2, parity);

	forallsitesofparity(x, parity) if(x(mu)==t) {

	   deltaSU3 = make_change(U, x, scale);
	   newU = deltaSU3 * U(x,mu);
	   
	   // Old action: [beta4 and beta5 are in staples]
	   oldaction=0.333333*real(trace( U(x,mu) * S(x) ));
	   // add (-h tr_A P)
	   tmpmat = U(x,mu) * Pls(x);
	   trPf = trace(tmpmat);
	   trPA = real(trPf)*real(trPf) + imag(trPf)*imag(trPf) - 1.0;
	   oldaction -= h*trPA;
	   
	   // New action
	   newaction=0.333333*real(trace( newU * S(x) ));
	   // add (-h tr_A P)
	   tmpmat = newU * Pls(x);
	   trPf = trace(tmpmat);
	   trPA = real(trPf)*real(trPf) + imag(trPf)*imag(trPf) - 1.0;
	   newaction -= h*trPA;
	   	   
	   /**
	      cout << " x: " << x << " mu: " << mu << endl;
	      cout << "gen: " << gen << endl;
	      cout << "deltaSU3: " << endl << deltaSU3 << endl;
	      cout << "det(deltaSU3) = " << det(deltaSU3) << endl;
	      cout << "S: " << endl << S(x) << endl;
	      cout << "U: " << endl << U(x,mu) << endl;
	      cout << "newU: " endl << newU << endl;
	      cout << "oldaction: " << oldaction << " newaction: " << newaction << endl;
	      //	cout << "r: " << U.lattice().random(x).plain() << endl;
	   **/
	   
	   if( newaction > oldaction ){
	      U(x,mu)=newU;
	      //	   cout << "new action higher: accepted" << endl;
	      accept++;
	   }
	   else {
	      r = U.lattice().random(x).plain();
	      //	   cout << "exp(new-old): " << exp(newaction - oldaction) << endl;
	      if( r < exp(newaction - oldaction) ) {
		 U(x,mu) = newU;
		 accept++;
		 //	      cout << "r < exp: accepted" << endl;
	      }
	      //	} else{
	      //	   reject++;
	   }
	}
	U.update();
     }
  }
     
  // Report P-loop average from metropolis update [Doesn't work]
  //  zP /= 8;
  //  if(ME==0) cout <<"METPLOOP "<< real(zP) <<" "<< imag(zP) <<endl;
     
  return accept;
}




/////////////////////////////////////////////////////////
// Main
/////////////////////////////////////////////////////////

int main(int argc, char** argv) {

  int verbose=false;
  int ndim=4,nc=3;
  int L[10]={4,4,4,4,4,4,4,4,4,4};
  mdp_real beta=6.0;
  int steps=0;
  char input[1024]="";
  char output[1024]="";
  int mode=0;
  int warms=0;
  int trajecs=1;
  int meas=1;
  register int i, j;
  double plaq;
  mdp_complex zP4, zP5;
  float Pre, Pim;
  float h=0;
  double staple_plaq=0;
  int mu, parity;
  long seed=-666;
  mdp_complex Pslices;
  float beta4=6.0, beta5=6.0, gamma=1.0;


  // //////////////////////////////
  // Parsing command line arguments
  // //////////////////////////////
  for(int i=1; i<argc; i++) {
    if(strncmp(argv[i],"-verbose",8)==0)      verbose=true;
    else if(strncmp(argv[i],"-hot",4)==0)     mode=1;
    else if(strncmp(argv[i],"-cold",5)==0)    mode=0;
    else if(strncmp(argv[i],"-seed",5)==0)    sscanf(argv[i+1],"%ld",&seed);
    else if(strncmp(argv[i],"-input",6)==0)   {mode=2; sscanf(argv[i+1],"%s",input); }
    else if(strncmp(argv[i],"-readmilcascii",14)==0) {mode=3; sscanf(argv[i+1],"%s",input); }  
    else if(strncmp(argv[i],"-output",7)==0)  sscanf(argv[i+1],"%s",output);
    else if(strncmp(argv[i],"-warms",6)==0)   sscanf(argv[i+1],"%i",&warms);  
    else if(strncmp(argv[i],"-trajecs",8)==0) sscanf(argv[i+1],"%i",&trajecs);  
    else if(strncmp(argv[i],"-meas",5)==0)    sscanf(argv[i+1],"%i",&meas);  
    else if(strncmp(argv[i],"-nc",3)==0)      sscanf(argv[i+1],"%i",&nc);
    else if(strncmp(argv[i],"-steps",6)==0)   sscanf(argv[i+1],"%i",&steps);
    else if(strncmp(argv[i],"-beta4",6)==0)    sscanf(argv[i+1],"%f",&beta4);    
    else if(strncmp(argv[i],"-beta5",6)==0)    sscanf(argv[i+1],"%f",&beta5);    
    else if(strncmp(argv[i],"-gamma",6)==0)    sscanf(argv[i+1],"%f",&gamma);    
    else if(strncmp(argv[i],"-beta",5)==0)    sscanf(argv[i+1],"%f",&beta);    
    else if(strncmp(argv[i],"-H",2)==0)       sscanf(argv[i+1],"%f",&h);    
    else if(strncmp(argv[i],"-L",2)==0)       ndim=sscanf(argv[i+1],"%ix%ix%ix%ix%ix%ix%ix%ix%ix%i",
							  L,L+1,L+2,L+3,L+4,L+5,L+6,L+7,L+8,L+9);
    else if(strncmp(argv[i],"-help",5)==0) {
       mdp << "SU(3) Metropolis updates\n";
       mdp << "usage:\n";
       mdp << "-hot / -cold start [default -cold]\n";
       mdp << "-L NxNxNxN\n";
       mdp << "-output filename\n";
       mdp << "-input filename\n";
       //       mdp << "-nc Ncolors\n";
       mdp << "-beta \n";
       mdp << "-beta4 \n";
       mdp << "-beta5 \n";
       mdp << "-gamma \n";
       mdp << "-H \n";
       mdp << "-warms \n";
       mdp << "-trajecs \n";
       mdp << "-seed \n";
       mdp << "-meas  [trajecs between measurements]\n";
       //       mdp << "works for any nc and up to 10 dimensions.\n";
       exit(1);
    }
  }

  // open communications
  mdp.open_wormholes(argc,argv);


  // Setup Lattice and initialize /////////////////////////
  if(mode==2) {
    mdp_field_file_header header;
    if(is_file(input)) header=get_info(input);
    else error("Unable to access input gauge configuration\n");
    ndim=header.ndim;
    nc=(int) sqrt((double) header.bytes_per_site/(ndim*sizeof(mdp_complex)));
    for(i=0; i<ndim; i++) L[i]=header.box[i];
  }
  
  if(!verbose)  mdp.print=false;  // eventualy print off
  mdp_lattice   lattice(ndim,L); // declare lattice
  mdp_site      x(lattice);      // declare site variable
  gauge_field	U(lattice,nc);   // declare SU(3) field
  coefficients	gauge;		 // declare coefficients


  if(seed == -666) {
     seed = (long)(time(NULL)-1495056447);
     lattice.initialize_random(seed);
     if(ME==0) { cout << "Initialized lattice random number generator with seed = "
		      << seed << endl; }
  } else {
     lattice.initialize_random(seed);
     if(ME==0) { cout << "Initialized lattice random number generator with seed = "
		      << seed << endl; }
  }

  // in case of drand48()
  // srand48(seed);
  //    rcount=0;
  //    cout <<"r1 = "<< drand48() << endl;

  // Anisotropy parameters
  // either enter parameters as:
  //    -beta 5.8 -gamma 1.5
  // or
  //    -beta4 5.8 -beta5 6.3
  //
  // any mixure of these will cause chaos

  if(gamma!=1) {
     beta4 = beta/gamma;
     beta5 = beta*gamma;
  }
  if(beta4!=beta5) { 
     gamma = sqrt(beta5/beta4);
     beta = sqrt(beta5*beta4);
  }


  
  // //////////////////////////////
  // Output parameters
  // //////////////////////////////
  mdp << "=============================================\n";
  mdp << "Number of colors = " << nc << '\n';
  mdp << "Lattice size = " << L[0]; 
  for(i=1; i<ndim; i++) mdp << "x" << L[i]; mdp << '\n';
  mdp << "Performing Metropolis updates with h*Tr_A(P) term" endl;
  if(ME==0) {
     cout <<"L ";
     for(i=0; i<ndim; i++) cout <<" "<< L[i]; cout << endl;
     cout <<"Nc "<< nc <<endl;
     cout <<"startmode "<< mode <<endl;
     cout <<"beta "<< beta <<endl;
     cout <<"beta4 "<< beta4 <<endl;
     cout <<"beta5 "<< beta5 <<endl;
     cout <<"gamma "<< gamma <<endl;
     cout <<"H "<< h <<endl;
     cout <<"warms "<< warms <<endl;
     cout <<"trajecs "<< trajecs <<endl;
     cout <<"meas "<< meas <<endl;
     cout <<"seed "<< seed <<endl;
     cout <<"input/output ["<< input <<"/"<< output <<"]"<<endl;
  }
  switch(mode) {
    case 0: mdp << "Creating a cold configuration\n"; break;
    case 1: mdp << "Creating a hot configuration\n"; break;
    case 2: mdp << "Reading MDP input file = " << input << '\n'; break;
    case 3: mdp << "Reading MILC ascii input file = " << input << '\n'; break;
  }
  mdp << "Saving output in file = " << output << '\n';
  mdp << "=============================================\n";



  
  // Field to hold the 6 Plaquettes /////
  mdp_nmatrix_field  P(lattice,6, U.nc,U.nc);
  forallsites(x) {
     for(i=0; i<6; i++) {
	P(x,i) = mdp_identity(U.nc);
     }
  }
  // Field to hold the Staple sum /////
  mdp_matrix_field  S(lattice,U.nc,U.nc);
  //  forallsites(x)  S(x) = mdp_identity(U.nc);

  // for Ploop_less_slice
  mdp_matrix_field  Pls(lattice,U.nc,U.nc);


  // tmp matrix fields
  mdp_matrix_field  T1(lattice,U.nc,U.nc);
  mdp_matrix_field  T2(lattice,U.nc,U.nc);
  mdp_matrix_field  T3(lattice,U.nc,U.nc);



  gauge["beta"]=beta;		 // set beta
  gauge["beta4"]=beta4;		 // set beta4
  gauge["beta5"]=beta5;		 // set beta5
  gauge["gamma"]=gamma;		 // set gamma
  gauge["h"]=h;		         // set h


  // initial gauge field
  switch(mode) {  
    case 0: set_cold(U); if(ME==0)cout<<"Initialized lattice: COLD\n"; break;
    case 1: set_hot(U); if(ME==0)cout<<"Initialized lattice: HOT\n"; break;       
    case 2: U.load(input);  break;
       // requires readmilcascii2.cpp at top
    case 3: read_ascii_lat(input, U); break;
  }



  // DEBUGGING ///////////////////////////////
  ////////////////////////////////////////////

  // Special Debugging lattices?  
  // diagonal integer time-like links for Ploop testing
  /**
  cout <<"Setting diagonal integer time-like links\n";
  mdp_matrix id3 = mdp_identity(3);
  forallsites(x) {
     U(x,3) = (1 + x(3)) * id3;
  }
  forallsites(x) {
     cout <<"U("<< x <<")\n";
     cout << U(x,3);
  }
  **/

  // Dump Ploop-less-slice
  /**
  int t, nt;
  mu = U.ndim-1;
  nt = U.lattice().size(mu);

  for(int parity=1; parity>=0; parity--) {
     int opp_parity;
     opp_parity = parity==0 ? 1 : 1;

     for(t=0; t<nt; t++) {
	zP4 = ploop_less_slice_t(U, Pls, mu, t, T1, T2, parity);

	//	forallsitesofparity(x, parity) if(x(mu)==t) {
	//	   cout <<"TEST:PLS x = "<< x <<" "<< real(Pls(x)[0,0]) <<endl;
	//      }
     }
  }
  **/

  /**
    cout << "size: "<< lattice.size() <<endl;
    cout << "U.ndim: "<< U.ndim << endl;
    cout << "nt = "<< lattice.size(U.ndim-1) <<endl;
    cout << "--" <<endl;
  **/

  //  cout << "Dump Staples" << endl;
  /**
  double splaqsum=0;
  for(parity=1; parity>=0; parity--)
     for(mu=0; mu<4; mu++) {
	cout << "p=" << parity <<" mu="<< mu <<endl;
	splaqsum += staple_plaq = load_staples(U, S, mu, parity);
	cout <<"staple_plaq: "<< staple_plaq << endl;
     }
  cout << "staple_plaq mean: "<< splaqsum/6 << endl;
  **/

  //  cout << "nvol_gl = " << U.lattice().nvol_gl << endl;

  // Print out topology and Node assignment
  /**
  forallsites(x) {
     cout <<"x= "<< x(0) <<" y= "<< x(1) <<" z= "<< x(2) <<" t= "<< x(3) <<" ME= "<< ME <<endl;
  }
  **/






  ///////////////////////////////////////////////////////////
  // Initial MEASUREMENTS
  ///////////////////////////////////////////////////////////

  plaq = average_plaquette(U);
  zP4 = ave_polyakov_loop(U, 3, T1, T2, T3);
  zP5 = ave_polyakov_loop(U, 4, T1, T2, T3);
  if(ME==0) {
     cout << "INIT "<< beta4 <<" "<< beta5 <<" "<< h <<" "<< plaq <<" "
	  << zP4.real() <<" "<< zP4.imag() <<" "<< abs(zP4) <<" "
	  << zP5.real() <<" "<< zP5.imag() <<" "<< abs(zP5) <<endl; 
  }


  ///////////////////////////////////////////////////////////
  // Thermalize /////////////////////////////////////////////
  ///////////////////////////////////////////////////////////

  for(i=0; i<warms; i++) {
     metropolis_trP(U, S, Pls, T1, T2, beta4, beta5, h);
  }



  ///////////////////////////////////////////////////////////
  // Do the measurement trajectories ////////////////////////
  ///////////////////////////////////////////////////////////

  for(i=0; i<trajecs; i++) {
     metropolis_trP(U, S, Pls, T1, T2, beta4, beta5, h);

     // Measurements //
     if(i%meas==0) {
	//	cout << "Trajectory " << i << endl;
	plaq = average_plaquette(U);
	zP4 = ave_polyakov_loop(U, 3, T1, T2, T3);
	zP5 = ave_polyakov_loop(U, 4, T1, T2, T3);
	if(ME==0) { 
	   cout << "MEAS "<< beta4 <<" "<< beta5 <<" "<< h <<" "<< plaq <<" "
		<< zP4.real() <<" "<< zP4.imag() <<" "<< abs(zP4) <<" "
		<< zP5.real() <<" "<< zP5.imag() <<" "<< abs(zP5) <<endl;
	}

     }
  }



  // Output lattice //////////////////////////////////////
  if(string(output)!="") U.save(output);      // save file


  // Shutdown gracefully /////////////////////////////////
  mdp.close_wormholes();
  return 0;
}



		  
