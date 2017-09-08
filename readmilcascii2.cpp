// Utility to read in a MILC save_ascii file

int read_ascii_lat(char *fname, gauge_field &U) {
   char line[100];
   FILE *fp;
   int x,y,z,t,i,a,b,dir, Nt,Nx,Ny,Nz;
   mdp_site siteX(U.lattice());
   mdp_matrix M(U.nc,U.nc);
   int version, mnx,mny,mnz,mnt;
   float Ur, Ui;
   long vol;

   Nt = U.lattice().size(0);
   Nx = U.lattice().size(1);
   Ny = U.lattice().size(2);
   Nz = U.lattice().size(3);

   fp = fopen(fname, "r");

   cout << "reading ASCII MILC lattice from " << fname << endl;
   // process header
   fscanf(fp,"%d\n", &version); cout <<"milc ascii file: version "<< version <<endl;
   if( fgets (line, 100 , fp)!=NULL ) {
      cout <<"milc ascii file: timestamp "<< line <<endl;
   }
   //   fscanf(fp,"%s\n", line); cout <<"milc ascii file: timestamp "<< line <<endl;
   fscanf(fp,"%d\t%d\t%d\t%d\n", &mnx,&mny,&mnz,&mnt); 
   cout <<"milc ascii file: size = "<< mnx <<" "<< mny <<" "<< mnz <<" "<< mnt <<endl;

   if( (Nx!=mnx)||(Ny!=mny)||(Nz!=mnz)||(Nt!=mnt) ) {
      cout <<"Size mismatch: "<<endl;
      cout <<"milc ascii file: size = "<< mnx <<" "<< mny <<" "<< mnz <<" "<< mnt <<endl;
      cout <<"U.lattice: size = "<< Nx <<" "<< Ny <<" "<< Nz <<" "<< Nt <<endl;
      exit(0);
   }

   vol = Nx*Ny*Nz*Nt;

   //   ### Define Nt, Nx, etc.

   for(i=0; i<4*vol; i++) {
      // get x= y= z= t= dir=line
      fscanf(fp,"x=%d y=%d z=%d t=%d dir=%d\n", &x,&y,&z,&t,&dir);
      siteX.set(x,y,z,t);
      //      cout <<  "set siteX: " <<t<<" "<<x<<" "<<y<<" "<<z<<" dir="<< dir <<endl;
      //read a link
      for(a=0; a<3; a++) for(b=0; b<3; b++) {
	    fscanf(fp,"%e %e\n", &Ur, &Ui);
	    M(a,b) = mdp_complex(Ur, Ui);
      }
      U(siteX, dir) = M;
      //      cout << "U: siteX=" << siteX <<"\n"<< U(siteX, dir);
   }

   fclose(fp);
}
		  
