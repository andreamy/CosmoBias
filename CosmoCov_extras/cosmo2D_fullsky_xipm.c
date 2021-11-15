double w_tomo_fullsky(int nt, int ni, int nj);
// galaxy clustering tomography 2PCF galaxies in bins ni, nj, computed as sum P_l(cos(like.theta[nt])*C(l,ni,nj)
double w_gamma_t_fullsky(int nt,int ni, int nj); //G-G lensing, lens bin ni, source bin nj, including IA contamination if like.IA = 3
double xi_pm_fullsky_FINEBINNING(int pm, double *theta, int nt, int ni, int nj); //shear tomography correlation functions, including IA contamination if like.IA = 3
double xi_pm_fullsky(int pm, int nt, int ni, int nj);
double xi_pm_tomo(int pm, double theta, int ni, int nj); //shear tomography correlation functions, including IA contamination if like.IA = 3
double xi_pm_rebin(int pm, double thetamin_i, double thetamax_i, int ni,int nj);//xi_pm averaged over large bins

/**************** these routines are only used internally ************/
typedef double (*C_tomo_pointer)(double l, int n1, int n2);
void xipm_via_hankel(double **xi, double *logthetamin, double *logthetamax,  C_tomo_pointer C_tomo,int ni, int nj);
void hankel_kernel_FT(double x, fftw_complex *res, double *arg, int argc);
/*************** look-up tables for angular correlation functions ***************/
/******************** all angles in radian! *******************/

double w_tomo_fullsky(int nt, int ni, int nj){
	static int LMAX = 100000;
	static int NTHETA = 0;
	static double ** Pl =0;
	static double *Cl =0;
	static double *w_vec =0;
	static cosmopara C;
	static nuisancepara N;
	static galpara G;
	int i,l,nz;
	if (like.Ntheta ==0){
		printf("cosmo2D_real.c:w_tomo_exact: like.Ntheta not initialized\nEXIT\n"); exit(1);
	}
	if (ni != nj){
		printf("cosmo2D_real.c:w_tomo_exact: ni != nj tomography not supported\nEXIT\n"); exit(1);    
	}
	if (Pl ==0){
		Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
		Cl = create_double_vector(0,LMAX-1);
		w_vec = create_double_vector(0,tomo.clustering_Nbin*like.Ntheta-1);
		NTHETA = like.Ntheta;
		double *xmin, *xmax, *Pmin, *Pmax;
		xmin= create_double_vector(0, like.Ntheta-1);
		xmax= create_double_vector(0, like.Ntheta-1);
		double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
		Pmin= create_double_vector(0, LMAX+1);
		Pmax= create_double_vector(0, LMAX+1);
		for(i=0; i<like.Ntheta ; i++){
			xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
			xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
		}

		for (i = 0; i<NTHETA; i ++){
			printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
			gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
			gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
			//gsl_sf_legendre_Pl_array(LMAX, cos(like.theta[i]),Pmin);			
			for (int l = 1; l < LMAX; l ++){
		        //Pl[i][l] = (2*l+1.)/(4.*M_PI)*Pmin[l];
				//Pl[i][l] = (2*l+1.)/(4.*M_PI)*gsl_sf_legendre_Pl(l,cos(like.theta[i]));
				Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
			}
		}
		free_double_vector(xmin,0,like.Ntheta-1);
		free_double_vector(xmax,0,like.Ntheta-1);
		free_double_vector(Pmin,0,LMAX+1);
		free_double_vector(Pmax,0,LMAX+1);
		}
		if (recompute_clustering(C,G,N,ni,nj)){
			for (nz = 0; nz <tomo.clustering_Nbin; nz ++){
				for (l = 1; l < LMAX; l++){	
//					if (l < 20){Cl[l]=C_cl_RSD_nointerp(l,nz,nz);}
					if (l < 20){Cl[l]=C_cl_tomo_nointerp(l,nz,nz);}
					else Cl[l]=C_cl_tomo(1.0*l,nz,nz);
				}
				for (i = 0; i < NTHETA; i++){
					w_vec[nz*like.Ntheta+i] =0;
					for (l = 1; l < LMAX; l++){
						w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
					}
				}
			}
		update_cosmopara(&C);
		update_galpara(&G);
		update_nuisance(&N);
	}
	return w_vec[ni*like.Ntheta+nt];  
}

double w_gamma_t_fullsky(int nt, int ni, int nj){
	static int LMAX = 100000;
	static int NTHETA = 0;
	static double ** Pl =0;
	static double *Cl =0;
	static double *w_vec =0;
	static cosmopara C;
	static nuisancepara N;
	static galpara G;
	int i,l,nz;
	if (like.Ntheta ==0){
		printf("cosmo2D_fullsky.c:w_gamma_t_tomo: like.Ntheta not initialized\nEXIT\n"); exit(1);
	}
	if (Pl ==0){
		Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
		Cl = create_double_vector(0,LMAX-1);
		w_vec = create_double_vector(0,tomo.ggl_Npowerspectra*like.Ntheta-1);
		NTHETA = like.Ntheta;
		double *xmin, *xmax, *Pmin, *Pmax, *dP;
		xmin= create_double_vector(0, like.Ntheta-1);
		xmax= create_double_vector(0, like.Ntheta-1);
		double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
		for(i=0; i<like.Ntheta ; i++){
			xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
			xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
		}
		Pmin= create_double_vector(0, LMAX+1);
		Pmax= create_double_vector(0, LMAX+1);

		for (i = 0; i<NTHETA; i ++){
			printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
			gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
			gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
			for (int l = 1; l < LMAX; l ++){
				//Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1))*gsl_sf_legendre_Plm(l,2,cos(like.theta[i]));	
				Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
				*((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
				+(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
				-2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
			}
		}
		free_double_vector(xmin,0,like.Ntheta-1);
		free_double_vector(xmax,0,like.Ntheta-1);
		free_double_vector(Pmin,0,LMAX+1);
		free_double_vector(Pmax,0,LMAX+1);
	}
	if (recompute_ggl(C,G,N,ni)){
		if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: w_gamma_t_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}
		C_tomo_pointer C_gl_pointer = &C_gl_tomo;
		if (like.IA ==3 || like.IA ==4) C_gl_pointer = &C_ggl_IA_tab;

		for (nz = 0; nz <tomo.ggl_Npowerspectra; nz ++){
			for (l = 1; l < LMAX; l++){
				Cl[l]=C_ggl_IA_tab(1.0*l,ZL(nz),ZS(nz));
			}
			for (i = 0; i < NTHETA; i++){
				w_vec[nz*like.Ntheta+i] =0;
				for (l = 2; l < LMAX; l++){
					w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
				}
			}
		}
		update_cosmopara(&C);
		update_galpara(&G);
		update_nuisance(&N);
	}
	return w_vec[N_ggl(ni,nj)*like.Ntheta+nt];  
}

double xi_pm_fullsky(int pm, int nt, int ni, int nj) //shear tomography correlation functions
{
    static int LMAX = 100000;
    static double **Glplus =0;
    static double **Glminus =0;
    static double *Cl =0;
    static double *xi_vec_plus =0;
    static double *xi_vec_minus =0;
    static cosmopara C;
    static nuisancepara N;

    int i,l,nz;
    if (like.Ntheta == 0){
        printf("cosmo2D_fullsky.c:xi_pm_tomo: like.theta not initialized\nEXIT\n"); exit(1);
    }
    if (Glplus ==0){
        Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
        Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
        Cl = create_double_vector(0,LMAX-1);
        xi_vec_plus = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
        xi_vec_minus = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
        double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
        xmin= create_double_vector(0, like.Ntheta-1);
        xmax= create_double_vector(0, like.Ntheta-1);

        double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
        double *log_theta_left, *theta_center;
        log_theta_left = create_double_vector(0, like.Ntheta);
        theta_center = create_double_vector(0, like.Ntheta);

        for(i=0; i<=like.Ntheta ; i++){
            log_theta_left[i] = log(like.vtmin)+(i+0.0)*logdt;
            theta_center[i] = exp(log_theta_left[i] + logdt*0.5);
        }

        for(i=0; i<like.Ntheta ; i++){
            xmin[i]=cos(theta_center[i]);
            xmax[i]=cos(theta_center[i+1]);
        }

        Pmin= create_double_vector(0, LMAX+1);
        Pmax= create_double_vector(0, LMAX+1);
        dPmin= create_double_vector(0, LMAX+1);
        dPmax= create_double_vector(0, LMAX+1);
        for (i = 0; i<like.Ntheta; i ++){
            //These functions compute an array of Legendre polynomials P_l(x), and optionally their derivatives dP_l(x)/dx, for l = 0, \dots, lmax, |x| <= 1
            gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
            gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
            for (int l = 3; l < LMAX; l ++){
                /*double plm = gsl_sf_legendre_Plm(l,2,x);
                double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
                Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
                *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
                +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


                Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
                *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
                +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

                Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

                        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
                        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
                        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])

                        +(4-l)   * (dPmin[l]-dPmax[l])
                        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

                        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
                        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

                )/(xmin[i]-xmax[i]);

                Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

                        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
                        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
                        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])

                        +(4-l)   * (dPmin[l]-dPmax[l])
                        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

                        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
                        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

                )/(xmin[i]-xmax[i]);

            }
        }
        free_double_vector(xmin,0,like.Ntheta-1);
        free_double_vector(xmax,0,like.Ntheta-1);
        free_double_vector(Pmin,0,LMAX+1);
        free_double_vector(Pmax,0,LMAX+1);
        free_double_vector(dPmin,0,LMAX+1);
        free_double_vector(dPmax,0,LMAX+1);

    }
    if (recompute_shear(C,N)){ //True with current settings, enters in conditional
        //if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: xi_pm_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}

        //C_tomo_pointer C_pointer = &C_shear_tomo;
        //if (like.IA ==3 || like.IA ==4) {C_pointer = &C_shear_shear_IA_tab;}

        update_cosmopara(&C); update_nuisance(&N);
        for (nz = 0; nz <tomo.shear_Npowerspectra; nz ++){
            for (l = 2; l < LMAX; l++){
                Cl[l]=C_shear_shear_IA_tab(1.0*l,Z1(nz),Z2(nz));
            }
            for (i = 0; i < like.Ntheta; i++){
                xi_vec_plus[nz*like.Ntheta+i] =0;
                xi_vec_minus[nz*like.Ntheta+i] =0;
                for (l = 2; l < LMAX; l++){
                    xi_vec_plus[nz*like.Ntheta+i]+=Glplus[i][l]*Cl[l];
                    xi_vec_minus[nz*like.Ntheta+i]+=Glminus[i][l]*Cl[l];
                }
            }
        }
    }
    if (pm> 0) return xi_vec_plus[N_shear(ni,nj)*like.Ntheta + nt];
    return xi_vec_minus[N_shear(ni,nj)*like.Ntheta + nt];
}


double xi_pm_fullsky_FINEBINNING(int pm, double *theta, int nt, int ni, int nj) //shear tomography correlation functions
{
    static int LMAX = 100000;
    static double **Glplus =0;
    static double **Glminus =0;
    static double *Cl =0;
    static double *xi_vec_plus =0;
    static double *xi_vec_minus =0;
    static cosmopara C;
    static nuisancepara N;

    int i,l,nz;
    if (like.Ntheta == 0){
        printf("cosmo2D_fullsky.c:xi_pm_tomo: like.theta not initialized\nEXIT\n"); exit(1);
    }
    if (Glplus ==0){
        Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
        Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
        Cl = create_double_vector(0,LMAX-1);
        xi_vec_plus = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
        xi_vec_minus = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
        double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
        xmin= create_double_vector(0, like.Ntheta-1);
        xmax= create_double_vector(0, like.Ntheta-1);

        /*double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
        double *log_theta_left, *theta_center;
        log_theta_left = create_double_vector(0, like.Ntheta);
        theta_center = create_double_vector(0, like.Ntheta);

        for(i=0; i<=like.Ntheta ; i++){
            log_theta_left[i] = log(like.vtmin)+(i+0.0)*logdt;
            theta_center[i] = exp(log_theta_left[i] + logdt*0.5);
        }

        for(i=0; i<like.Ntheta ; i++){
            xmin[i]=cos(theta_center[i]);
            xmax[i]=cos(theta_center[i+1]);
        }*/

        for(i=0; i<like.Ntheta ; i++){
            xmin[i]=cos(theta[i]);
            xmax[i]=cos(theta[i+1]);
        }


        Pmin= create_double_vector(0, LMAX+1);
        Pmax= create_double_vector(0, LMAX+1);
        dPmin= create_double_vector(0, LMAX+1);
        dPmax= create_double_vector(0, LMAX+1);
        for (i = 0; i<like.Ntheta; i ++){
            //These functions compute an array of Legendre polynomials P_l(x), and optionally their derivatives dP_l(x)/dx, for l = 0, \dots, lmax, |x| <= 1
            gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
            gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
            for (int l = 3; l < LMAX; l ++){
                /*double plm = gsl_sf_legendre_Plm(l,2,x);
                double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
                Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
                *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
                +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


                Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
                *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
                +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

                Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

                        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
                        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
                        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])

                        +(4-l)   * (dPmin[l]-dPmax[l])
                        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

                        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
                        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

                )/(xmin[i]-xmax[i]);

                Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

                        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
                        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
                        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])

                        +(4-l)   * (dPmin[l]-dPmax[l])
                        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

                        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
                        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

                )/(xmin[i]-xmax[i]);

            }
        }
        free_double_vector(xmin,0,like.Ntheta-1);
        free_double_vector(xmax,0,like.Ntheta-1);
        free_double_vector(Pmin,0,LMAX+1);
        free_double_vector(Pmax,0,LMAX+1);
        free_double_vector(dPmin,0,LMAX+1);
        free_double_vector(dPmax,0,LMAX+1);

    }
    if (recompute_shear(C,N)){ //True with current settings, enters in conditional
        //if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: xi_pm_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}

        //C_tomo_pointer C_pointer = &C_shear_tomo;
        //if (like.IA ==3 || like.IA ==4) {C_pointer = &C_shear_shear_IA_tab;}

        update_cosmopara(&C); update_nuisance(&N);
        for (nz = 0; nz <tomo.shear_Npowerspectra; nz ++){
            for (l = 2; l < LMAX; l++){
                Cl[l]=C_shear_shear_IA_tab(1.0*l,Z1(nz),Z2(nz));
            }
            for (i = 0; i < like.Ntheta; i++){
                xi_vec_plus[nz*like.Ntheta+i] =0;
                xi_vec_minus[nz*like.Ntheta+i] =0;
                for (l = 2; l < LMAX; l++){
                    xi_vec_plus[nz*like.Ntheta+i]+=Glplus[i][l]*Cl[l];
                    xi_vec_minus[nz*like.Ntheta+i]+=Glminus[i][l]*Cl[l];
                }
            }
        }
    }
    if (pm> 0) return xi_vec_plus[N_shear(ni,nj)*like.Ntheta + nt];
    return xi_vec_minus[N_shear(ni,nj)*like.Ntheta + nt];
}



/******************************************************************************************/
/**************************** HANKEL ROUTINES, worse results though **********************/
/******************************************************************************************/

double xi_pm_rebin(int pm, double thetamin_i, double thetamax_i, int ni,int nj) {
    int Nsub_G = 10;
    int ii;
    double ti, dti;
    double xi = 0.;
    dti = (thetamax_i - thetamin_i) / (double) Nsub_G;
    for (ii = 0; ii < Nsub_G; ii++) {
        ti = 2. / 3. * (pow(thetamin_i + (ii + 1.) * dti, 3.) - pow(thetamin_i + (ii + 0.) * dti, 3.)) /
             (pow(thetamin_i + (ii + 1.) * dti, 2.) - pow(thetamin_i + (ii + 0.) * dti, 2.));
        xi += xi_pm_tomo(pm, ti, ni, nj) *
              (pow(thetamin_i + (ii + 1.) * dti, 2.) - pow(thetamin_i + (ii + 0.) * dti, 2.));
    }
    return xi / (pow(thetamax_i, 2.) - pow(thetamin_i, 2.));
}

double xi_pm_tomo(int pm, double theta, int ni, int nj) //shear tomography correlation functions
{
    static cosmopara C;
    static nuisancepara N;
    static double **table;
    static double dlogtheta, logthetamin, logthetamax;
    if (recompute_shear(C,N)){
        //if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: xi_pm_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}
        C_tomo_pointer C_pointer = &C_shear_tomo;
        if (like.IA ==3 || like.IA ==4) {C_pointer = &C_shear_shear_IA_tab;}
        update_cosmopara(&C); update_nuisance(&N);
        double **tab;
        int i,k;
        tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
        if (table==0) table   = create_double_matrix(0, 2*tomo.shear_Npowerspectra-1, 0, Ntable.N_thetaH-1);
        for (i = 0; i < tomo.shear_Npowerspectra; i++){
            xipm_via_hankel(tab, &logthetamin, &logthetamax,C_pointer, Z1(i),Z2(i));
            for (k = 0; k < Ntable.N_thetaH; k++){
                table[2*i][k] = tab[0][k];
                table[2*i+1][k] = tab[1][k];
            }
        }
        dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
        free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
    }
    return interpol(table[2*N_shear(ni,nj)+(1-pm)/2], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);
}

void xipm_via_hankel(double **xi, double *logthetamin, double *logthetamax,  C_tomo_pointer C_tomo,int ni, int nj)
{
    const double l_min = 0.0001;
    const double l_max = 5.0e6;
    static double loglmax = -123.0, loglmin, dlnl,  lnrc, arg[2];
    static int nc;

    double        l, kk, *lP, t;
    fftw_plan     plan1,plan;
    fftw_complex *f_lP,*conv;
    fftw_complex  kernel;
    int           i, count;
    lP   = fftw_malloc(Ntable.N_thetaH*sizeof(double));
    f_lP = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
    conv = fftw_malloc((Ntable.N_thetaH/2+1)*sizeof(fftw_complex));
    plan  = fftw_plan_dft_r2c_1d(Ntable.N_thetaH, lP, f_lP, FFTW_ESTIMATE);
    plan1 = fftw_plan_dft_c2r_1d(Ntable.N_thetaH, conv, lP, FFTW_ESTIMATE);
    if (loglmax==-123.0) {
        loglmax  = log(l_max);
        loglmin  = log(l_min);
        dlnl     = (loglmax-loglmin)/(1.0*Ntable.N_thetaH);
        lnrc     = 0.5*(loglmax+loglmin);
        nc       = Ntable.N_thetaH/2+1;
    }
    /* Power spectrum on logarithmic bins */
    for(i=0; i<Ntable.N_thetaH; i++) {
        l     = exp(lnrc+(i-nc)*dlnl);
        lP[i] = l*C_tomo(l,ni,nj);
    }
    /* go to log-Fourier-space */
    fftw_execute(plan);
    arg[0] = 0;   /* bias */
    for (count=0; count<=1; count++) {
        arg[1] = (count==0 ? 0 : 4);   /* order of Bessel function */
        /* perform the convolution, negative sign for kernel (complex conj.!) */
        for(i=0; i<Ntable.N_thetaH/2+1; i++) {
            kk = 2*constants.pi*i/(dlnl*Ntable.N_thetaH);
            hankel_kernel_FT(kk, &kernel, arg, 2);
            conv[i][0] = f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
            conv[i][1] = f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
        }
        /* force Nyquist- and 0-frequency-components to be double */
        conv[0][1] = 0;
        conv[Ntable.N_thetaH/2][1] = 0;
        /* go back to double space, i labels log-bins in theta */
        fftw_execute(plan1);
        for(i=0; i<Ntable.N_thetaH; i++) {
            t = exp((nc-i)*dlnl-lnrc);             /* theta=1/l */
            xi[count][Ntable.N_thetaH-i-1] = lP[i]/(t*2*constants.pi*Ntable.N_thetaH);
        }
    }

    *logthetamin = (nc-Ntable.N_thetaH+1)*dlnl-lnrc;
    *logthetamax = nc*dlnl-lnrc;
    /* clean up */
    fftw_free(conv);
    fftw_free(lP);
    fftw_free(f_lP);
    fftw_destroy_plan(plan);
    fftw_destroy_plan(plan1);
}

/****************  Hankel transform routine used in cosmo2D_xx ****************/
void hankel_kernel_FT(double x, fftw_complex *res, double *arg, int argc)
{
    fftw_complex a1, a2, g1, g2;
    int           mu;
    double        mod, xln2, si, co, d1, d2, pref, q;
    q = arg[0];
    mu = (int)(arg[1]+0.1);

    /* arguments for complex gamma */
    a1[0] = 0.5*(1.0+mu+q);
    a2[0] = 0.5*(1.0+mu-q);
    a1[1] = 0.5*x; a2[1]=-a1[1];
    cdgamma(a1,&g1);
    cdgamma(a2,&g2);
    xln2 = x*constants.ln2;
    si   = sin(xln2);
    co   = cos(xln2);
    d1   = g1[0]*g2[0]+g1[1]*g2[1]; /* Re */
    d2   = g1[1]*g2[0]-g1[0]*g2[1]; /* Im */
    mod  = g2[0]*g2[0]+g2[1]*g2[1];
    pref = exp(constants.ln2*q)/mod;

    (*res)[0] = pref*(co*d1-si*d2);
    (*res)[1] = pref*(si*d1+co*d2);
}



/*double theta[9] = {0.00020751, 0.0004224 , 0.00085981, 0.0017502 , 0.00356264,
                   0.00725196, 0.01476178, 0.03004846, 0.06116538};
double theta_plus_one[9] = {0.0004224 , 0.00085981, 0.0017502 , 0.00356264,
                            0.00725196, 0.01476178, 0.03004846, 0.06116538,0.10181087 };
// made up the last number to check things

for(i=0; i<like.Ntheta ; i++){
    theta[i] = exp(log(like.vtmin)+(i+0.0)*logdt);
    theta_plus_one[i] = exp(log(like.vtmin)+(i+1.0)*logdt); //equivalent to theta[i+1]

    xmin[i]=cos(theta[i]);
    xmax[i]=cos(theta_plus_one[i]);
    //printf("XIPM: theta[i]= %le, theta[i+1] = %le, xmin[i] = %le, xmax[i] = %le \n",
    //                        i, theta[i], theta_plus_one[i], xmin[i], xmax[i]);
}*/