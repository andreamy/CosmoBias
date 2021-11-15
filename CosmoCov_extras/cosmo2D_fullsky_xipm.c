double w_tomo_fullsky(int nt, int ni, int nj);
// galaxy clustering tomography 2PCF galaxies in bins ni, nj, computed as sum P_l(cos(like.theta[nt])*C(l,ni,nj)
double w_gamma_t_fullsky(int nt,int ni, int nj); //G-G lensing, lens bin ni, source bin nj, including IA contamination if like.IA = 3
double xi_pm_fullsky_FINEBINNING(int pm, double *theta, int nt, int ni, int nj); //shear tomography correlation functions, including IA contamination if like.IA = 3
double xi_pm_fullsky(int pm, int nt, int ni, int nj);
double xi_pm_tomo(int pm, double theta, int ni, int nj); //shear tomography correlation functions, including IA contamination if like.IA = 3
double xi_pm_rebin(int pm, double thetamin_i, double thetamax_i, int ni,int nj);//xi_pm averaged over large bins

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
