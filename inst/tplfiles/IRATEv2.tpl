//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--
//
// Instantaneous rates tag return models 
//  
//
//  Gary Nelson, Massachusetts Division of Marine Fisheries
//  Version 1.0.0
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--
TOP_OF_MAIN_SECTION
// set buffer sizes
 arrmblsize=5000000;
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);
 gradient_structure::set_MAX_NVAR_OFFSET(50000);
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);
 time(&start); //this is to see how long it takes to run
 cout << endl << "Start time : " << ctime(&start) << endl;
GLOBALS_SECTION
 #include <admodel.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
DATA_SECTION
// Model type 0= age-independent; 1=age-depedent
 init_int modeltype; 
// Tag recoveries model type 0= harvest only 1= harv,catch/release
 init_int tagtype;
// Starting and ending year of the release year 
 init_int styrR;
 init_int endyrR;
//Starting and ending year of recovery years
 init_int styr;
 init_int endyr;
  //Number of ages always 1 for ageindependent
 init_int nages;
 init_int firage; // Just a place holder
//Total Releases by Year
 init_matrix trelease(1,nages,styrR,endyrR);
//Recapture Matrix for harvest fish
 init_3darray rh(1,nages,styrR,endyrR,styr,endyr);
//Recapture Matrix for releases fish
 init_3darray rr(1,nages,styrR,endyrR,styr,endyr);
//Hooking Mortality
 init_vector h(styr,endyr); 
//Number of Natural Mortality Periods and Beginnng Years
 init_int mperiod;
 init_int mp;
 init_ivector mp_int(1,mp);
 init_vector mp_value(1,mp);
 int pp;
 //Estimate Nonmixing
 init_int mixperiod;
//Number of Fishing  Mortality Periods and Beginning Years
 init_int fperiod;
 init_int fp;
 init_vector forig(1,2);
 init_ivector fp_int(1,fp);
 init_vector fp_value(1,fp);
 int qq;
 //Number of Tag Mortality Periods
 init_int ftperiod;
 init_int fap;
 init_vector ftagorig(1,2);
 init_ivector fap_int(1,fap);
 init_vector fap_value(1,fap);
 int ss;
 int tp;
 init_int combl;
 //Estimate Reporting Rate Harvest
 init_int l1period;
 init_int lp1;
 init_ivector lp1_int(1,lp1);
 init_vector  lp1_value(1,lp1);
 int ll1;
 //Estimate Reporting Rate Release
 init_int l2period;
 init_int lp2;
 init_ivector lp2_int(1,lp2);
 init_vector  lp2_value(1,lp2);
 int ll2;
 //Estimate Phi  Harvest
 init_int combp;
 init_int phi1period;
 init_int phi1p;
 init_ivector phi1p_int(1,phi1p);
 init_vector  phi1p_value(1,phi1p);
 int phi1l;
 //Estimate Phi Release
 init_int phi2period;
 init_int phi2p;
 init_ivector phi2p_int(1,phi2p);
 init_vector  phi2p_value(1,phi2p);
 int phi2l;
 //Nonmixing values
 init_int nmph;
 init_ivector nonmix_int(1,nmph);
 init_vector  nonmix_value(1,nmph);
 int nonmixl;
 //Harvest selectivity
 init_int comsel;
 init_int sel1period;
 init_vector sel1_value(1,nages);
 //Release selectivity
 init_int sel2period;
 init_vector sel2_value(1,nages);
 //Natural mortality Age Selectivity
 //Estimate Phi Release
 init_int mselperiod;
 init_int mselp;
 init_ivector msel_int(1,mselp);
 init_vector  msel_value(1,mselp);
 int msell;
 int mvec;
 int bigmvec;
 int whichmperiod;
 // Use likeilhood constant
 init_int use_constant;
 //Read boundry limits
 init_number flow;
 init_number fup;
 init_number mlow;
 init_number mup;
 init_number ftlow;
 init_number ftup;
 init_number l1low;
 init_number l1up;
 init_number l2low;
 init_number l2up;
 init_number p1low;
 init_number p1up;
 init_number p2low;
 init_number p2up;
 init_number sel1low;
 init_number sel1up;
 init_number sel2low;
 init_number sel2up;
 init_number nonmixlow;
 init_number nonmixup;
 int nmpr;
 int mix
 int fp2;
 int nsel1;
 int nsel2;
 int nonmixh;
 int nonmixr;
  LOCAL_CALCS
   tp=0;
   mix=0;
  // Incomplete mixing
    if(mixperiod>0){
      if(nonmix_int(1)==styrR) mix=1;
       if(tagtype==0){
         tp=tp+nmph;
         nmpr=1;
         nonmixr=-1;
         nonmixh=mixperiod;
       }
     if(tagtype==1){
      nmpr=nmph;
      tp=tp+nmph+nmpr;
      nonmixr=mixperiod;
      nonmixh=mixperiod;
     }
    }
   if(mixperiod<=0){
     nmph=1;
     nmpr=1;
     nonmixh=mixperiod;
     nonmixr=mixperiod;
   }
   //Selectivity if agedependent
    nsel1=0;
    nsel2=0;
    for(t=1;t<=nages;t++){
       if (sel1_value(t)<1.) nsel1+=1;
       if (tagtype==1){
           if (sel2_value(t)<1.) nsel2+=1;
       }
    }
  if(mselperiod>0.) mvec=mselp;
  if(mselperiod<=0.) mvec=1;
  if(mperiod>0.) bigmvec=mp;
  if(mperiod<=0.) bigmvec=1;
  if(mselperiod>0 && mperiod>0.) whichmperiod=mperiod;
  if(mselperiod<=0 && mperiod>0.) whichmperiod=mperiod;
  if(mselperiod<=0 && mperiod<=0.) whichmperiod=mperiod;
  if(mselperiod>0 && mperiod<=0.) whichmperiod=mselperiod;
  if(sel1period>0) tp=tp+nsel1;
  if(sel2period>0) tp=tp+nsel2;
  pp=mp+1;
  qq=fp+1;
  ss=fap+1;
  ll1=lp1+1;
  ll2=lp2+1;
  phi1l=phi1p+1;
  phi2l=phi2p+1;
  msell=mselp+1;
  nonmixl=nmph+1;
   if(mperiod>0) tp=tp+mp;
   if(fperiod>0) tp=tp+fp;
  if(ftperiod>0) tp=tp+fap;
  if(l1period>0) tp=tp+lp1;
  if(l2period>0) tp=tp+lp2;
  if(phi1period>0) tp=tp+phi1p;
  if(phi2period>0) tp=tp+phi2p;
  if(mselperiod>0) tp=tp+mselp;
   //Survival, FM, FT, NM, L1, L2, P1, P2,Nonmx1, nonmix2
     tp+=nages*(endyr-styr+1);//Survival
     tp+=(endyr-styr+1);//FM
     tp+=(endyr-styr+1);//FT
     tp+=nages*(endyr-styr+1);//NM
     tp+=(endyr-styr+1);//LAMBDA1
     tp+=(endyr-styr+1);//LAMBDA2
     tp+=(endyrR-styrR+1);//PHI1
     tp+=(endyrR-styrR+1);//PHI2
     tp+=(endyrR-styrR+1);//FNMIXh
     tp+=(endyrR-styrR+1);//FNMIXr
 END_CALCS
 matrix sigma(1,tp,1,tp+1); 
 !! set_covariance_matrix(sigma); 
 //looping variables
 int y;
 int t;
 int a;
 int d;
 int cnt;
 int total;
 int Ntags;
 int looper;
 int df_r;
 int df_h;
 int hless;
 int rless;
 int kk;
 int age;
PARAMETER_SECTION
 number dodo;
 number dodo1;
 number probs;
 number K;
 number up_df;
 number up_count;
 number up_chi;
 number up_chat;
 number p_chi;
 number p_df;
 number p_chat;
 number constant;
 number nfact_out;
 number nfact_in;
 number ns_count;
 number rless;
 number hless;
 //-------------F estimates-------------------------------
 init_bounded_vector e_F(1,fp,flow,fup,fperiod);
 matrix F(styrR,endyrR,styr,endyr);
 vector fp_yr(1,qq);
  //-------------M estimates-------------------------------
 init_bounded_matrix e_M(1,bigmvec,1,mvec,mlow,mup,whichmperiod);
 matrix M(styr,endyr,1,nages);
 vector mp_yr(1,pp);
 //--------------Tag Mortality----------------------------
 init_bounded_vector e_FA(1,fap,ftlow,ftup,ftperiod);
 matrix FA(styrR,endyrR,styr,endyr);
 vector fap_yr(1,ss);
 //--------------Reporting Rate----------------------------
 init_bounded_vector LRR1(1,lp1,l1low,l1up,l1period);
 vector lr1(styr,endyr);
 vector lp1_yr(1,ll1);
 init_bounded_vector LRR2(1,lp2,l2low,l2up,l2period);
 vector lr2(styr,endyr);
 vector lp2_yr(1,ll2);
 //--------------Phi----------------------------
 init_bounded_vector PHI1R(1,phi1p,p1low,p1up,phi1period);
 vector phi1rr(styrR,endyrR);
 vector phi1p_yr(1,phi1l);
 init_bounded_vector PHI2R(1,phi2p,p2low,p2up,phi2period);
 vector phi2rr(styrR,endyrR);
 vector phi2p_yr(1,phi2l);
 //--------------NonMixing----------------------------
 init_bounded_vector e_nmixh(1,nmph,nonmixlow,nonmixup,nonmixh);
 init_bounded_vector e_nmixr(1,nmpr,nonmixlow,nonmixup,nonmixr);
  vector nonmix_yr(1,nonmixl);
//--------------Catch Age Selectivity Estimates-------------------
 init_bounded_vector sel1(1,nsel1,sel1low,sel1up,sel1period);
 init_bounded_vector sel2(1,nsel2,sel2low,sel2up,sel2period);
 //-----------------M age slecitivity-------------------------
 vector msel_age(1,msell);
 //------------------Number of Tags-------------------
 matrix tags(styrR,endyrR,1,nages);
 matrix N(styrR,endyrR,1,nages);
 //----------------Mortality Calculations------------
 3darray s(1,nages,styrR,endyrR,styr,endyr);
 3darray u_h(1,nages,styrR,endyrR,styr,endyr);
 3darray u_r(1,nages,styrR,endyrR,styr,endyr);
 //---------Predicted Matrices-------------
 matrix sum_prob_h(styrR,endyrR,1,nages);
 matrix sum_prob_r(styrR,endyrR,1,nages);
 3darray s_prob(1,nages,styrR,endyrR,styr,endyr);
 3darray exp_prob_h(1,nages,styrR,endyrR,styr,endyr);
 3darray ll_h(1,nages,styrR,endyrR,styr,endyr);
 3darray exp_prob_r(1,nages,styrR,endyrR,styr,endyr);
 3darray ll_r(1,nages,styrR,endyrR,styr,endyr);
 matrix  ll_ns(styrR,endyrR,1,nages);
 3darray exp_r_h(1,nages,styrR,endyrR,styr,endyr);
 3darray exp_r_r(1,nages,styrR,endyrR,styr,endyr);
 3darray pool_r(1,nages,styrR,endyrR,styr,endyr);
 3darray pool_h(1,nages,styrR,endyrR,styr,endyr);
 3darray pool_r_e(1,nages,styrR,endyrR,styr,endyr);
 3darray pool_h_e(1,nages,styrR,endyrR,styr,endyr);
 3darray chi_r(1,nages,styrR,endyrR,styr,endyr);
 3darray chi_h(1,nages,styrR,endyrR,styr,endyr);
 3darray p_chi_r(1,nages,styrR,endyrR,styr,endyr);
 3darray p_chi_h(1,nages,styrR,endyrR,styr,endyr);
 3darray pear_r(1,nages,styrR,endyrR,styr,endyr);
 3darray pear_h(1,nages,styrR,endyrR,styr,endyr);
 3darray stdres_r(1,nages,styrR,endyrR,styr,endyr);
 3darray stdres_h(1,nages,styrR,endyrR,styr,endyr);
 matrix exp_ns(styrR,endyrR,1,nages);
 matrix chi_ns(styrR,endyrR,1,nages);
 matrix pear_ns(styrR,endyrR,1,nages);
 matrix stdres_ns(styrR,endyrR,1,nages);
 vector hsel(1,nages);
 vector rsel(1,nages);
 sdreport_matrix S(styr,endyr,1,nages);
 sdreport_vector FM(styr,endyr);
 sdreport_vector FT(styr,endyr);
 sdreport_matrix NM(styr,endyr,1,nages);
 sdreport_vector LAMBDA1(styr,endyr);
 sdreport_vector LAMBDA2(styr,endyr);
 sdreport_vector PHI1OUT(styrR,endyrR);
 sdreport_vector PHI2OUT(styrR,endyrR);
 sdreport_vector FNMIXh(styrR,endyrR);
 sdreport_vector FNMIXr(styrR,endyrR);
 sdreport_vector SEL1OUT(1,nages);
 sdreport_vector SEL2OUT(1,nages);
 //----------Likelihood Values------------------------------------------
 number f_tag;
 objective_function_value f;
 
INITIALIZATION_SECTION

RUNTIME_SECTION
 maximum_function_evaluations 100, 500, 5000;
 convergence_criteria 1e-5, 1e-7, 1e-16;
PRELIMINARY_CALCS_SECTION
 F.initialize();
 FA.initialize();
 M.initialize();
 LAMBDA1.initialize();
 LAMBDA2.initialize();
 PHI1OUT.initialize();
 PHI2OUT.initialize();
 if(mperiod>0. && mselperiod>0.){
  for(t=1;t<=mp;t++){
    for(a=1;a<=mselp;a++){
      e_M(t,a)=log(mp_value(t)*msel_value(a));
    }
  }
 }
 if(mperiod>0. && mselperiod<=0.){
  for(t=1;t<=mp;t++){
    for(a=1;a<=mvec;a++){
      e_M(t,a)=log(mp_value(t));
    }
  }
 }
  if(mperiod<=0. && mselperiod>0.){
  for(t=1;t<=bigmvec;t++){
    for(a=1;a<=mselp;a++){
      e_M(t,a)=log(msel_value(a));
    }
  }
 }

 for(t=1;t<=fp;t++){
   e_F(t)=log(fp_value(t));
 }
 for(t=1;t<=fap;t++){
   e_FA(t)=log(fap_value(t));
 }
 for(t=1;t<=lp1;t++){
   LRR1(t)=lp1_value(t);
 }
 for(t=1;t<=lp2;t++){
   LRR2(t)=lp2_value(t);
 }
 for(t=1;t<=phi1p;t++){
   PHI1R(t)=phi1p_value(t);
 }
 for(t=1;t<=phi2p;t++){
   PHI2R(t)=phi2p_value(t);
 }
 for(t=1;t<=nmph;t++){
  e_nmixh(t)=log(nonmix_value(t));
  if (tagtype==1) e_nmixr(t)=log(nonmix_value(t));
 }
  cnt=0;
  for(t=1;t<=nages;t++){
       if (sel1_value(t)<1.){ 
          cnt+=1;
          sel1(cnt)=sel1_value(t);
        }
  }
  if(tagtype==1){
    cnt=0;
    for(t=1;t<=nages;t++){
       if (sel2_value(t)<1.){ 
          cnt+=1;
          sel2(cnt)=sel2_value(t);
        }
    }
  }
  N=trans(trelease);
  //if missing a year of data
 for(a=1;a<=nages;a++){
  for(t=styrR;t<=endyrR;t++){
    if(N(t,a)<0) N(t,a)=0;
  }
 }

 
 //Count tags
  if(tagtype==1){
    for(a=1;a<=nages;a++){
      cnt=0;
      for(t=styrR;t<=endyrR;t++){
        Ntags=0;
        for (y=styr+cnt;y<=endyr;y++){
           if(rh(a,t,y)>=0)  Ntags+=rh(a,t,y);
           if(rr(a,t,y)>=0)  Ntags+=rr(a,t,y);
        }
        tags(t,a)=Ntags;
        cnt+=1;
      }
    }
  }

 if(tagtype==0){
   for(a=1;a<=nages;a++){
     cnt=0;
     for (t=styrR;t<=endyrR;t++){
        Ntags=0;
        for(y=styr+cnt;y<=endyr;y++){
         if(rh(a,t,y)>=0)  Ntags+=rh(a,t,y);
        }
        tags(t,a)=Ntags;
       cnt+=1;
      }
   }
 }


 for(y=styrR;y<=endyrR;y++){
   FNMIXh(y)=-1;
   FNMIXr(y)=-1;
 }

 if(use_constant>0){
 constant=0.0;
  for(a=1;a<=nages;a++){
   for(t=styrR;t<=endyrR;t++){
    //N factorial
     nfact_in=N(t,a);
       nfact_out=0.0;
         if(nfact_in>=2){
           for(kk=2;kk<=nfact_in;kk++) nfact_out+=log(kk);
         }
     constant+=nfact_out;
   //never seen factorial
      nfact_in=N(t,a)-tags(t,a);
      nfact_out=0.0;
         if(nfact_in>=2){
           for(kk=2;kk<=nfact_in;kk++) nfact_out+=log(kk);
         }
     constant+=-1.0*nfact_out;
     for(y=styr;y<=endyr;y++){
       nfact_in=rh(a,t,y);
        nfact_out=0.0;
         if(nfact_in>=2){
           for(kk=2;kk<=nfact_in;kk++) nfact_out+=log(kk);
         }
       constant+=-1.0*nfact_out;

       if(tagtype==1){
         nfact_in=rr(a,t,y);
         nfact_out=0.0;
         if(nfact_in>=2){
           for(kk=2;kk<=nfact_in;kk++) nfact_out+=log(kk);
         }
         constant+=-1.0*nfact_out;
        
       }
    }
  }
 }
 }
PROCEDURE_SECTION
 fill_sdreports();
 calc_M_vector();
 calc_F_vector();
 calc_FA_vector();
 calc_RR1_vector();
 calc_RR2_vector();
 calc_PHI1_vector();
 calc_PHI2_vector();
 calc_selectivity();
 calc_fish_surv();
 calc_s();
 calc_s_prob();
 calc_u_h();
 calc_u_r();
 calc_exp_prob_h();
 calc_exp_prob_r();
 calc_LL();
 evaluate_the_objective_function();

FUNCTION fill_sdreports    
     for(y=styrR;y<=endyrR;y++){
      if(nonmixh<=0){
        if(mperiod>0) FNMIXh(y)=mfexp(e_M(1,1))*-1.;
        if(fperiod>0) FNMIXh(y)=mfexp(e_F(1))*-1.;
        if(ftperiod>0) FNMIXh(y)=mfexp(e_FA(1))*-1.;
        if(l1period>0) FNMIXh(y)=LRR1(1)*-1.;
        if(l2period>0) FNMIXh(y)=LRR2(1)*-1.;
        if(phi1period>0) FNMIXh(y)=PHI1R(1)*-1.;
        if(phi2period>0) FNMIXh(y)=PHI2R(1)*-1.;
        if(mselperiod>0) FNMIXh(y)=mfexp(e_M(1,1))*-1.;
        if(sel1period>0) FNMIXh(y)=mfexp(sel1(1))*-1.;
        if(sel2period>0) FNMIXh(y)=mfexp(sel2(1))*-1.;
      }
      if(nonmixr<=0){
        if(mperiod>0) FNMIXr(y)=mfexp(e_M(1,1))*-1.;
        if(fperiod>0) FNMIXr(y)=mfexp(e_F(1))*-1.;
        if(ftperiod>0) FNMIXr(y)=mfexp(e_FA(1))*-1.;
        if(l1period>0) FNMIXr(y)=LRR1(1)*-1.;
        if(l2period>0) FNMIXr(y)=LRR2(1)*-1.;
        if(phi1period>0) FNMIXr(y)=PHI1R(1)*-1.;
        if(phi2period>0) FNMIXr(y)=PHI2R(1)*-1.;
        if(mselperiod>0) FNMIXr(y)=mfexp(e_M(1,1))*-1.;
        if(sel1period>0) FNMIXr(y)=mfexp(sel1(1))*-1.;
        if(sel2period>0) FNMIXr(y)=mfexp(sel2(1))*-1.;
     }
  }
FUNCTION calc_M_vector
 for(t=1;t<=mp;t++) mp_yr(t)=mp_int(t);  
   mp_yr(pp)=endyr+1;
  for(t=1;t<=mselp;t++) msel_age(t)=msel_int(t);  
      msel_age(msell)=nages+1;
 
 for(a=1;a<=nages;a++){
   for(t=styr;t<=endyr;t++){
     for(d=1;d<=mp;d++){
      for(kk=1;kk<=mselp;kk++){
         if(t>=mp_yr(d) && t<mp_yr(d+1)){
           if(a>=msel_age(kk) && a<msel_age(kk+1)){
             if(mperiod>0. && mselperiod>0.){ 
              M(t,a)=mfexp(e_M(d,kk));
              NM(t,a)=M(t,a);
             }
           if(mperiod>0. && mselperiod<=0.){ 
             M(t,a)=mfexp(e_M(d,1))*msel_value(kk);
             NM(t,a)=M(t,a);
            }
           if(mperiod<=0. && mselperiod>0.){ 
             M(t,a)=mp_value(d)*mfexp(e_M(1,kk));
             NM(t,a)=M(t,a);
            }
           if(mperiod<=0. && mselperiod<=0.){ 
             M(t,a)=mp_value(d)*msel_value(kk);
             if(fperiod>0) NM(t,a)=mfexp(e_F(1))*-1.;
             if(ftperiod>0) NM(t,a)=mfexp(e_FA(1))*-1.;
             if(l1period>0) NM(t,a)=LRR1(1)*-1.;
             if(l2period>0) NM(t,a)=LRR2(1)*-1.;
             if(phi1period>0) NM(t,a)=PHI1R(1)*-1.;
             if(phi2period>0) NM(t,a)=PHI1R(1)*-1.;
             if(sel1period>0) NM(t,a)=sel1(1)*-1.;
             if(sel2period>0) NM(t,a)=sel2(1)*-1.;
           }
         }
       }
   }
 }
 }
 }

FUNCTION calc_F_vector
    for(t=1;t<=fp;t++) fp_yr(t)=fp_int(t);  
     fp_yr(qq)=endyr+1;
     for(t=styrR;t<=endyrR;t++){
       for(y=styr+mix;y<=endyr;y++){
       for(d=1;d<=fp;d++){
         if(y>=fp_yr(d) && y<fp_yr(d+1)){
           if(fperiod>0){
              F(t,y)=mfexp(e_F(d));
              FM(y)=mfexp(e_F(d));
          }
           if(fperiod<=0){
            F(t,y)=fp_value(d);  
             if(mperiod>0) FM(y)=mfexp(e_M(1,1))*-1.;
             if(ftperiod>0) FM(y)=mfexp(e_FA(1))*-1.;
             if(l1period>0)  FM(y)=LRR1(1)*-1.;
             if(l2period>0)  FM(y)=LRR2(1)*-1.;
             if(phi1period>0)  FM(y)=PHI1R(1)*-1.;
             if(phi2period>0)  FM(y)=PHI2R(1)*-1.;
             if(sel1period>0) FM(y)=sel1(1)*-1.;
             if(sel2period>0) FM(y)=sel2(1)*-1.;
             if(mselperiod>0) FM(y)=mfexp(e_M(1,1))*-1.;
            }
           }
         }
        }
      }
 if(fperiod>0 && mix==1) FM(styr)=mfexp(e_F(1))*-1.;
  if(fperiod<=0 && mix==1) {
    if(mperiod>0) FM(styr)=mfexp(e_M(1,1))*-1.;
    if(ftperiod>0) FM(styr)=mfexp(e_FA(1))*-1.;
    if(l1period>0) FM(styr)=LRR1(1)*-1.;
    if(l2period>0) FM(styr)=LRR2(1)*-1.;
    if(phi1period>0) FM(styr)=PHI1R(1)*-1.;
    if(phi2period>0) FM(styr)=PHI2R(1)*-1.;
    if(sel1period>0) FM(styr)=sel1(1)*-1.;
    if(sel2period>0) FM(styr)=sel2(1)*-1.;
    if(mselperiod>0) FM(styr)=mfexp(e_M(1,1))*-1.;
   }
 if(nonmixh>0){
     // Nonmiximg code
    // Fill diagonals in F matrix
    for(t=1;t<=nmph;t++) nonmix_yr(t)=nonmix_int(t);  
     nonmix_yr(nonmixl)=endyrR+1;
     for(t=styrR;t<=endyrR;t++){
       if(endyr==endyrR){
        for(y=styr;y<=endyr;y++){
          for(d=1;d<=nmph;d++){
           if(y>=nonmix_yr(d) && y<nonmix_yr(d+1)){
             if(t==y){
                F(t,y)=mfexp(e_nmixh(d));
                FNMIXh(y)=mfexp(e_nmixh(d));
             }
           }
        }
       }
     }
    if(endyr>endyrR){
       for(y=styr;y<=endyrR;y++){
          for(d=1;d<=nmph;d++){
           if(y>=nonmix_yr(d) && y<nonmix_yr(d+1)){
             if(t==y){
                F(t,y)=mfexp(e_nmixh(d));
                FNMIXh(y)=mfexp(e_nmixh(d));
             }
           }
        }
     }
   }
 }
 }//if mixperiod>0
FUNCTION calc_FA_vector
 
   for(t=1;t<=fap;t++) fap_yr(t)=fap_int(t);
    fap_yr(ss)=endyr+1;
    for(t=styrR;t<=endyrR;t++){
      for(y=styr+mix;y<=endyr;y++){
        for(d=1;d<=fap;d++){
         if(y>=fap_yr(d) && y<fap_yr(d+1)){
           if(ftperiod>0){
                FA(t,y)=mfexp(e_FA(d));
                FT(y)=FA(t,y);
          }
           if(ftperiod<=0){
             FA(t,y)=fap_value(d);
             if(mperiod>0) FT(y)=mfexp(e_M(1,1))*-1.;
             if(fperiod>0) FT(y)=mfexp(e_F(1))*-1.;
             if(l1period>0) FT(y)=LRR1(1)*-1.;
             if(l2period>0) FT(y)=LRR2(1)*-1.;
             if(phi1period>0) FT(y)=PHI1R(1)*-1.;
             if(phi2period>0) FT(y)=PHI2R(1)*-1.;
             if(sel1period>0) FT(y)=sel1(1)*-1.;
             if(sel2period>0) FT(y)=sel2(1)*-1.;
              if(mselperiod>0) FT(y)=mfexp(e_M(1,1))*-1.;
           }
          }
       }
      }
   }
  if(ftperiod>0 && mix==1) FT(styr)=mfexp(e_FA(1))*-1.;
  if(ftperiod<=0 && mix==1) {
    if(mperiod>0) FT(styr)=mfexp(e_M(1,1))*-1.;
    if(fperiod>0) FT(styr)=mfexp(e_F(1))*-1.;
    if(l1period>0) FT(styr)=LRR1(1)*-1.;
    if(l2period>0) FT(styr)=LRR2(1)*-1.;
    if(phi1period>0) FT(styr)=PHI1R(1)*-1.;
    if(phi2period>0) FT(styr)=PHI2R(1)*-1.;
    if(sel1period>0) FT(styr)=sel1(1)*-1.;
    if(sel2period>0) FT(styr)=sel2(1)*-1.;
    if(mselperiod>0) FT(styr)=mfexp(e_M(1,1))*-1.;
   }

 if(nonmixr>0){
       // Fill diagonals in F matrix
    for(t=1;t<=nmpr;t++) nonmix_yr(t)=nonmix_int(t);  
     nonmix_yr(nonmixl)=endyrR+1;
     for(t=styrR;t<=endyrR;t++){
       if(endyr==endyrR){
        for(y=styr;y<=endyr;y++){
          for(d=1;d<=nmpr;d++){
           if(y>=nonmix_yr(d) && y<nonmix_yr(d+1)){
             if(t==y){
                FA(t,y)=mfexp(e_nmixr(d));
                FNMIXr(y)=mfexp(e_nmixr(d));
             }
           }
        }
       }
     }
    if(endyr>endyrR){
       for(y=styr;y<=endyrR;y++){
          for(d=1;d<=nmpr;d++){
           if(y>=nonmix_yr(d) && y<nonmix_yr(d+1)){
             if(t==y){
                FA(t,y)=mfexp(e_nmixr(d));
                FNMIXr(y)=mfexp(e_nmixr(d));
             }
           }
        }
     }
   }
 }
 }

FUNCTION calc_RR1_vector
  for(t=1;t<=lp1;t++) lp1_yr(t)=lp1_int(t);
    lp1_yr(ll1)=endyr+1;
 for(t=styr;t<=endyr;t++) {
     for(d=1;d<=lp1;d++){
         if(t>=lp1_yr(d) && t<lp1_yr(d+1)){
            if(l1period>0){
                lr1(t)=LRR1(d);
               LAMBDA1(t)=lr1(t);
             }
            if(l1period<=0){
              lr1(t)=lp1_value(d);
              if(mperiod>0) LAMBDA1(t)=mfexp(e_M(1,1))*-1.;
              if(fperiod>0) LAMBDA1(t)=mfexp(e_F(1))*-1.;
              if(ftperiod>0) LAMBDA1(t)=mfexp(e_FA(1))*-1.;
              if(l2period>0) LAMBDA1(t)=LRR2(1)*-1.;
              if(phi1period>0) LAMBDA1(t)=PHI1R(1)*-1.;
              if(phi2period>0) LAMBDA1(t)=PHI2R(1)*-1.;
              if(sel1period>0) LAMBDA1(t)=sel1(1)*-1.;
             if(sel2period>0)  LAMBDA1(t)=sel2(1)*-1.;
              if(mselperiod>0) LAMBDA1(t)=mfexp(e_M(1,1))*-1.;
            }
         }
       }
   }

FUNCTION calc_RR2_vector
 if(combl<=0){
   for(t=1;t<=lp2;t++) lp2_yr(t)=lp2_int(t);
   lp2_yr(ll2)=endyr+1;
   for(t=styr;t<=endyr;t++) {
     for(d=1;d<=lp2;d++){
         if(t>=lp2_yr(d) && t<lp2_yr(d+1)){
            if(l2period>0){
                lr2(t)=LRR2(d);
               LAMBDA2(t)=lr2(t);
             }
            if(l2period<=0){
              lr2(t)=lp2_value(d);
              if(mperiod>0) LAMBDA2(t)=mfexp(e_M(1,1))*-1.;
              if(fperiod>0) LAMBDA2(t)=mfexp(e_F(1))*-1.;
              if(ftperiod>0) LAMBDA2(t)=mfexp(e_FA(1))*-1.;
              if(l1period>0) LAMBDA2(t)=LRR1(1)*-1.;
              if(phi1period>0) LAMBDA2(t)=PHI1R(1)*-1.;
              if(phi2period>0) LAMBDA2(t)=PHI2R(1)*-1.;
             if(sel1period>0) LAMBDA2(t)=sel1(1)*-1.;
             if(sel2period>0)  LAMBDA2(t)=sel2(1)*-1.;
             if(mselperiod>0) LAMBDA2(t)=mfexp(e_M(1,1))*-1.;
            }
         }
       }
   }
  }
 if(combl>0){
  for(t=styr;t<=endyr;t++){
      lr2(t)=lr1(t);  
        if(mperiod>0) LAMBDA2(t)=mfexp(e_M(1,1))*-1.;
        if(fperiod>0) LAMBDA2(t)=mfexp(e_F(1))*-1.;
        if(ftperiod>0) LAMBDA2(t)=mfexp(e_FA(1))*-1.;
        if(l1period>0) LAMBDA2(t)=LRR1(1)*-1.;
        if(phi1period>0) LAMBDA2(t)=PHI1R(1)*-1.;
        if(phi2period>0) LAMBDA2(t)=PHI2R(1)*-1.;
        if(sel1period>0) LAMBDA2(t)=sel1(1)*-1.;
        if(sel2period>0)  LAMBDA2(t)=sel2(1)*-1.;
         if(mselperiod>0) LAMBDA2(t)=mfexp(e_M(1,1))*-1.;
   }
 }
FUNCTION calc_PHI1_vector
    for(t=1;t<=phi1p;t++) phi1p_yr(t)=phi1p_int(t);
    phi1p_yr(phi1l)=endyrR+1;
    for(t=styrR;t<=endyrR;t++) {
      for(d=1;d<=phi1p;d++){
         if(t>=phi1p_yr(d) && t<phi1p_yr(d+1)){
            if(phi1period>0){
                phi1rr(t)=PHI1R(d);
               PHI1OUT(t)= phi1rr(t);
             }
            if(phi1period<=0){
              phi1rr(t)=phi1p_value(d);
              if(mperiod>0) PHI1OUT(t)=mfexp(e_M(1,1))*-1.;
              if(fperiod>0)  PHI1OUT(t)=mfexp(e_F(1))*-1.;
              if(ftperiod>0) PHI1OUT(t)=mfexp(e_FA(1))*-1.;
              if(l1period>0) PHI1OUT(t)=LRR1(1)*-1.;
              if(l2period>0) PHI1OUT(t)=LRR2(1)*-1.;
               if(phi2period>0) PHI1OUT(t)=PHI2R(1)*-1.;
               if(sel1period>0) PHI1OUT(t)=sel1(1)*-1.;
               if(sel2period>0) PHI1OUT(t)=sel2(1)*-1.;
               if(mselperiod>0) PHI1OUT(t)=mfexp(e_M(1,1))*-1.;
            }
         }
       }
   }
  
FUNCTION calc_PHI2_vector
 if(combp<=0){
    for(t=1;t<=phi2p;t++) phi2p_yr(t)=phi2p_int(t);
    phi2p_yr(phi2l)=endyrR+1;
    for(t=styrR;t<=endyrR;t++) {
      for(d=1;d<=phi2p;d++){
         if(t>=phi2p_yr(d) && t<phi2p_yr(d+1)){
            if(phi2period>0){
                phi2rr(t)=PHI2R(d);
               PHI2OUT(t)= phi2rr(t);
             }
            if(phi2period<=0){
              phi2rr(t)=phi2p_value(d);
              if(mperiod>0) PHI2OUT(t)=mfexp(e_M(1,1))*-1.;
              if(fperiod>0)  PHI2OUT(t)=mfexp(e_F(1))*-1.;
              if(ftperiod>0) PHI2OUT(t)=mfexp(e_FA(1))*-1.;
              if(l2period>0) PHI2OUT(t)=LRR2(1)*-1.;
              if(l1period>0) PHI2OUT(t)=LRR1(1)*-1.;
               if(phi1period>0) PHI2OUT(t)=PHI1R(1)*-1.;
              if(sel1period>0) PHI2OUT(t)=sel1(1)*-1.;
              if(sel2period>0) PHI2OUT(t)=sel2(1)*-1.;
              if(mselperiod>0) PHI2OUT(t)=mfexp(e_M(1,1))*-1.;
            }
         }
       }
     }
  }
 if(combp>0){
   for(t=styrR;t<=endyrR;t++) {
       phi2rr(t)=phi1rr(t);
        if(mperiod>0) PHI2OUT(t)=mfexp(e_M(1,1))*-1.;
        if(fperiod>0)  PHI2OUT(t)=mfexp(e_F(1))*-1.;
        if(ftperiod>0) PHI2OUT(t)=mfexp(e_FA(1))*-1.;
        if(l2period>0) PHI2OUT(t)=LRR2(1)*-1.;
        if(l1period>0) PHI2OUT(t)=LRR1(1)*-1.;
         if(phi1period>0) PHI2OUT(t)=PHI1R(1)*-1.;
         if(sel1period>0) PHI2OUT(t)=sel1(1)*-1.;
          if(sel2period>0) PHI2OUT(t)=sel2(1)*-1.;
          if(mselperiod>0) PHI2OUT(t)=mfexp(e_M(1,1))*-1.;
     }
  }
FUNCTION calc_selectivity
 if(modeltype==1){
    cnt=0;
   for(a=1;a<=nages;a++){
         if(sel1period<=0){
              hsel(a)=sel1_value(a);
              if(mperiod>0) SEL1OUT(a)=mfexp(e_M(1,1))*-1.;
              if(fperiod>0)  SEL1OUT(a)=mfexp(e_F(1))*-1.;
              if(ftperiod>0) SEL1OUT(a)=mfexp(e_FA(1))*-1.;
              if(l2period>0) SEL1OUT(a)=LRR2(1)*-1.;
              if(l1period>0) SEL1OUT(a)=LRR1(1)*-1.;
              if(phi1period>0) SEL1OUT(a)=PHI1R(1)*-1.;
              if(phi2period>0) SEL1OUT(a)=PHI2R(1)*-1.;
              if(sel2period>0) SEL1OUT(a)=sel2(1)*-1.;
              if(mselperiod>0) SEL1OUT(a)=mfexp(e_M(1,1))*-1.;
              if(tagtype==0){
                if(mperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
                if(fperiod>0)  SEL2OUT(a)=mfexp(e_F(1))*-1.;
                if(ftperiod>0) SEL2OUT(a)=mfexp(e_FA(1))*-1.;
                if(l2period>0) SEL2OUT(a)=LRR2(1)*-1.;
                if(l1period>0) SEL2OUT(a)=LRR1(1)*-1.;
                if(phi1period>0) SEL2OUT(a)=PHI1R(1)*-1.;
                if(phi2period>0) SEL2OUT(a)=PHI2R(1)*-1.;
                if(mselperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
             }  
         }
        if(sel1period>0){
           if(sel1_value(a)<1.){
              cnt+=1.;
              hsel(a)=sel1(cnt);
              SEL1OUT(a)=sel1(cnt);
            }
           if(sel1_value(a)>=1.){
            hsel(a)=1.;
            SEL1OUT(a)=sel1(cnt)*-1;
           }
            if(tagtype==0){
                if(mperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
                if(fperiod>0)  SEL2OUT(a)=mfexp(e_F(1))*-1.;
                if(ftperiod>0) SEL2OUT(a)=mfexp(e_FA(1))*-1.;
                if(l2period>0) SEL2OUT(a)=LRR2(1)*-1.;
                if(l1period>0) SEL2OUT(a)=LRR1(1)*-1.;
                if(phi1period>0) SEL2OUT(a)=PHI1R(1)*-1.;
                if(phi2period>0) SEL2OUT(a)=PHI2R(1)*-1.;
                if(sel1period>0) SEL2OUT(a)=sel1(1)*-1.;
                if(mselperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
             }  
         }
      }
             
   if(tagtype==1){
     cnt=0;
     for(a=1;a<=nages;a++){
          if(sel2period<=0){
            if(comsel<=0){
              rsel(a)=sel2_value(a);
               if(mperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
              if(fperiod>0)  SEL2OUT(a)=mfexp(e_F(1))*-1.;
              if(ftperiod>0) SEL2OUT(a)=mfexp(e_FA(1))*-1.;
              if(l2period>0) SEL2OUT(a)=LRR2(1)*-1.;
              if(l1period>0) SEL2OUT(a)=LRR1(1)*-1.;
              if(phi1period>0) SEL2OUT(a)=PHI1R(1)*-1.;
              if(phi2period>0) SEL2OUT(a)=PHI2R(1)*-1.;
              if(sel1period>0) SEL2OUT(a)=sel1(1)*-1.;
              if(mselperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
            }
            if(comsel>0){
               rsel(a)=hsel(a);
               if(mperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
              if(fperiod>0)  SEL2OUT(a)=mfexp(e_F(1))*-1.;
              if(ftperiod>0) SEL2OUT(a)=mfexp(e_FA(1))*-1.;
              if(l2period>0) SEL2OUT(a)=LRR2(1)*-1.;
              if(l1period>0) SEL2OUT(a)=LRR1(1)*-1.;
              if(phi1period>0) SEL2OUT(a)=PHI1R(1)*-1.;
              if(phi2period>0) SEL2OUT(a)=PHI2R(1)*-1.;
              if(sel1period>0) SEL2OUT(a)=sel1(1)*-1.;
              if(mselperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
            }
       }     
         if(sel2period>0){
           if(sel2_value(a)<1.){
              cnt+=1.;
              rsel(a)=sel2(cnt);
              SEL2OUT(a)=sel2(cnt);
            }
           if(sel2_value(a)>=1.){
             rsel(a)=1.;
             SEL2OUT(a)=sel2(cnt)*-1;
           }
         
         }
     }
   }//tagtype=1
 }//modetype==1
 if(modeltype==0){
   for(a=1;a<=nages;a++){
     hsel(a)=1.;
     rsel(a)=1.;
     if(mperiod>0) SEL1OUT(a)=mfexp(e_M(1,1))*-1.;
     if(fperiod>0)  SEL1OUT(a)=mfexp(e_F(1))*-1.;
     if(ftperiod>0) SEL1OUT(a)=mfexp(e_FA(1))*-1.;
     if(l2period>0) SEL1OUT(a)=LRR2(1)*-1.;
     if(l1period>0) SEL1OUT(a)=LRR1(1)*-1.;
     if(phi1period>0) SEL1OUT(a)=PHI1R(1)*-1.;
     if(phi2period>0) SEL1OUT(a)=PHI2R(1)*-1.;
     if(mselperiod>0) SEL1OUT(a)=mfexp(e_M(1,1))*-1.;

     if(mperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
     if(fperiod>0)  SEL2OUT(a)=mfexp(e_F(1))*-1.;
     if(ftperiod>0) SEL2OUT(a)=mfexp(e_FA(1))*-1.;
     if(l2period>0) SEL2OUT(a)=LRR2(1)*-1.;
     if(l1period>0) SEL2OUT(a)=LRR1(1)*-1.;
     if(phi1period>0) SEL2OUT(a)=PHI1R(1)*-1.;
     if(phi2period>0) SEL2OUT(a)=PHI2R(1)*-1.;
     if(mselperiod>0) SEL2OUT(a)=mfexp(e_M(1,1))*-1.;
   }
 }//modetype==0

FUNCTION calc_fish_surv
  if(modeltype==0){
    if(tagtype==1){
      for(a=1;a<=nages;a++){
      for(y=styr;y<=endyr;y++) {
         S(y,a)=mfexp(-1*(hsel(a)*F(styrR,y)+h(y)*(rsel(a)*FA(styrR,y))+M(y,a)));
     }
    }
   }
  if(tagtype==0){
    for(a=1;a<=nages;a++){
   for(y=styr;y<=endyr;y++) {
      S(y,a)=mfexp(-1*(hsel(a)*F(styrR,y)+M(y,a)));
   }
  }
 }
 }
 if(modeltype==1){
    if(tagtype==1){
      for(a=1;a<=nages;a++){
      for(y=styr;y<=endyr;y++) {
       age=a+y-styrR;
       if(age<=nages) S(y,a)=mfexp(-1*(hsel(age)*F(styrR,y)+h(y)*(rsel(age)*FA(styrR,y))+M(y,a)));
       if(age>nages) S(y,a)=mfexp(-1*(hsel(nages)*F(styrR,y)+h(y)*(rsel(nages)*FA(styrR,y))+M(y,a)));
     }
    }
   }
  if(tagtype==0){
    for(a=1;a<=nages;a++){
   for(y=styr;y<=endyr;y++) {
       age=a+y-styrR;
       if(age<=nages) S(y,a)=mfexp(-1*(hsel(age)*F(styrR,y)+M(y,a)));
       if(age>nages)  S(y,a)=mfexp(-1*(hsel(nages)*F(styrR,y)+M(y,a)));
   }
  }
 }
 }
FUNCTION calc_s
  // Age-indpenedent
 if(modeltype==0){
   for(a=1;a<=nages;a++){
     if(tagtype==1){
       cnt=0;
       for(t=styrR;t<=endyrR;t++){
         for(y=styr+cnt;y<=endyr;y++){
           if(t==y) s(a,t,y)=1.;
           if(t!=y) s(a,t,y)=mfexp(-(hsel(a)*F(t,y-1))-(rsel(a)*FA(t,y-1))-M(y-1,a));
         }		   
         cnt+=1;
       }
     }
   if(tagtype==0){
     cnt=0;
     for(t=styrR;t<=endyrR;t++){
       for(y=styr+cnt;y<=endyr;y++){
         if(t==y) s(a,t,y)=1.;
         if(t!=y) s(a,t,y)=mfexp(-(hsel(a)*F(t,y-1))-M(y-1,a));
       }		   
       cnt+=1;
    }
   }
 }
 }

 if(modeltype==1){
   for(a=1;a<=nages;a++){
     if(tagtype==1){
       cnt=0;
       for(t=styrR;t<=endyrR;t++){
         for(y=styr+cnt;y<=endyr;y++){
           if(t==y) s(a,t,y)=1.;
           if(t!=y){
             age=a+y-t;
             if(age<=nages) s(a,t,y)=mfexp(-(hsel(age-1)*F(t,y-1))-(rsel(age-1)*FA(t,y-1))-M(y-1,age-1));
             if(age>nages) s(a,t,y)=mfexp(-(hsel(nages)*F(t,y-1))-(rsel(nages)*FA(t,y-1))-M(y-1,nages)); 
          }
           }		   
         cnt+=1;
       }
     }
   
   if(tagtype==0){
     cnt=0;
     for(t=styrR;t<=endyrR;t++){
       for(y=styr+cnt;y<=endyr;y++){
         if(t==y) s(a,t,y)=1.;
         if(t!=y) {
          age=a+y-t;
          if(age<=nages) s(a,t,y)=mfexp(-(hsel(age-1)*F(t,y-1))-M(y-1,age-1));
          if(age>nages)  s(a,t,y)=mfexp(-(hsel(nages)*F(t,y-1))-M(y-1,nages));    
          }
       }		   
       cnt+=1;
    }
   }
 }
 }

FUNCTION calc_u_h
 if(modeltype==0){
   for(a=1;a<=nages;a++){
     cnt=0;
    if(tagtype==1){
       for(t=styrR;t<=endyrR;t++){
          for(y=styr+cnt;y<=endyr;y++){
          u_h(a,t,y)=((hsel(a)*F(t,y))/((hsel(a)*F(t,y))+(rsel(a)*FA(t,y))+M(y,a)))*(1-mfexp(-(hsel(a)*F(t,y))-(rsel(a)*FA(t,y))-M(y,a)));
        }		   
       cnt+=1;
     }
   }
   if(tagtype==0){
    for(t=styrR;t<=endyrR;t++){
      for(y=styr+cnt;y<=endyr;y++){
        u_h(a,t,y)=((hsel(a)*F(t,y))/((hsel(a)*F(t,y))+M(y,a)))*(1-mfexp(-(hsel(a)*F(t,y))-M(y,a)));
       }		   
      cnt+=1;
    }
  }
  }
 }

 if(modeltype==1){
  for(a=1;a<=nages;a++){
     cnt=0;
    if(tagtype==1){
       for(t=styrR;t<=endyrR;t++){
          for(y=styr+cnt;y<=endyr;y++){
            age=a+y-t;
            if(age<=nages) u_h(a,t,y)=((hsel(age)*F(t,y))/((hsel(age)*F(t,y))+(rsel(age)*FA(t,y))+M(y,age)))*(1-mfexp(-(hsel(age)*F(t,y))-(rsel(age)*FA(t,y))-M(y,age)));
            if(age>nages) u_h(a,t,y)=((hsel(nages)*F(t,y))/((hsel(nages)*F(t,y))+(rsel(nages)*FA(t,y))+M(y,nages)))*(1-mfexp(-(hsel(nages)*F(t,y))-(rsel(nages)*FA(t,y))-M(y,nages)));
          }   
       cnt+=1;
     }
   }
  if(tagtype==0){
     for(t=styrR;t<=endyrR;t++){
      for(y=styr+cnt;y<=endyr;y++){
        age=a+y-t;
        if(age<=nages) u_h(a,t,y)=((hsel(age)*F(t,y))/((hsel(age)*F(t,y))+M(y,age)))*(1-mfexp(-(hsel(age)*F(t,y))-M(y,age)));
        if(age>nages) u_h(a,t,y)=((hsel(nages)*F(t,y))/((hsel(nages)*F(t,y))+M(y,nages)))*(1-mfexp(-(hsel(nages)*F(t,y))-M(y,nages)));
       }		   
      cnt+=1;
    }
  }
  }
 }
FUNCTION calc_u_r
 if(modeltype==0){
 if(tagtype==1){
   for(a=1;a<=nages;a++){
     cnt=0;
     for(t=styrR;t<=endyrR;t++){
       for(y=styr+cnt;y<=endyr;y++) {
         u_r(a,t,y)=((rsel(a)*FA(t,y))/((hsel(a)*F(t,y))+(rsel(a)*FA(t,y))+M(y,a)))*(1-mfexp(-(hsel(a)*F(t,y))-(rsel(a)*FA(t,y))-M(y,a)));
        }		   
        cnt+=1;
      }
    }
  }
 }
 if(modeltype==1){
  if(tagtype==1){
   for(a=1;a<=nages;a++){
     cnt=0;
     for(t=styrR;t<=endyrR;t++){
       for(y=styr+cnt;y<=endyr;y++) {
         age=a+y-t;
         if(age<=nages) u_r(a,t,y)=((rsel(age)*FA(t,y))/((hsel(age)*F(t,y))+(rsel(age)*FA(t,y))+M(y,age)))*(1-mfexp(-(hsel(age)*F(t,y))-(rsel(age)*FA(t,y))-M(y,age)));
         if(age>nages)  u_r(a,t,y)=((rsel(nages)*FA(t,y))/((hsel(nages)*F(t,y))+(rsel(nages)*FA(t,y))+M(y,nages)))*(1-mfexp(-(hsel(nages)*F(t,y))-(rsel(nages)*FA(t,y))-M(y,nages)));
        }		   
        cnt+=1;
      }
    }
  }
 }

FUNCTION calc_s_prob
 for(a=1;a<=nages;a++){
   cnt=0;
   for(t=styrR;t<=endyrR;t++){
     looper=0;
     for(y=styr+cnt;y<=endyr;y++) {
	probs=1;
	for(kk=y-looper;kk<=y;kk++){
           probs=probs*s(a,t,kk);
          }
          s_prob(a,t,y)=probs;
          looper+=1;
      }		   
    cnt+=1;
   }
 }

FUNCTION calc_exp_prob_h
  for(a=1;a<=nages;a++){
    cnt=0;
    for(t=styrR;t<=endyrR;t++){
      dodo=0;
      for(y=styr+cnt;y<=endyr;y++){
         exp_prob_h(a,t,y)=lr1(y)*phi1rr(t)*s_prob(a,t,y)*u_h(a,t,y);
 	 dodo+=exp_prob_h(a,t,y);
       }	
       sum_prob_h(t,a)=dodo;   
       cnt+=1;
   }
 }

FUNCTION calc_exp_prob_r
 if(tagtype==1){
    for(a=1;a<=nages;a++){
      cnt=0;
      for(t=styrR;t<=endyrR;t++){
        dodo=0;
        for(y=styr+cnt;y<=endyr;y++){
          exp_prob_r(a,t,y)=lr2(y)*phi2rr(t)*s_prob(a,t,y)*u_r(a,t,y);
	  dodo+=exp_prob_r(a,t,y);
        }	
        sum_prob_r(t,a)=dodo;   
        cnt+=1;
     }
   }
 }

FUNCTION calc_LL

  for(a=1;a<=nages;a++){
    cnt=0;
   for(t=styrR;t<=endyrR;t++) {
     for(y=styr+cnt;y<=endyr;y++){
       ll_h(a,t,y)=0.;
       ll_r(a,t,y)=0.;
       if(rh(a,t,y)>0.) ll_h(a,t,y)=rh(a,t,y)*log(exp_prob_h(a,t,y));  
       if(tagtype==1){
         if(rr(a,t,y)>0.) ll_r(a,t,y)=rr(a,t,y)*log(exp_prob_r(a,t,y));  
       }
    }		   
    cnt+=1;
   }
 }
 if(tagtype==0){
   for(a=1;a<=nages;a++){
    for (t=styrR;t<=endyrR;t++){
     ll_ns(t,a)=0.;
     if(N(t,a)>0.) ll_ns(t,a)=(N(t,a)-tags(t,a))*log(1-(sum_prob_h(t,a)));
    }
  }
 }
 if(tagtype==1){
  for(a=1;a<=nages;a++){
   for (t=styrR;t<=endyrR;t++){
      ll_ns(t,a)=0.;
      if(N(t,a)>0.) ll_ns(t,a)=(N(t,a)-tags(t,a))*log(1-(sum_prob_h(t,a)+sum_prob_r(t,a)));
     }
   }
 }


FUNCTION evaluate_the_objective_function
 f_tag=0;
 for(a=1;a<=nages;a++){
   cnt=0;
   for(t=styrR;t<=endyrR;t++){
    for(y=styr+cnt;y<=endyr;y++) {
      if(tagtype==1) f_tag+=ll_h(a,t,y)+ll_r(a,t,y);
      if(tagtype==0) f_tag+=ll_h(a,t,y);
      }		 
    cnt+=1;
   }
 }
  for(a=1;a<=nages;a++){
   for(t=styrR;t<=endyrR;t++){
      f_tag+=ll_ns(t,a);
    }
  }
  if(use_constant>0) f=-1.0*(f_tag+constant);
  if(use_constant<=0) f=-1.0*f_tag;

 
REPORT_SECTION
 //Count cells for degrees of freedom
  up_count=0;
  for(a=1;a<=nages;a++){
    cnt=0;
    for(t=styrR;t<=endyrR;t++){
      for(y=styr+cnt;y<=endyr;y++){
         up_count+=1; 
        if(tagtype==1){
         up_count+=1; 
        }
      }
      cnt+=1;
  }
 }
 
  for(a=1;a<=nages;a++){
   cnt=0;
   for (t=styrR;t<=endyrR;t++){
     for (y=styr+cnt;y<=endyr;y++){
        exp_r_h(a,t,y)=exp_prob_h(a,t,y)*N(t,a);
       if(tagtype==1) exp_r_r(a,t,y)=exp_prob_r(a,t,y)*N(t,a);
      }	
    cnt+=1;
   }
 }

 for(a=1;a<=nages;a++){
  cnt=0;
  for(t=styrR;t<=endyrR;t++) {
    for(y=styr+cnt;y<=endyr;y++) {
         chi_h(a,t,y)=0.;
         pear_h(a,t,y)=0.;
         stdres_h(a,t,y)=0.;
        if(exp_r_h(a,t,y)>0.){
         chi_h(a,t,y)=square(rh(a,t,y)-exp_r_h(a,t,y))/exp_r_h(a,t,y);
         pear_h(a,t,y)=(rh(a,t,y)-exp_r_h(a,t,y))/sqrt(exp_r_h(a,t,y));
         stdres_h(a,t,y)=(rh(a,t,y)-exp_r_h(a,t,y))/sqrt(exp_r_h(a,t,y)*(1.-exp_r_h(a,t,y)/N(t,a)));
       }
         if(tagtype==1){
           chi_r(a,t,y)=0.;
           pear_r(a,t,y)=0.;
           stdres_r(a,t,y)=0.;
          if(exp_r_r(a,t,y)>0.){
           chi_r(a,t,y)=square(rr(a,t,y)-exp_r_r(a,t,y))/exp_r_r(a,t,y);
           pear_r(a,t,y)=(rr(a,t,y)-exp_r_r(a,t,y))/sqrt(exp_r_r(a,t,y));
           stdres_r(a,t,y)=(rr(a,t,y)-exp_r_r(a,t,y))/sqrt(exp_r_r(a,t,y)*(1.-exp_r_r(a,t,y)/N(t,a)));
         }
        }
     }	
      cnt+=1;
   }
  }
    for(a=1;a<=nages;a++){
      for (t=styrR;t<=endyrR;t++) {
         if(tagtype==1) exp_ns(t,a)=N(t,a)*(1-(sum_prob_h(t,a)+sum_prob_r(t,a)));
          if(tagtype==0) exp_ns(t,a)=N(t,a)*(1-(sum_prob_h(t,a)));
      }
    }
  //Not seen chi
 for(a=1;a<=nages;a++){
  for (t=styrR;t<=endyrR;t++) {
        chi_ns(t,a)=0.;
        pear_ns(t,a)=0.;
        stdres_ns(t,a)=0.;
       if(exp_ns(t,a)>0.){
        ns_count+=1;
        chi_ns(t,a)=square((N(t,a)-tags(t,a))-exp_ns(t,a))/exp_ns(t,a);
        pear_ns(t,a)=((N(t,a)-tags(t,a))-exp_ns(t,a))/sqrt(exp_ns(t,a));
        stdres_ns(t,a)=((N(t,a)-tags(t,a))-exp_ns(t,a))/sqrt(exp_ns(t,a)*(1.-exp_ns(t,a)/N(t,a)));
       }
   }
 }
 //total chi square
 if(tagtype==1) up_chi=sum(chi_r)+sum(chi_h)+sum(chi_ns);
 if(tagtype==0) up_chi=sum(chi_h)+sum(chi_ns);
  K=0; 
  
  if(fperiod>0) K=K+fp;
  if(ftperiod>0) K=K+fap;
  if(l1period>0) K=K+lp1;
  if(l2period>0) K=K+lp2;
  if(phi1period>0) K=K+phi1p;
  if(phi2period>0) K=K+phi2p;
  if(nonmixh>0) K=K+nmph;
  if(nonmixr>0) K=K+nmpr;
  if(sel1period>0) K=K+nsel1;
  if(sel2period>0) K=K+nsel2;
  if(mperiod>0 && mselperiod>0) K=K+mp*mselp;
  if(mperiod>0 && mselperiod<=0) K=K+mp;
  if(mperiod<=0 && mselperiod>0) K=K+mselp;
  up_df=up_count+ns_count-K;
  up_chat=up_chi/up_df;
 
 //calc_pooled_cells
 // Pool harvested cells
  for(a=1;a<=nages;a++){
    cnt=0;
    for(t=styrR;t<=endyrR;t++) { 
        for(y=styr+cnt;y<=endyr;y++){
            pool_h_e(a,t,y)=0.;
            pool_h(a,t,y)=0.;
            pool_h_e(a,t,y)=exp_r_h(a,t,y);
            if(rh(a,t,y)>=0.) pool_h(a,t,y)=rh(a,t,y); 
        }
      cnt+=1;
    }
  }
 hless=0;
 for(a=1;a<=nages;a++){
   cnt=0;
   for(t=styrR;t<=endyrR;t++) {    
     for(y=endyr;y>=styr+cnt;y--){
          if(pool_h_e(a,t,y)>=2.) 
            {
             pool_h(a,t,y)=pool_h(a,t,y);
             pool_h_e(a,t,y)=pool_h_e(a,t,y);
            }
          if(pool_h_e(a,t,y)>=0. && pool_h_e(a,t,y)<2.){
             if (y!=styr+cnt){
               hless+=1;
               pool_h_e(a,t,y-1)=pool_h_e(a,t,y-1)+pool_h_e(a,t,y);
                pool_h(a,t,y-1)=pool_h(a,t,y-1)+pool_h(a,t,y);
                pool_h(a,t,y)=0.;
                pool_h_e(a,t,y)=0.;
               }
               if (y==styr+cnt) break;
            }
         }//for
         cnt+=1;
     }//for
  }
// Pool released cells
 if(tagtype==1){
   for(a=1;a<=nages;a++){
     cnt=0;
     for(t=styrR;t<=endyrR;t++) { 
        for(y=styr+cnt;y<=endyr;y++) {
            pool_r_e(a,t,y)=0.;
            pool_r(a,t,y)=0.;
            pool_r_e(a,t,y)=exp_r_r(a,t,y);
            if(rr(a,t,y)>=0.) pool_r(a,t,y)=rr(a,t,y);
        }
      cnt+=1;
   }
 }
  rless=0;
  for(a=1;a<=nages;a++){
     cnt=0;
     for(t=styrR;t<=endyrR;t++){    
      for(y=endyr;y>=styr+cnt;y--){
          if(pool_r_e(a,t,y)>=2.){
             pool_r(a,t,y)=pool_r(a,t,y);
             pool_r_e(a,t,y)=pool_r_e(a,t,y);
            }
          if(pool_r_e(a,t,y)>=0. && pool_r_e(a,t,y)<2.){
              if (y!=styr+cnt){
                rless+=1;
                pool_r_e(a,t,y-1)=pool_r_e(a,t,y-1)+pool_r_e(a,t,y);
                pool_r(a,t,y-1)=pool_r(a,t,y-1)+pool_r(a,t,y);
                pool_r(a,t,y)=0.;
                pool_r_e(a,t,y)=0.;
               }
               if (y==styr+cnt) break;
            }
         }//for
         cnt+=1;
     }//for
  }
 }//Tagtype==1
 //Calculate pooled degrees of freedom
  p_df=up_count+ns_count-rless-hless-K;
 
  //Pooled Chi-square
  for(a=1;a<=nages;a++){
     cnt=0;
     for (t=styrR;t<=endyrR;t++) {
       for (y=styr+cnt;y<=endyr;y++) {
        p_chi_h(a,t,y)=0.;
        p_chi_r(a,t,y)=0.;
        if(pool_h_e(a,t,y)>0.){
          p_chi_h(a,t,y)=square(pool_h(a,t,y)-pool_h_e(a,t,y))/pool_h_e(a,t,y);
         }
       if(tagtype==1){
        if(pool_r_e(a,t,y)>0.){
          p_chi_r(a,t,y)=square(pool_r(a,t,y)-pool_r_e(a,t,y))/pool_r_e(a,t,y);
         }
       }
      }	
      cnt+=1;
   }
 }

  if(tagtype==1) p_chi=sum(p_chi_h)+sum(p_chi_r)+sum(chi_ns);
  if(tagtype==0) p_chi=sum(p_chi_h)+sum(chi_ns);
  p_chat=p_chi/p_df;

 report<<"# Instantaneous Rates Tag Return Models (IRATE) Version 1.0"<<endl;
 report<<"# Start time for run: "<<ctime(&start)<<"";
 if (modeltype==0) report<<"# Age-Independent Model"<<endl;
 if (modeltype==1)  report<<"# Age-Dependent Model"<<endl;
 if (tagtype==0)  report<<"# Harvest Only Tag Returns"<<endl;
 if (tagtype==1) report<<"# Harvest and Catch/Release Tag Returns"<<endl;
 report<<"# Log-L"<<endl;
 report<<-1*f<<endl;
 report<<"# K"<<endl;
 report<<K<<endl;
 report<<"# AIC"<<endl;
 if(use_constant>0){
 report<<-1.*2*(-1.0*f)+2*K<<endl;
 report<<"# AICc"<<endl;
 report<<-1.*2*(-1.0*f)+2*K+(2*K*(K+1))/(sum(N)-K-1)<<endl;
 }
 if(use_constant<=0){
 report<<-1.*2*f_tag+2*K<<endl;
 report<<"# AICc"<<endl;
 report<<-1.*2*f_tag+2*K+(2*K*(K+1))/(sum(N)-K-1)<<endl;
 }
 report<<"# Effective Sample Size"<<endl;
 report<<sum(N)<<endl;
 report<<"# Unpooled Chi-square"<<endl;
 report<<up_chi<<endl;
 report<<"# Upooled df"<<endl;
 report<<up_df<<endl;
 report<<"# Unpooled c-hat"<<endl;
 report<<up_chat<<endl;
 report<<"# Pooled Chi-square"<<endl;
 report<<p_chi<<endl;
 report<<"# Pooled df"<<endl;
 report<<p_df<<endl;
 report<<"# Pooled c-hat"<<endl;
 report<<p_chat<<endl;
 for(a=1;a<=nages;a++){
  if (modeltype==0) report<<"# Harvest Residuals"<<endl;
  if (modeltype==1) report<<"# Harvest Residuals of Age Class "<<a<<endl;
   for(y=styrR;y<=endyrR;y++){
     for(t=styr;t<=endyr;t++){
       if(t<endyr) report<<stdres_h(a,y,t)<<" ";
       if(t==endyr) report<<stdres_h(a,y,t)<<endl;
     }
   }
  }

 for(a=1;a<=nages;a++){
    if (modeltype==0) report<<"# Release Residuals"<<endl;
    if (modeltype==1) report<<"# Release Residuals of Age Class "<<a<<endl;
    for(y=styrR;y<=endyrR;y++){
      for(t=styr;t<=endyr;t++){
         if(t<endyr) report<<stdres_r(a,y,t)<<" ";
         if(t==endyr) report<<stdres_r(a,y,t)<<endl;
      }
   }
 }

   for(a=1;a<=nages;a++){
     if (modeltype==0) report<<"# Not seen Residuals"<<endl;
     if (modeltype==1) report<<"# Not seen Residuals of Age Class "<<a<<endl;
     for(y=styrR;y<=endyrR;y++){
        if(y<endyrR) report<<stdres_ns(y,a)<<" ";
        if(y==endyrR) report<<stdres_ns(y,a)<<endl;
     }
   }
 report<<"#*************************Observed and Calculated Data***************************************"<<endl;
 report << "# Obs Recoveries of harvest fish "<< endl;
 report<<rh<<endl;
 report<<" "<<endl;
 report << "# Obs Recoveries of release fish "<< endl;
 report<<rr<<endl;
 report<<" "<<endl;
 report << "# Total Released "<< endl;
 report<<N<<endl;
 report<<" "<<endl;
 report <<"# Total Recovered Tags"<<endl;
 report <<tags<<endl;
 report<<" "<<endl;
 report << "# s matrix" << endl;
 report <<s<<endl;
 report<<" "<<endl;
 report << "# S_prob matrix" << endl;
 report <<s_prob<<endl;
 report<<" "<<endl;
 report << "# Exploitation Rate of harvested fish" << endl;
 report <<u_h<<endl;
 report<<" "<<endl;
 if(tagtype==1){
    report << "# Exploitation Rate of released fish" << endl;
    report <<u_r<<endl;
    report<<" "<<endl;
  }
 report <<"# Expected Probability of harvested fish"<<endl;
 report<<exp_prob_h<<endl;
 report<<" "<<endl;
  if(tagtype==1){
    report <<"# Expected Probability of released fish"<<endl;
    report<<exp_prob_r<<endl;
    report<<" "<<endl;
   }
 report<<"# Not Seen Probability"<<endl;
 if(tagtype==1) report<<1-(sum_prob_h+sum_prob_r)<<endl;
 if(tagtype==0) report<<1-(sum_prob_h)<<endl;
 report<<" "<<endl;
 report <<"# Expected Number of harvested fish"<<endl;
 report<<exp_r_h<<endl;
 report<<" "<<endl;
 if(tagtype==1){
   report <<"# Expected Number of released fish"<<endl;
   report<<exp_r_r<<endl;
   report<<" "<<endl;
  }
 report <<"# Expected Number of not seen"<<endl;
 report<<exp_ns<<endl;
 report<<" "<<endl;
 report <<"# Cell Likelihoods of harvested fish"<<endl;
 report<<ll_h<<endl;
 report<<" "<<endl;
 if(tagtype==1){
   report <<"# Cell Likelihoods of released fish"<<endl;
   report<<ll_r<<endl;
   report<<" "<<endl;
 }
 report <<"# Cell Likelihoods of not seen"<<endl;
 report<<ll_ns<<endl;
 report<<" "<<endl;
 report <<"# Unpooled Chi-squares of Harvested Fish"<<endl;
 report<<chi_h<<endl;
 report<<" "<<endl;
 if(tagtype==1){
   report <<"# Unpooled Chi-squares of Released Fish"<<endl;
   report<<chi_r<<endl;
   report<<" "<<endl;
 }
 report <<"# Chi-squares of Not Seen"<<endl;
 report<<chi_ns<<endl;
 report<<" "<<endl;
 report <<"# Pooled Cells of Harvested Fish"<<endl;
 report<<pool_h<<endl;
 report<<" "<<endl;
 report <<"# Pooled Expected Cells of Harvested Fish"<<endl;
 report<<pool_h_e<<endl;
 report<<" "<<endl;
 if(tagtype==1){
   report <<"# Pooled Cells of Released Fish"<<endl;
   report<<pool_r<<endl;
   report<<" "<<endl;
 }
 if(tagtype==1){
   report <<"# Pooled Expected Cells of Released Fish"<<endl;
   report<<pool_r_e<<endl;
   report<<" "<<endl;
 }
 report <<"# Pooled Chi-squares of Harvested Fish"<<endl;
 report<<p_chi_h<<endl;
 report<<" "<<endl;
 if(tagtype==1){
    report <<"# Pooled Chi-squares of Released Fish"<<endl;
   report<<p_chi_r<<endl;
   report<<" "<<endl;
 }
 report <<"# Pearson Residuals for harvested fish"<<endl;
 report<<pear_h<<endl;
 report<<" "<<endl;
 if(tagtype==1){
   report <<"# Pearson Residuals for released fish"<<endl;
   report<<pear_r<<endl;
   report<<" "<<endl;
 }

 report <<"# Pearson Residuals for not seen"<<endl;
 report<<pear_ns<<endl;
 report<< " "<<endl;
FINAL_SECTION
 time(&finish);
 elapsed_time = difftime(finish,start);
  hour = long(elapsed_time)/3600;
  minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout<<endl<<" starting time: "<<ctime(&start);
 cout<<" finishing time: "<<ctime(&finish);
 cout<<" This run took: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds."<<endl<<endl<<endl;

