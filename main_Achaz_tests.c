//
//  main.c
//  Achaz_stats
//
//  Created by Sebastian E. Ramos Onsins on 02/05/2016.
//  Copyright © 2016 Sebastian E. Ramos Onsins. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define ACHAZSTATS "Achaz_tests: version 20160922.\n"
#define SAMPLE_LARGE 4000
#define NITER 1e6
#define BOOTSTRAP1 1

struct stats {
    /*standard test*/
    double Dtaj;
    double Dfl;
    double Ffl;
    double Hnfw;
    double Ez;
    double Yach;

    /*standard test without singletons*/
    double Dtaj0;
    double Dfl0;
    double Ffl0;
    double Hnfw0;
    double Ez0;
}statistics;

struct thetas {
    /*standard test*/
    double theta1_Dtaj;
    double theta2_Dtaj;
    double theta1_Dfl;
    double theta2_Dfl;
    double theta1_Ffl;
    double theta2_Ffl;
    double theta1_Hnfw;
    double theta2_Hnfw;
    double theta1_Ez;
    double theta2_Ez;
    double theta1_Yach;
    double theta2_Yach;
    
    /*standard test without singletons*/
    double theta1_Dtaj0;
    double theta2_Dtaj0;
    double theta1_Dfl0;
    double theta2_Dfl0;
    double theta1_Ffl0;
    double theta2_Ffl0;
    double theta1_Hnfw0;
    double theta2_Hnfw0;
    double theta1_Ez0;
    double theta2_Ez0;
    double theta1_Yach0;
    double theta2_Yach0;
}Th;

double freqtesto_achaz(long int sample_size,long int *fr,long int singleton,double *w1,double *w2,
                       double *theta1,double *theta2);
double freqtestn_achaz(long int sample_size,long int *fr,long int singleton,double *w1,double *w2,
                       double *theta1,double *theta2);

int main(int argc, const char * argv[]) {
    int arg;
    void usage(void);
    int outgroup_presence=0;
    long int nsam = 0;
    long int nfreq = 0;
    unsigned long *freq;
    unsigned long segsites=0;
 
    long int *sfreq;
    long int *snfreq;
    double *w1;
    double *w2;
    long int *sfreqn;
    
    long int x,fr,sfr;
    char *sf;
    long int seed=0;
    void init_seed1(long int);
    
    printf(ACHAZSTATS);

    /*read stdin and include data in variables and vectors: nsam, outg, freqs_i*/
    if((freq = (unsigned long *) calloc( (unsigned long)1, sizeof(unsigned long) )) == 0) {
        printf("Error allocating memory");
        exit(1);
    }

    if( argc > 1 )
    {
        arg = 1;
        while(arg < argc)
        {
            if( argv[arg][0] != '-' )
            {
                if(argv[arg][0] == '>')
                    break;
                printf(" argument should be -%s ?\n", argv[arg]);
                usage();
                exit(1);
            }
            
            switch (argv[arg][1])
            {
                case 'f' : /*outgroup_presence (1/0) and nsam followed by nsam-1 or (nsam-1)/2N numbers indicating the frequency spectra*/
                    arg++;
                    outgroup_presence = (int)atoi(argv[arg]);
                    arg++;
                    nsam = atoi(argv[arg]);
                    
                    if(outgroup_presence == 1)
                        nfreq = nsam - 1;
                    if(outgroup_presence == 0) {
                        if(nsam/2.0 == nsam/2)
                            nfreq = nsam/2;
                        else
                            nfreq = (nsam-1)/2;
                    }
                    
                    if((freq = (unsigned long *) calloc( (unsigned long)nfreq, sizeof(unsigned long) )) == 0) {
                        printf("Error allocating memory");
                        exit(1);
                    }
                    segsites = 0;
                    for(x=0;x<nfreq;x+=1)
                    {
                        arg++;
                        freq[x] = atol(argv[arg]);
                        segsites += freq[x];
                    }
                    break;
                case 's' : /*outgroup_presence (1/0) and nsam followed by indertemined number of freqs,':' plus the value of the freq together*/
                    arg++;
                    outgroup_presence = (int)atoi(argv[arg]);
                    arg++;
                    nsam = atoi(argv[arg]);
                    sf = (char *)calloc(100,sizeof(char));
                    arg++;
                    
                    if(outgroup_presence == 1)
                        nfreq = nsam - 1;
                    if(outgroup_presence == 0) {
                        if(nsam/2.0 == nsam/2)
                            nfreq = nsam/2;
                        else
                            nfreq = (nsam-1)/2;
                    }
                    
                    if((freq = (unsigned long *) calloc( (unsigned long)nfreq, sizeof(unsigned long) )) == 0) {
                        printf("Error allocating memory");
                        exit(1);
                    }
                    segsites = 0;
                    while((arg < argc) && (argv[arg][0] != '>')) {
                        fr=0;
                        while((arg < argc) && (argv[arg][fr] != '>') &&
                              (argv[arg][fr] >= '0') && (argv[arg][fr] <= '9') &&
                              (argv[arg][fr] != ':')) {
                            sf[fr]=argv[arg][fr];
                            fr++;
                        }
                        sf[fr] = '\0';
                        x=atol(sf);
                        fr++;
                        sfr=fr;
                        while((arg < argc) && (argv[arg][fr] != '>') &&
                              (argv[arg][fr] >= '0') && (argv[arg][fr] <= '9') &&
                              (argv[arg][fr] != ' ')) {
                            sf[fr-sfr]=argv[arg][fr];
                            fr++;
                        }
                        sf[fr-sfr]= '\0';
                        freq[x-1] = atol(sf);
                        segsites += freq[x-1];
                               
                        arg++;
                    }
                    free(sf);
                    break;
                case 'd' : /*seed in case sample_size is bigger than SAMPLE_LARGE*/
                    arg++;
                    seed = atol(argv[arg]);
                    break;
                case 'h' : /* h HEEEEEEEEEEELPPPPPPPPPPPPPP!!! */
                    usage();
                    exit(0);
                    break;
                default :
                    usage();
                    exit(0);
                    break;
           }
            arg++;
        }
        
        if(nsam > 1) {
            /*calculate all statistics*/
            sfreq  = (long int *)calloc(nsam+1,sizeof(long int));
            snfreq = (long int *)calloc(nsam+1,sizeof(long int));
            
            w1 		= (double *) calloc(nsam,sizeof(double));
            w2 		= (double *) calloc(nsam,sizeof(double));
            sfreqn 	= (long int *)calloc(nsam,sizeof(long int));
            
            for(x=0;x<nfreq;x++) {
                sfreq[x+1] = freq[x];
            }
            
            if(seed == 0)
                seed = 123456;
            /*OUTGROUP*/
            init_seed1(seed);
            if(outgroup_presence == 1) {
                /*Tajima's D with outgroup....*/
                for(x=1;x<nsam;x++) {
                    w1[x] = (double)(nsam - x);
                    w2[x] = 1.0/(double)x;
                }
                statistics.Dtaj  = freqtesto_achaz(nsam,sfreq,1,w1,w2,&Th.theta1_Dtaj,&Th.theta2_Dtaj);
                for(x=1;x<nsam;x++) {
                    if(x==1) w1[x] = 0.0;
                    else w1[x] = (double)(nsam - x);
                    if(x==1) w2[x] = 0.0;
                    else w2[x] = 1.0/(double)x;
                }
                statistics.Dtaj0 = freqtesto_achaz(nsam,sfreq,0,w1,w2,&Th.theta1_Dtaj0,&Th.theta2_Dtaj0);
                
                /*Fu and Li's D*/
                for(x=1;x<nsam;x++) {
                    w1[x] = 1.0/(double)x;
                    if(x==1) w2[x] = 1.0;
                    else w2[x] = 0.0;
                }
                statistics.Dfl  = freqtesto_achaz(nsam,sfreq,1,w1,w2,&Th.theta1_Dfl,&Th.theta2_Dfl);
                
                /*Fu and Li's F*/
                for(x=1;x<nsam;x++) {
                    w1[x] = (double)(nsam - x);
                    if(x==1) w2[x] = 1.0;
                    else w2[x] = 0.0;
                }
                statistics.Ffl  = freqtesto_achaz(nsam,sfreq,1,w1,w2,&Th.theta1_Ffl,&Th.theta2_Ffl);
                
                /*Fay and Wu's H*/
                for(x=1;x<nsam;x++) {
                    w1[x] = (double)(nsam - x);
                    w2[x] = (double)x;
                }
                statistics.Hnfw  = freqtesto_achaz(nsam,sfreq,1,w1,w2,&Th.theta1_Hnfw,&Th.theta2_Hnfw);
                for(x=1;x<nsam;x++) {
                    if(x==1) w1[x] = 0.0;
                    else w1[x] = (double)(nsam - x);
                    if(x==1) w2[x] = 0.0;
                    else w2[x] = (double)x;
                }
                statistics.Hnfw0 = freqtesto_achaz(nsam,sfreq,0,w1,w2,&Th.theta1_Hnfw0,&Th.theta2_Hnfw0);
                
                /*Zengs' E*/
                for(x=1;x<nsam;x++) {
                    w1[x] = (double)1;
                    w2[x] = 1.0/(double)x;
                }
                statistics.Ez  = freqtesto_achaz(nsam,sfreq,1,w1,w2,&Th.theta1_Ez,&Th.theta2_Ez);
                for(x=1;x<nsam;x++) {
                    if(x==1) w1[x] = 0.0;
                    else w1[x] = (double)1;
                    if(x==1) w2[x] = 0.0;
                    else w2[x] = 1.0/(double)x;
                }
                statistics.Ez0 = freqtesto_achaz(nsam,sfreq,0,w1,w2,&Th.theta1_Ez0,&Th.theta2_Ez0);
     
                /*Achaz's Y*/
                for(x=1;x<nsam;x++) {
                    if(x==1) w1[x] = 0.0;
                    else w1[x] = (double)(nsam - x);
                    if(x==1) w2[x] = 0.0;
                    else w2[x] = 1.0/(double)x;
                }
                statistics.Yach = freqtesto_achaz(nsam,sfreq,0,w1,w2,&Th.theta1_Yach,&Th.theta2_Yach);

                /*write statistics with the names of each*/
                printf("\n");
                printf("TESTS USING UNFOLDED SFS:\n");
                /*
                if(nsam > SAMPLE_LARGE) {
                    printf("**************************************************************************************\n");
                    printf("WARNING: for sample sizes > %ld THE VARIANCES OF TESTS ARE APPROXIMATED\n",SAMPLE_LARGE);
                   printf("**************************************************************************************\n");
                }
                */
                printf("DATA: nsam= %ld S=%ld SFS (unfolded):",nsam,segsites);
                for(x=1;x<nsam;x++) {
                    if(sfreq[x]>0) printf(" f[%ld]= %ld",x,sfreq[x]);
                }
                printf("\n");
                if(nsam > SAMPLE_LARGE) printf("seed= %ld\n",seed);
                if(nsam > 1e7 || segsites == 0 || nsam < 3) {
                    /*printf("TAJIMA_D: thetaT= %.5f thetaW= %.5f TEST= %.5f \n",Th.theta1_Dtaj,Th.theta2_Dtaj,statistics.Dtaj);*/
                    printf("FU&LI_D: thetaW= %.5f thetaFL= %.5f TEST= NA \n",Th.theta1_Dfl,Th.theta2_Dfl);
                    printf("FU&LI_F: thetaT= %.5f thetaFL= %.5f TEST= NA \n",Th.theta1_Ffl,Th.theta2_Ffl);
                    printf("FAY&WU_H: thetaT= %.5f thetaH= %.5f TEST= NA \n",Th.theta1_Hnfw,Th.theta2_Hnfw);
                    printf("ZENG_E: thetaE= %.5f thetaW= %.5f TEST= NA \n",Th.theta1_Ez,Th.theta2_Ez);
                    if(nsam < 3) {
                        printf("ACHAZ_Y: thetaTns= NA thetaWns= NA TEST= NA \n");
                        /*printf("TAJIMA_D_NO_SINGL: thetaTns= NA thetaWns= NA TEST= NA \n");*/
                        printf("FAY&WU_H_NO_SINGL: thetaTns= NA thetaHns= NA TEST= NA \n");
                        printf("ZENG_E_NO_SINGL: thetaEns= NA thetaW= NA TEST= NA \n");
                    }
                    else {
                        printf("ACHAZ_Y: thetaTns= %.5f thetaWns= %.5f TEST= NA \n",Th.theta1_Yach,Th.theta2_Yach);
                        /*printf("TAJIMA_D_NO_SINGL: thetaTns= %.5f thetaWns= %.5f TEST= NA \n",Th.theta1_Dtaj0,);*/
                        printf("FAY&WU_H_NO_SINGL: thetaTns= %.5f thetaHns= %.5f TEST= NA \n",Th.theta1_Hnfw0,Th.theta2_Hnfw0);
                        printf("ZENG_E_NO_SINGL: thetaEns= %.5f thetaW= %.5f TEST= NA \n",Th.theta1_Ez0,Th.theta2_Ez0);
                    }
                }
                else {
                    /*printf("TAJIMA_D: thetaT= %.5f thetaW= %.5f TEST= %.5f \n",Th.theta1_Dtaj,Th.theta2_Dtaj,statistics.Dtaj);*/
                    printf("FU&LI_D: thetaW= %.5f thetaFL= %.5f TEST= %.5f \n",Th.theta1_Dfl,Th.theta2_Dfl,statistics.Dfl);
                    printf("FU&LI_F: thetaT= %.5f thetaFL= %.5f TEST= %.5f \n",Th.theta1_Ffl,Th.theta2_Ffl,statistics.Ffl);
                    printf("FAY&WU_H: thetaT= %.5f thetaH= %.5f TEST= %.5f \n",Th.theta1_Hnfw,Th.theta2_Hnfw,statistics.Hnfw);
                    printf("ZENG_E: thetaE= %.5f thetaW= %.5f TEST= %.5f \n",Th.theta1_Ez,Th.theta2_Ez,statistics.Ez);
                    if(statistics.Yach == -10000)
                        printf("ACHAZ_Y: thetaTns= %.5f thetaWns= %.5f TEST= NA \n",Th.theta1_Yach,Th.theta2_Yach);
                    else
                        printf("ACHAZ_Y: thetaTns= %.5f thetaWns= %.5f TEST= %.5f \n",Th.theta1_Yach,Th.theta2_Yach,statistics.Yach);
                    
                    /*printf("TAJIMA_D_NO_SINGL: thetaTns= %.5f thetaWns= %.5f TEST= %.5f \n",Th.theta1_Dtaj0,Th.theta2_Dtaj0,statistics.Dtaj0);*/
                    if(statistics.Hnfw0 == -10000)
                        printf("FAY&WU_H_NO_SINGL: thetaTns= %.5f thetaHns= %.5f TEST= NA \n",Th.theta1_Hnfw0,Th.theta2_Hnfw0);
                    else
                        printf("FAY&WU_H_NO_SINGL: thetaTns= %.5f thetaHns= %.5f TEST= %.5f \n",Th.theta1_Hnfw0,Th.theta2_Hnfw0,statistics.Hnfw0);
                    if(statistics.Ez0 == -10000)
                        printf("ZENG_E_NO_SINGL: thetaEns= %.5f thetaW= %.5f TEST= NA \n",Th.theta1_Ez0,Th.theta1_Ez0);
                    else
                        printf("ZENG_E_NO_SINGL: thetaEns= %.5f thetaW= %.5f TEST= %.5f \n",Th.theta1_Ez0,Th.theta2_Ez0,statistics.Ez0);
                }
            }
            
            /*NO OUTGROUP*/
            init_seed1(seed);
            if(outgroup_presence == 1) {
                if((double)nsam/(double)2 == (int)nsam/(int)2)
                    nfreq = nsam/2;
                else
                    nfreq = (nsam-1)/2;
            }
            for(x=1;x<=floor(nsam/2);x++) {
                if(x == nsam-x) sfreqn[x] = (int)((sfreq[x] + sfreq[nsam-x])/2.0);
                else sfreqn[x] = sfreq[x] + sfreq[nsam-x];
            }
            /*check analysis using Achaz approach*/
            /*Tajima's D*/
            for(x=1;x<=floor(nsam/2);x++) {
                if(x == nsam-x) w1[x] = (double)nsam/2.0;
                else w1[x] = (double)nsam/1.0;
                if(x == nsam-x) w2[x] = (double)nsam/((double)x*(double)(nsam - x)*2.0);
                else w2[x] = (double)nsam/((double)x*(double)(nsam - x)*1.0);
            }
            statistics.Dtaj  = freqtestn_achaz(nsam,sfreqn,1,w1,w2,&Th.theta1_Dtaj,&Th.theta2_Dtaj);
            for(x=1;x<=floor(nsam/2);x++) {
                if(x == 1) w1[x] = 0.0;
                else {
                    if(x == nsam-x) w1[x] = (double)nsam/2.0;
                    else w1[x] = (double)nsam/1.0;
                }
                if(x == 1) w2[x] = 0.0;
                else {
                    if(x == nsam-x) w2[x] = (double)nsam/((double)x*(double)(nsam - x)*2.0);
                    else w2[x] = (double)nsam/((double)x*(double)(nsam - x)*1.0);
                }
            }
            statistics.Dtaj0 = freqtestn_achaz(nsam,sfreqn,0,w1,w2,&Th.theta1_Dtaj0,&Th.theta2_Dtaj0);

            /*Fu and Li's D* */
            for(x=1;x<=floor(nsam/2);x++) {
                if(x == nsam-x) w1[x] = (double)nsam/((double)x*(double)(nsam - x)*2.0);
                else w1[x] = (double)nsam/((double)x*(double)(nsam - x)*1.0);
                if(x==1) w2[x] = nsam;
                else w2[x] = 0.0;
            }
            statistics.Dfl  = freqtestn_achaz(nsam,sfreqn,1,w1,w2,&Th.theta1_Dfl,&Th.theta2_Dfl);

            /*Fu and Li's F* */
            for(x=1;x<=floor(nsam/2);x++) {
                if(x == nsam-x) w1[x] = (double)nsam/2.0;
                else w1[x] = (double)nsam/1.0;
                if(x==1) w2[x] = nsam;
                else w2[x] = 0.0;
            }
            statistics.Ffl  = freqtestn_achaz(nsam,sfreqn,1,w1,w2,&Th.theta1_Ffl,&Th.theta2_Ffl);


            /*Achaz's Y* */
            for(x=1;x<=floor(nsam/2);x++) {
                if(x == 1) w1[x] = 0.0;
                else {
                    if(x == nsam-x) w1[x] = (double)nsam/2.0;
                    else w1[x] = (double)nsam/1.0;
                }
                if(x == 1) w2[x] = 0.0;
                else {
                    if(x == nsam-x) w2[x] = (double)nsam/((double)x*(double)(nsam - x)*2.0);
                    else w2[x] = (double)nsam/((double)x*(double)(nsam - x)*1.0);
                }
            }
            statistics.Yach = freqtestn_achaz(nsam,sfreqn,0,w1,w2,&Th.theta1_Yach,&Th.theta2_Yach);

            /*write statistics with the names of each*/
            printf("\n");
            printf("TESTS USING FOLDED SFS:\n");
            /*
            if(nsam > SAMPLE_LARGE) {
                printf("*************************************************************************************\n");
                printf("WARNING: for sample sizes > %ld THE VARIANCES OF THE TEST ARE APPROXIMATED\n",SAMPLE_LARGE);
                printf("*************************************************************************************\n");
            }
            */
            printf("DATA: nsam= %ld S=%ld SFS (folded):",nsam,segsites);
            for(x=1;x<=nfreq;x++) {
                if(sfreq[x]>0) printf(" f[%ld]= %ld",x,sfreqn[x]);
            }
            printf("\n");
            if(nsam > SAMPLE_LARGE) printf("seed= %ld\n",seed);
            if(nsam > 1e7 || segsites == 0 || nsam < 4) {
                printf("TAJIMA_D: thetaT= %.5f thetaW= %.5f TEST= NA \n",Th.theta1_Dtaj,Th.theta2_Dtaj);
                printf("FU&LI_D*: thetaW= %.5f thetaFL*= %.5f TEST= NA \n",Th.theta1_Dfl,Th.theta2_Dfl);
                printf("FU&LI_F*: thetaT= %.5f thetaFL*= %.5f TEST= NA \n",Th.theta1_Ffl,Th.theta2_Ffl);
                if(nsam <= 4)
                    printf("ACHAZ_Y*: thetaTns*= %.5f thetaWns*= %.5f TEST= NA \n",Th.theta1_Yach,Th.theta2_Yach);
                else
                    printf("ACHAZ_Y*: thetaTns*= NA thetaWns*= NA TEST= NA \n");
            }
            else {
                printf("TAJIMA_D: thetaT= %.5f thetaW= %.5f TEST= %.5f \n",Th.theta1_Dtaj,Th.theta2_Dtaj,statistics.Dtaj);
                printf("FU&LI_D*: thetaW= %.5f thetaFL*= %.5f TEST= %.5f \n",Th.theta1_Dfl,Th.theta2_Dfl,statistics.Dfl);
                printf("FU&LI_F*: thetaT= %.5f thetaFL*= %.5f TEST= %.5f \n",Th.theta1_Ffl,Th.theta2_Ffl,statistics.Ffl);
                if(statistics.Yach == -10000)
                    printf("ACHAZ_Y*: thetaTns*= %.5f thetaWns*= %.5f TEST= NA \n",Th.theta1_Yach,Th.theta2_Yach);
                else
                    printf("ACHAZ_Y*: thetaTns*= %.5f thetaWns*= %.5f TEST= %.5f \n",Th.theta1_Yach,Th.theta2_Yach,statistics.Yach);
            }
            /*printf("TAJIMA_D_NO_SINGL: thetaTns= %.5f thetaWns= %.5f TEST= %.5f \n",Th.theta1_Dtaj0,Th.theta2_Dtaj0,statistics.Dtaj0);*/
            
            printf("\n");
            fflush(stdin);
        }
    }
    else {
        usage();
        exit(0);
    }
    
    printf("Program finished.\n");
    return 0;
}

void usage(void)
{
    printf("Flags:\n");
    printf("      -f [outgroup_presence(1/0)] [#nsam] [freq1] [freq2] ... \n");
    printf("      -s [outgroup_presence(1/0)] [#nsam] [#freq:value] [#freq:value] ... \n");
    printf("      -d [seed] (use only in case large sample sizes, nsam > %d) \n",SAMPLE_LARGE);
    printf("      -h [help] \n");
    
}

/*Frequency tests from Achaz Genetics 2009: outgroup*/
double freqtesto_achaz(long int sample_size,long int *fr,long int singleton,double *w1,double *w2,double *theta1,double *theta2) /* nomes outgroup */
{
    long int i,j,ss,it,nit;
    double Th1,Th2,Test,Thw,Thw2,*ww;
    double sumw1,sumw2,sumww,omi;
    double alfan,betan,alfat,betat;
    double omegai(long int,long int,double *,double *,double,double);
    double sigmaii(long int,long int),sigmaij(long int,long int,long int);
    double an(long int),bn(long int,long int);
    double ran1();
    double betan2,betat2;
    long int max(double,double);
    
    if(sample_size < 2) return(-10000);
    
    ww = (double *)calloc(sample_size,sizeof(double));
    
    i=an(sample_size);/*just calculate all an*/
    
    Th1 = 0.;
    sumw1 = 0.;
    for(i=1;i<sample_size;i++) {
        Th1 += ((double)*(fr+i))*((double)i*(double)w1[i]);
        sumw1 += w1[i];
    }
    Th1 /= sumw1;
    *theta1 = Th1;
    
    Th2 = 0.;
    sumw2 = 0.;
    for(i=1;i<sample_size;i++) {
        Th2 += ((double)*(fr+i))*((double)i*(double)w2[i]);
        sumw2 += w2[i];
    }
    Th2 /= sumw2;
    *theta2 = Th2;
    
    if(Th1 == 0. && Th2 == 0.) return(-10000);
    
    Thw = 0.;
    sumww = 0.;
    for(i=1;i<sample_size;i++) {
        if(i==1) ss = singleton;
        else ss = 1;
        ww[i] = (double)1/(double)i * (double)ss;
        Thw += (double)*(fr+i)*(double)i*ww[i];
        sumww += ww[i];
    }
    Thw /= sumww;
    
    sumw1=0.0;
    for(i=1;i<sample_size;i++)
        sumw1 += w1[i];
    sumw2=0.0;
    for(i=1;i<sample_size;i++)
        sumw2 += w2[i];

    alfan = 0.;
    for(i=1;i<sample_size;i++) {
        omi = omegai(sample_size,i,w1,w2,sumw1,sumw2);
        alfan += i*(omi*omi);
    }
    
    if(sample_size > SAMPLE_LARGE){
#if BOOTSTRAP1==1
        betan = 0.;
        betan2 = 0;
        nit = max(1,(double)NITER/((double)sample_size-1));
        for(i=1;i<sample_size;i++) {
            omi = omegai(sample_size,i,w1,w2,sumw1,sumw2);
            betan += i*i * (omi*omi) * sigmaii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{
                    j=(long int)floor((double)ran1() * (double)(sample_size-1))+1;
                }while(i==j);
                betan2 += 2.0 * i*j * omegai(sample_size,i,w1,w2,sumw1,sumw2) * omegai(sample_size,j,w1,w2,sumw1,sumw2) * sigmaij(sample_size,j,i);
            }
        }
        betan2 /= nit*(sample_size-1);
        betan2 *= ((sample_size-1)*(sample_size-2)/2.0)-sample_size-1;
        betan  += betan2;
#else
        betan = 0.;
        for(it=0;it<NITER;it++) {
         do {
            i=(long int)floor(ran1() * (sample_size-1))+1;
            j=(long int)floor(ran1() * (sample_size-1))+1;
         }while(i==j);
            betan += 2.0 * i*j * omegai(sample_size,i,w1,w2,sumw1,sumw2) * omegai(sample_size,j,w1,w2,sumw1,sumw2) * sigmaij(sample_size,j,i);
        }
        betan /= NITER;
        betan *= ((sample_size-1)*(sample_size-2)/2.0)-sample_size-1;
        for(i=1;i<sample_size;i++) {
            omi = omegai(sample_size,i,w1,w2,sumw1,sumw2);
            betan += i*i * (omi*omi) * sigmaii(sample_size,i);
        }
#endif
    }
    else {
        betan = 0.;
        for(i=1;i<sample_size;i++) {
            omi = omegai(sample_size,i,w1,w2,sumw1,sumw2);
            betan += i*i * (omi*omi) * sigmaii(sample_size,i);
            /*printf("\nbetan=%f\tsigmaii[n=%d,i=%d]=%f",betan,sample_size,i,sigmaii(sample_size,i));*/
            for(j=i+1;j<sample_size;j++) {
                betan += 2.0 * i*j * omegai(sample_size,i,w1,w2,sumw1,sumw2) * omegai(sample_size,j,w1,w2,sumw1,sumw2) * sigmaij(sample_size,j,i);
                /*printf("\nbetan=%f\tsigmaij[n=%d,i=%d,j=%d]=%f",betan,sample_size,i,j,sigmaij(sample_size,j,i));*/
            }
        }
    }
    
    /*Theta2*/
    alfat = 0.;
    for(i=1;i<sample_size;i++) {
        alfat += (ww[i]/sumww * ww[i]/sumww)*i;
    }
    if(sample_size > SAMPLE_LARGE){
#if BOOTSTRAP1==1
        betat = 0.;
        betat2 = 0.;
        nit = max(1,(double)NITER/((double)sample_size-1));
        for(i=1;i<sample_size;i++) {
            betat += i*i * (ww[i]/sumww * ww[i]/sumww) * sigmaii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{
                    j=(long int)floor(ran1() * (sample_size-1))+1;
                }while(i==j);
                betat2 += 2.0 * i*j * ww[i]/sumww * ww[j]/sumww * sigmaij(sample_size,j,i);
            }
        }
        betat2 /= nit*(sample_size-1);
        betat2 *= ((sample_size-1)*(sample_size-2)/2.0)-sample_size-1;
        betat  += betat2;
#else
        betat = 0.;
        for(it=0;it<NITER;it++) {
         do {
            i=(long int)floor(ran1() * (sample_size-1))+1;
            j=(long int)floor(ran1() * (sample_size-1))+1;
         }while(i==j);
            betat += 2.0 * i*j * ww[i]/sumww * ww[j]/sumww * sigmaij(sample_size,j,i);
        }
        betat /= NITER;
        betat *= ((sample_size-1)*(sample_size-2)/2.0)-sample_size-1;
        for(i=1;i<sample_size;i++) {
            betat += i*i * (ww[i]/sumww * ww[i]/sumww) * sigmaii(sample_size,i);
        }
#endif
    }
    else {
        betat = 0.;
        for(i=1;i<sample_size;i++) {
            betat += i*i * (ww[i]/sumww * ww[i]/sumww) * sigmaii(sample_size,i);
            for(j=i+1;j<sample_size;j++) {
                betat += 2.0 * i*j * ww[i]/sumww * ww[j]/sumww * sigmaij(sample_size,j,i);
            }
        }
    }
    
    Thw2 = (Thw*Thw - alfat*Thw)/(1.0 + betat);
    if((sqrt(alfan*Thw + betan*Thw2)) == 0)
        return(-10000);
    
    /*Test*/
    Test = (Th1 - Th2)/(sqrt(alfan*Thw + betan*Thw2));
    
    free(ww);
    if (fabs(Test) < 1.0E-15)
        return 0.0;
    
    return Test;
}
/*Frequency tests from Achaz Genetics 2009: NO outgroup*/
double freqtestn_achaz(long int sample_size,long int *fr,long int singleton,double *w1,double *w2,double *theta1,double *theta2) /* NO outgroup */
{
    long int i,j,ss,it,nit;
    double Th1,Th2,Test,Thw,Thw2,*ww;
    double sumw1,sumw2,sumww,omi,omj,psi,psj;
    double psii(long int,long int),rhoii(long int,long int),rhoij(long int,long int,long int);
    double alfan,betan,alfat,betat;
    double omegain(long int,long int,double *,double *,double,double);
    double sigmaii(long int,long int),sigmaij(long int,long int,long int);
    double an(long int),bn(long int,long int);
    double betan2,betat2;
    long int max(double,double);
    double ran1();
   
    if(sample_size < 2) return(-10000);
    
    ww = (double *)calloc(sample_size,sizeof(double));
    i=an(sample_size);/*just calculate all an*/
    
    Th1 = 0.;
    sumw1 = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
        Th1 += ((double)*(fr+i))*(double)w1[i]/((double)psii(sample_size,i));
        sumw1 += w1[i];
    }
    Th1 /= sumw1;
    *theta1 = Th1;
    
    Th2 = 0.;
    sumw2 = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
        Th2 += ((double)*(fr+i))*(double)w2[i]/((double)psii(sample_size,i));
        sumw2 += w2[i];
    }
    Th2 /= sumw2;
    *theta2 = Th2;
    
    if(Th1 == 0. && Th2 == 0.) return(-10000);

    Thw = 0.;
    sumww = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
        if(i==1) ss = singleton;
        else ss = 1;
        if(i == sample_size-i) ww[i] = (double)sample_size/((double)i*(double)(sample_size - i)*2.0)*(double)ss;
        else ww[i] = (double)sample_size/((double)i*(double)(sample_size - i)*1.0)*(double)ss;
        Thw += ((double)*(fr+i))*ww[i]/((double)psii(sample_size,i));
        sumww += ww[i];
    }
    Thw /= sumww;
    
    sumw1=0.0;
    for(i=1;i<=floor(sample_size/2);i++)
        sumw1 += w1[i];
    sumw2=0.0;
    for(i=1;i<=floor(sample_size/2);i++)
        sumw2 += w2[i];

    alfan = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
        omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
        psi = psii(sample_size,i);
        alfan += (omi*omi)/psi;
    }
    
    if(floor(sample_size/*/2.*/) > SAMPLE_LARGE) {
#if BOOTSTRAP1==1
        betan = 0.;
        betan2 = 0.;
        nit = max(1,(double)NITER/((double)sample_size/2.));
        for(i=1;i<=floor(sample_size/2);i++) {
            omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
            psi = psii(sample_size,i);
            betan += omi/psi * omi/psi * rhoii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{j=(long int)floor(ran1() * (sample_size/2.))+1;}while(i==j);
                omj = omegain(sample_size,j,w1,w2,sumw1,sumw2);
                psj = psii(sample_size,j);
                betan2 += 2.0 * omi/psi * omj/psj * rhoij(sample_size,j,i);
            }
        }
        betan2 /= nit*(sample_size/2.);
        betan2 *= ((sample_size/2.)*(sample_size/2.-1)/2.0)-sample_size/2.;
        betan  += betan2;
#else
        betan = 0.;
        for(it=0;it<NITER;it++) {
         do {
            i=(long int)floor(ran1() * (sample_size/2.))+1;
            j=(long int)floor(ran1() * (sample_size/2.))+1;
         }while(i==j);
            omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
            omj = omegain(sample_size,j,w1,w2,sumw1,sumw2);
            psi = psii(sample_size,i);
            psj = psii(sample_size,j);
            betan += 2.0 * omi/psi * omj/psj * rhoij(sample_size,j,i);
        }
        betan /= NITER;
        betan *= ((sample_size/2.)*(sample_size/2.-1)/2.0)-sample_size/2.;
        for(i=1;i<=floor(sample_size/2);i++) {
            omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
            psi = psii(sample_size,i);
            betan += omi/psi * omi/psi * rhoii(sample_size,i);
        }
#endif
    }
    else {
        betan = 0.;
        for(i=1;i<=floor(sample_size/2);i++) {
            omi = omegain(sample_size,i,w1,w2,sumw1,sumw2);
            psi = psii(sample_size,i);
            betan += omi/psi * omi/psi * rhoii(sample_size,i);
            for(j=i+1;j<=floor(sample_size/2);j++) {
                omj = omegain(sample_size,j,w1,w2,sumw1,sumw2);
                psj = psii(sample_size,j);
                betan += 2.0 * omi/psi * omj/psj * rhoij(sample_size,j,i);
                /*printf("\nrhoij[n=%d,j=%d,i=%d] = %f",sample_size,j,i,rhoij(sample_size,j,i));*/
            }
        }
    }
    
    /*Theta2*/
    alfat = 0.;
    for(i=1;i<=floor(sample_size/2);i++) {
        psi = psii(sample_size,i);
        alfat += (ww[i]/sumww * ww[i]/sumww)/psi;
    }
    if(floor(sample_size/*/2.*/) > SAMPLE_LARGE) {
#if BOOTSTRAP1==1
        betat = 0.;
        betat2 = 0.;
        nit = max(1,(double)NITER/((double)sample_size/2.));
        for(i=1;i<=floor(sample_size/2);i++) {
            psi = psii(sample_size,i);
            betat += (ww[i]/sumww)/psi * (ww[i]/sumww)/psi * rhoii(sample_size,i);
            for(it=0;it<nit;it++) {
                do{j=(long int)floor(ran1() * (sample_size/2.))+1;}while(i==j);
                psj = psii(sample_size,j);
                betat2 += 2.0 * (ww[i]/sumww)/psi * (ww[j]/sumww)/psj * rhoij(sample_size,j,i);
            }
        }
        betat2 /= nit*(sample_size/2.);
        betat2 *= ((sample_size/2.)*(sample_size/2.-1)/2.0)-sample_size/2.;
        betat  += betat2;
#else
        betat = 0.;
        for(it=0;it<NITER;it++) {
         do {
            i=(long int)floor(ran1() * (sample_size/2.))+1;
            j=(long int)floor(ran1() * (sample_size/2.))+1;
         }while(i==j);
            psi = psii(sample_size,i);
            psj = psii(sample_size,j);
            betat += 2.0 * (ww[i]/sumww)/psi * (ww[j]/sumww)/psj * rhoij(sample_size,j,i);
        }
        betat /= NITER;
        betat *= ((sample_size/2.)*(sample_size/2.-1)/2.0)-sample_size/2.;
        for(i=1;i<=floor(sample_size/2);i++) {
            psi = psii(sample_size,i);
            betat += (ww[i]/sumww)/psi * (ww[i]/sumww)/psi * rhoii(sample_size,i);
        }
#endif
    }
    else {
        betat = 0.;
        for(i=1;i<=floor(sample_size/2);i++) {
            psi = psii(sample_size,i);
            betat += (ww[i]/sumww)/psi * (ww[i]/sumww)/psi * rhoii(sample_size,i);
            for(j=i+1;j<=floor(sample_size/2);j++) {
                psj = psii(sample_size,j);
                betat += 2.0 * (ww[i]/sumww)/psi * (ww[j]/sumww)/psj * rhoij(sample_size,j,i);
            }
        }
    }
    Thw2 = (Thw*Thw - alfat*Thw)/(1.0 + betat);
    if((sqrt(alfan*Thw + betan*Thw2)) == 0)
        return(-10000);
        
    /*Test*/
    Test = (Th1 - Th2)/(sqrt(alfan*Thw + betan*Thw2));
    
    free(ww);
    if (fabs(Test) < 1.0E-15)
        return 0.0;
    
    return Test;
}
double an(long int n)
{
    static double *an_m=0;
    static long int nsam=0;
    long int i;
    
    if(nsam == 0) {
        an_m = (double *)calloc(n+2, sizeof(double));
        for(i=1;i<n+1;i++)
            an_m[i+1] = an_m[i] + 1.0/(double)i;
        nsam = n;
    } else {
        if(nsam < n) {
            an_m = (double *)realloc(an_m,(n+2)*sizeof(double));
            for(i=nsam+1;i<n+1;i++)
                an_m[i+1] = an_m[i] + 1.0/(double)i;
            nsam = n;
        }
    }
    return an_m[n];
}
double a2n(long int n)
{
    double a2ni = 0.;
    int i;
    for(i=1;i<n;i++)
        a2ni += (1.0/((double)i*(double)i));
    return a2ni;
}
double bn(long int n,long int i)
{
    double an(long int);
    double bni;
    bni = (2.0*(double)n)/((double)(n-i+1)*(double)(n-i)) * (an(n+1) - an(i)) - 2.0/((double)(n-i));
    return bni;
}
double sigmaii(long int n,long int i)
{
    double an(long int),bn(long int,long int);
    double sii=0;
    if(2*i<n)
        sii = bn(n,i+1);
    if(2*i==n)
        sii = 2.0 * (an(n)-an(i))/((double)(n-i)) - 1.0/((double)i*(double)i);
    if(2*i>n)
        sii = bn(n,i) - 1.0/((double)i*(double)i);
    return sii;
}
double sigmaij(long int n,long int i,long int j)
{
    double an(long int),bn(long int,long int);
    double sigmaii(long int,long int);
    double sij=0;
    long int ii,jj;
    if(i < j) {
        ii = j;
        jj = i;
    }
    else {
        if(i==j) {
            return(sigmaii(n,i));
        }
        else {
            ii = i;
            jj = j;
        }
    }
    if((ii+jj)<n) 
        sij = (bn(n,ii+1) - bn(n,ii))/2.0;
    if((ii+jj)==n) 
        sij = (an(n)-an(ii))/((double)(n-ii)) + (an(n)-an(jj))/((double)(n-jj)) - (bn(n,ii) + bn(n,jj+1))/2.0 - 1.0/((double)ii*(double)jj);
    if((ii+jj)>n) 
        sij = (bn(n,jj) - bn(n,jj+1))/2.0 - 1.0/((double)ii*(double)jj);
    return sij;
}
double omegai(long int n,long int i,double *w1,double *w2,double sumw1,double sumw2)
{
    double omi;
    /*
    long int x;
    double sumw1,sumw2;
    
    sumw1=0.0;
    for(x=1;x<n;x++) sumw1 += w1[x];
    sumw2=0.0;
    for(x=1;x<n;x++) sumw2 += w2[x];
    */
    omi = w1[i]/sumw1 - w2[i]/sumw2;
    return omi;
}
double psii(long int n,long int i)
{
    double psi;
    long int krond;
    
    if((long int)i==(long int)(n-i)) krond = 1;
    else krond = 0;
    psi= (double)n/((double)(1.+krond)*(double)i*(double)(n-i));
    return psi;
}
double rhoii(long int n,long int i)
{
    double sigmaii(long int,long int);
    double sigmaij(long int,long int,long int);
    double rhoi;
    long int krond;
    
    if((long int)i==(long int)(n-i)) 
        krond = 1;
    else 
        krond = 0;
    
    rhoi  = (sigmaii(n,i)+sigmaii(n,n-i)+2.0*sigmaij(n,i,n-i));
    rhoi /= ((1.0 + krond) * (1.0 + krond));
    return rhoi;
}
double rhoij(long int n,long int i,long int j)
{
    double sigmaii(long int,long int);
    double sigmaij(long int,long int,long int);
    double rhoj;
    long int krondi,krondj;
    
    if(i==(n-i)) krondi = 1;
    else krondi = 0;
    if(j==(n-j)) krondj = 1;
    else krondj = 0;
    
    rhoj  = (sigmaij(n,i,j)+sigmaij(n,i,n-j)+sigmaij(n,n-i,j)+sigmaij(n,n-i,n-j));
    rhoj /= ((1.0 + krondi) * (1.0 + krondj));
    return rhoj;
}
double omegain(long int n,long int i,double *w1,double *w2,double sumw1,double sumw2)
{
    double omi;
    /*
    long int x;
    double sumw1,sumw2;
    
    sumw1=0.0;
    for(x=1;x<=floor(n/2);x++) sumw1 += w1[x];
    sumw2=0.0;
    for(x=1;x<=floor(n/2);x++) sumw2 += w2[x];
    */
    omi = w1[i]/sumw1 - w2[i]/sumw2;
    return omi;
}

long int max(double x,double y) {
    if(x > y) return x;
    else return y;
}

