#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "kd.h"
#include "smooth.h"


#define BIGNOTQUITEMAXFLOAT		((float)1.0e+37)


void usage(void)
{
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"smooth [-s <nSmooth>[dgs]] [-b <nBucket>] [-g]\n");
	fprintf(stderr,"   [-o <Output Name>] [-p <xyzPeriod>]\n");
	fprintf(stderr,"   [-px <xPeriod>] [-py <yPeriod>] [-pz <zPeriod>]\n");
	fprintf(stderr,"   [-do <MarkFile>]\n");
	fprintf(stderr,"   [density] [meanvel] [speed] [veldisp] [mach]\n");
	fprintf(stderr,"   [phase] [hsmooth] [divv] [all] [null]\n\n");
	fprintf(stderr,"Input taken from stdin in tipsy binary format.\n");
	fprintf(stderr,"SEE MAN PAGE: smooth(1) for more information.\n");
	exit(1);
	}

int main(int argc,char **argv)
{
	KD kd;
	SMX smx;
	int nBucket,nSmooth,i,j;
	FILE *fp;
	char ach[80],achFile[80],achMark[80];
	float fPeriod[3];
	int bDensity,bMeanVel,bVelDisp,bPhase,bMach,bSpeed,bNull,bSym;
	int bHsmooth;
	int bDark,bGas,bStar;
	int bMark;
	int bDivv;
	char *p,*q;
	
   	nBucket = 16;
	nSmooth = 64;
	bDensity = 0;
	bMeanVel = 0;
	bVelDisp = 0;
	bSpeed = 0;
	bMach = 0;
	bPhase = 0;
	bHsmooth = 0;
	bDivv = 0;
	bNull = 0;
	bSym = 1;
	bDark = 1;
	bGas = 1;
	bStar = 1;
	bMark = 0;
	strcpy(achFile,"smooth");
	i = 1;
	for (j=0;j<3;++j) fPeriod[j] = BIGNOTQUITEMAXFLOAT;
	while (i < argc) {
		if (!strcmp(argv[i],"-b")) {
			++i;
			if (i >= argc) usage();
			nBucket = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-s")) {
			++i;
			if (i >= argc) usage();
			p = argv[i];
			while (isdigit(*p)) ++p;
			q = p;
			if (isalpha(*p)) {
				bDark = 0;
				bGas = 0;
				bStar = 0;
				}
			while (isalpha(*p)) {
				switch (*p) {
				case 'd':
					bDark = 1;
					break;
				case 'g':
					bGas = 1;
					break;
				case 's':
					bStar = 1;
					break;
				default:
					usage();
					}
				++p;
				}
			*q = 0;
			nSmooth = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-o")) {
			++i;
			if (i >= argc) usage();
			strcpy(achFile,argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-do")) {
			++i;
			if (i >= argc) usage();
			strcpy(achMark,argv[i]);
			bMark = 1;
			bSym = 0;	/* Symmetrical kernal is inconsistent here! */
			++i;
			}
		else if (!strcmp(argv[i],"-g")) {
			bSym = 0;
			++i;
			}
		else if (!strcmp(argv[i],"-p")) {
			++i;
			if (i >= argc) usage();
			fPeriod[0] = atof(argv[i]);
			fPeriod[1] = atof(argv[i]);
			fPeriod[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-px")) {
			++i;
			if (i >= argc) usage();
			fPeriod[0] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-py")) {
			++i;
			if (i >= argc) usage();
			fPeriod[1] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-pz")) {
			++i;
			if (i >= argc) usage();
		    fPeriod[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"density")) {
			bDensity |= 3;
		    ++i;
			}
		else if (!strcmp(argv[i],"meanvel")) {
			bDensity |= 1;
			bMeanVel |= 3;
			++i;
			}
		else if (!strcmp(argv[i],"veldisp")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bVelDisp |= 3;
			++i;
			}
		else if (!strcmp(argv[i],"phase")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bVelDisp |= 1;
			bDivv |= 1;
			bPhase |= 2;
			++i;
			}
		else if (!strcmp(argv[i],"mach")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bVelDisp |= 1;
			bMach |= 2;
			++i;
			}
		else if (!strcmp(argv[i],"speed")) {
			bDensity |= 1;
			bMeanVel |= 1;
			bSpeed |= 2;
			++i;
			}
		else if (!strcmp(argv[i],"divv")) {
			bDensity |= 1;
			bDivv |= 3;
			++i;
			}
		else if (!strcmp(argv[i],"hsmooth")) {
			bNull |= 1;
			bHsmooth |= 2;
			++i;
			}
		else if (!strcmp(argv[i],"null")) {
			bNull |= 1;
			++i;
			}
		else if (!strcmp(argv[i],"all")) {
			bDensity |= 3;
			bMeanVel |= 3;
			bVelDisp |= 3;
			bDivv |= 3;
			bPhase |= 2;
			bMach |= 2;
			bSpeed |= 2;
			bHsmooth |= 2;
			++i;
			}
		else usage();
		}
	kdInit(&kd,nBucket);
	kdReadTipsy(kd,stdin,bDark,bGas,bStar);
	if (bMark) kdInMark(kd,achMark);
	kdBuildTree(kd);
	smInit(&smx,kd,nSmooth,fPeriod);
	if (bNull&1) {
		smSmooth(smx,smNull);
		}
	if (bSym) {
		if (bDensity&1) {
			if (bNull&1) smReSmooth(smx,smDensitySym);
			else smSmooth(smx,smDensitySym);
			}
		if (bMeanVel&1) smReSmooth(smx,smMeanVelSym);
		if (bDivv&1) smReSmooth(smx,smDivvSym);
		if (bVelDisp&1) smReSmooth(smx,smVelDispNBSym);
		}
	else {
		if (bDensity&1) {
			if (bNull&1) smReSmooth(smx,smDensity);
			else smSmooth(smx,smDensity);
			}
		if (bMeanVel&1) smReSmooth(smx,smMeanVel);
		if (bDivv&1) smReSmooth(smx,smDivv);
		if (bVelDisp&1) smReSmooth(smx,smVelDispNB);
		}
	kdOrder(kd);
	if (bDensity&2) {
		strcpy(ach,achFile);
		strcat(ach,".den");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutDensity(smx,fp);
		fclose(fp);
		}
	if (bMeanVel&2) {
		strcpy(ach,achFile);
		strcat(ach,".mvl");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutMeanVel(smx,fp);
		fclose(fp);
		}
	if (bSpeed&2) {
		strcpy(ach,achFile);
		strcat(ach,".spd");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutSpeed(smx,fp);
		fclose(fp);
		}
	if (bVelDisp&2) {
		strcpy(ach,achFile);
		strcat(ach,".dsp");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutVelDisp(smx,fp);
		fclose(fp);
		}
	if (bDivv&2) {
		strcpy(ach,achFile);
		strcat(ach,".dvv");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutDivv(smx,fp);
		fclose(fp);
		}
	if (bMach&2) {
		strcpy(ach,achFile);
		strcat(ach,".mch");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutMach(smx,fp);
		fclose(fp);
		}
	if (bPhase&2) {
		strcpy(ach,achFile);
		strcat(ach,".phs");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutPhase(smx,fp);
		fclose(fp);
		}
	if (bHsmooth&2) {
		strcpy(ach,achFile);
		strcat(ach,".hsm");
		fp = fopen(ach,"w");
		assert(fp != NULL);
		smOutHsmooth(smx,fp);
		fclose(fp);
		}
	smFinish(smx);
	kdFinish(kd);
	return 0;
	}
	

