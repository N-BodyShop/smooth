#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include "kd.h"
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>


int xdr_header(XDR *xdrs, struct dump *header)
{
	int pad=0;
  
	if (xdr_double(xdrs,&header->time) != TRUE) return 0;
	if (xdr_int(xdrs,&header->nbodies) != TRUE) return 0;
	if (xdr_int(xdrs,&header->ndim) != TRUE) return 0;
	if (xdr_int(xdrs,&header->nsph) != TRUE) return 0;
	if (xdr_int(xdrs,&header->ndark) != TRUE) return 0;
	if (xdr_int(xdrs,&header->nstar) != TRUE) return 0;
	if (xdr_int(xdrs,&pad) != TRUE) return 0;
	return 1;
	}


int xdr_gas(XDR *xdrs,struct gas_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->rho) != TRUE) return 0;
	if (xdr_float(xdrs,&p->temp) != TRUE) return 0;
	if (xdr_float(xdrs,&p->hsmooth) != TRUE) return 0;
	if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


int xdr_dark(XDR *xdrs,struct dark_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


int xdr_star(XDR *xdrs,struct star_particle *p)
{
	if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
	if (xdr_float(xdrs,&p->metals) != TRUE) return 0;
	if (xdr_float(xdrs,&p->tform) != TRUE) return 0;
	if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
	if (xdr_float(xdrs,&p->phi) != TRUE) return 0;
	return 1;
	}  


#define MAX_ROOT_ITTR	32


void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket)
{
	KD kd;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	kd->p = NULL;
	kd->kdNodes = NULL;
	*pkd = kd;
	return(1);
	}


int kdReadTipsy(KD kd,FILE *fp,int bNative,int bDark,int bGas,int bStar)
{
	int i,j,nCnt,ret;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;
	XDR xdr;

	if (bNative) ret = fread(&h,sizeof(struct dump),1,fp);
	else {
		xdrstdio_create(&xdr,fp,XDR_DECODE);
		ret = xdr_header(&xdr,&h);
		}	
	if (ret != 1) {
		fprintf(stderr,"Error reading header from input file.\n");
		exit(1);
		}
	kd->nParticles = h.nbodies;
	kd->nDark = h.ndark;
	kd->nGas = h.nsph;
	kd->nStar = h.nstar;
	kd->fTime = h.time;
	kd->nActive = 0;
	if (bDark) kd->nActive += kd->nDark;
	if (bGas) kd->nActive += kd->nGas;
	if (bStar) kd->nActive += kd->nStar;
	kd->bDark = bDark;
	kd->bGas = bGas;
	kd->bStar = bStar;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<h.nsph;++i) {
		if (bNative) ret = fread(&gp,sizeof(struct gas_particle),1,fp);
		else ret = xdr_gas(&xdr,&gp);
		if (ret != 1) {
			fprintf(stderr,"Error reading gas particle from input file.\n");
			exit(1);
			}
		if (bGas) {
			kd->p[nCnt].fMass = gp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = gp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = gp.vel[j];
			kd->p[nCnt].fSmooth = 0;
			++nCnt;
			}
		}
	for (i=0;i<h.ndark;++i) {
		if (bNative) ret = fread(&dp,sizeof(struct dark_particle),1,fp);
		else ret = xdr_dark(&xdr,&dp);
		if (ret != 1) {
			fprintf(stderr,"Error reading dark particle from input file.\n");
			exit(1);
			}
		if (bDark) {
			kd->p[nCnt].fMass = dp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = dp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = dp.vel[j];
			kd->p[nCnt].fSmooth = 0;
			++nCnt;
			}
		}
	for (i=0;i<h.nstar;++i) {
		if (bNative) ret = fread(&sp,sizeof(struct star_particle),1,fp);
		else ret = xdr_star(&xdr,&sp);
		if (ret != 1) {
			fprintf(stderr,"Error reading star particle from input file.\n");
			exit(1);
			}
		if (bStar) {
			kd->p[nCnt].fMass = sp.mass;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = sp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = sp.vel[j];
			kd->p[nCnt].fSmooth = 0;
			++nCnt;
			}
		}
	if (!bNative) xdr_destroy(&xdr);
	return(kd->nParticles);
	}


void kdInMark(KD kd,char *pszFile)
{
	FILE *fp;
	char ach[80];
	int i,iCnt,iDum;

	fp = fopen(pszFile,"r");
	if (!fp) {
		fprintf(stderr,"Could not open mark array, %s\n",pszFile);
		exit(1);
		}
	fgets(ach,80,fp);	/* ignore the array header! */
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		kd->p[iCnt++].iMark = 0;
		}
	for (i=0;i<kd->nDark;++i) {
		kd->p[iCnt++].iMark = 0;
		}
	for (i=0;i<kd->nStar;++i) {
		kd->p[iCnt++].iMark = 0;
		}
	while(1) {
	    fscanf(fp, "%d", &iDum);
	    if(feof(fp))
		break;
	    iDum--;
	    assert(iDum < kd->nParticles);
	    kd->p[iDum].iMark = 1;
	    }
	fclose(fp);
	}


void kdSelect(KD kd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	double v;
	int i,j;

	p = kd->p;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void kdUpPass(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		kdUpPass(kd,l);
		kdUpPass(kd,u);
		kdCombine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->p[u].r[j];
			c[iCell].bnd.fMax[j] = kd->p[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->p[pj].r[j];
				if (kd->p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->p[pj].r[j];
				}
			}
		}
	}


void kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes != NULL) free(kd->kdNodes);
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->p[0].r[j];
		bnd.fMax[j] = kd->p[0].r[j];
		}
	for (i=1;i<kd->nActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->p[i].r[j]) 
				bnd.fMin[j] = kd->p[i].r[j];
			else if (bnd.fMax[j] < kd->p[i].r[j])
				bnd.fMax[j] = kd->p[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelect(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->p[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	kdUpPass(kd,ROOT);
	}


int cmpParticles(const void *v1,const void *v2)
{
	PARTICLE *p1=(PARTICLE *)v1,*p2=(PARTICLE *)v2;
	
	return(p1->iOrder - p2->iOrder);
	}


void kdOrder(KD kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLE),cmpParticles);
	}


void kdFinish(KD kd)
{
	free(kd->p);
	free(kd->kdNodes);
	free(kd);
	}

