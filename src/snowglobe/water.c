/*
 * Compiz cube snowglobe plugin
 *
 * snowglobe.c
 *
 * This is an example plugin to show how to render something inside
 * of the transparent cube
 *
 * Copyright : (C) 2007 by Dennis Kasprzyk
 * E-mail    : onestone@opencompositing.org
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 */

/* Uses water amplitude calculation by David Mikos */

/*
 * At the moment dealing with deformation has introduced lots of
 * redundant/repeated sections code and generally messy code.
 * This needs some tidying up.
 */ 

#include "snowglobe-internal.h"
#include "math.h"
#include "snowglobe_options.h"

static void
genTriMesh (Vertex       *vertices,
	    unsigned int *indices,
	    unsigned int idxBase,
	    unsigned int subdiv,
	    Vertex       a,
	    Vertex       b,
	    Vertex       c)
{
    int          nRow;
    Vertex       *v;
    unsigned int *idx;
    int          i, j, k;
    int          tr, br;

    float    vab[3], vac[3], rb[3], re[3], ri[3];

    if (subdiv < 0)
	return;
    if (!vertices || !indices)
	return;

    nRow = (subdiv)?(2 << (subdiv - 1)) + 1 : 2;

    v =   vertices;
    idx = indices;


    for (i = 1; i < nRow; i++)
    {
	tr = (i * (i - 1)) / 2;
	br = (i * (i + 1)) / 2;

	for (j = 0; j < (2 * i) - 1; j++)
	{
	    if (j & 1)
	    {
		k = (j - 1) / 2;
		idx[(j * 3)] = idxBase + tr + k + 1;
		idx[(j * 3) + 1] = idxBase + tr + k;
		idx[(j * 3) + 2] = idxBase + br + k + 1;
	    }
	    else
	    {
		k = j / 2;
		idx[(j * 3)] = idxBase + tr + k;
		idx[(j * 3) + 1] = idxBase + br + k;
		idx[(j * 3) + 2] = idxBase + br + k + 1;
	    }
	}
	idx += ((2 * i) - 1) * 3;
    }

    for (i = 0; i < 3; i++)
    {
	vab[i] = b.v[i] - a.v[i];
	vab[i] /= nRow - 1.0;
	vac[i] = c.v[i] - a.v[i];
	vac[i] /= nRow - 1.0;
    }

    v[0] = a;

    for (i = 1; i < nRow; i++)
    {
	br = (i * (i + 1)) / 2;
	for (k = 0; k < 3; k++)
	{
	    rb[k] = a.v[k] + (i * vab[k]);
	    re[k] = a.v[k] + (i * vac[k]);
	    ri[k] = re[k] - rb[k];
	    ri[k] /= i;
	}
	for (j = 0; j <= i; j++)
	    for (k = 0; k < 3; k++)
		v[br + j].v[k] = (rb[k] + (j * ri[k]));
    }

}

static void
genTriWall (Vertex       *lVer,
	    Vertex       *hVer,
	    unsigned int *indices,
	    unsigned int idxBaseL,
	    unsigned int idxBaseH,
	    int          subdiv,
	    Vertex       a,
	    Vertex       b,
	    Vertex       c,
	    Vertex       d)
{
    int   nRow;
    int   i, k;
    float vab[3], vcd[3];

    if (subdiv < 0)
	return;
    if (!lVer || !hVer || !indices)
	return;

    nRow = pow (2, subdiv);

    for (i = 0; i < 3; i++)
    {
	vab[i] = b.v[i] - a.v[i];
	vab[i] /= nRow;
	vcd[i] = d.v[i] - c.v[i];
	vcd[i] /= nRow;
    }

    for (i = 0; i <= nRow; i++)
    {
	for (k = 0; k < 3; k++)
	{
	    lVer[i].v[k] = a.v[k] + (i * vab[k]);
	    hVer[i].v[k] = c.v[k] + (i * vcd[k]);
	}
    }

    for (i = 0; i < nRow; i++)
    {
	indices[(i * 6)] = idxBaseL + i;
	indices[(i * 6) + 1] = idxBaseH + i;
	indices[(i * 6) + 2] = idxBaseH + i + 1;
	indices[(i * 6) + 3] = idxBaseL + i + 1;
	indices[(i * 6) + 4] = idxBaseL + i;
	indices[(i * 6) + 5] = idxBaseH + i + 1;
    }
}

static Water *
genWater (int size, int sDiv, float distance, float bottom)
{

    Water  *w;
    int    i;
    float  ang, r, aStep;
    int    nVer, nRow, nIdx, nWVer, nWIdx;;
    float ratioRadiusToSideDist;

    Vertex a = {{ 0.0, 0.0, 0.0 }};
    Vertex b = {{ 0.0, 0.0, 0.0 }};
    Vertex c = {{ 0.0, 0.0, 0.0 }};
    Vertex d = {{ 0.0, bottom, 0.0 }};
    Vertex e = {{ 0.0, bottom, 0.0 }};

    Vertex       *wv;
    unsigned int *wi;

    if (sDiv < 0)
	return NULL;

    if (size < 3)
	return NULL;

    w = malloc (sizeof (Water));
    if (!w)
	return NULL;

    nRow = (sDiv) ? (2 << (sDiv - 1)) + 1 : 2;
    nVer = (nRow * (nRow + 1)) / 2;
    nIdx = pow (4, sDiv) * 3;

    nWIdx = pow (2, sDiv + 1) * 3;
    nWVer = pow (2, sDiv + 1) + 2;

    w->nVertices  = (nVer + nWVer) * size;
    w->nIndices   = (nIdx + nWIdx) * size;

    w->nSVer = nVer * size;
    w->nSIdx = nIdx * size;
    w->nWVer = nWVer * size;
    w->nWIdx = nWIdx * size;

    w->size     = size;
    w->distance = distance;
    w->sDiv     = sDiv;

    w->wave1 = 0.0;
    w->wave2 = 0.0;

    w->vertices = calloc (1,sizeof (Vertex) * w->nVertices);
    if (!w->vertices)
    {
	free (w);
	return NULL;
    }

    w->indices = calloc (1,sizeof (int) * w->nIndices);
    if (!w->indices)
    {
	free (w->vertices);
	free (w);
	return NULL;
    }

    w->vertices2 = NULL;
    w->indices2  = NULL;

    r = distance / cos (M_PI / size);
    ang = M_PI / size;
    aStep = 2 * M_PI / size;

    wv = w->vertices + (size * nVer);
    wi = w->indices + (size * nIdx);

    for (i = 0; i < size; i++)
    {
	d.v[0] = b.v[0] = sin (ang - aStep) * r;
	d.v[2] = b.v[2] = cos (ang - aStep) * r;

	e.v[0] = c.v[0] = sin (ang) * r;
	e.v[2] = c.v[2] = cos (ang) * r;

	genTriMesh (w->vertices + (i * nVer), w->indices + (i * nIdx),
		    i * nVer, sDiv, a, b, c);

	genTriWall (wv + (i * nWVer / 2), wv + ((i + size) * nWVer / 2),
		    wi + (i * nWIdx), (size * nVer) + (i * nWVer / 2),
		    (size * nVer) + ((i + size) * nWVer / 2), sDiv, b, c, d, e);

	ang += aStep;
    }

    return w;
}

void
freeWater (Water *w)
{
    if (!w)
	return;

    if (w->vertices)
	free (w->vertices);
    if (w->indices)
	free (w->indices);
    if (w->vertices2)
	free (w->vertices2);
    if (w->indices2)
	free (w->indices2);
}

static void
setAmplitude (Vertex *v,
	      float  bh,
	      float  wt,
	      float  swt,
	      float  wa,
	      float  swa,
	      float  wf,
	      float  swf)
{
    float  dx, dz, d, c1, c2;
    Vertex a, b;

    v->v[1] = bh + (wa * sinf (wt + wf * v->v[0] * v->v[2])) +
		   (swa * sinf (swt + swf * v->v[0] * v->v[2]));
    v->v[1] = MIN (0.5, MAX (-0.5, v->v[1]));

    c1 = wa * cosf (wt + wf * v->v[0] * v->v[2]) * wf;
    c2 = swa * cosf (swt + swf * v->v[0] * v->v[2]) * swf;

    dx = (c1 * v->v[2]) + (c2 * v->v[2]);
    dz = (c1 * v->v[0]) + (c2 * v->v[0]);

    a.v[0] = 10;
    a.v[2] = 0;
    b.v[0] = 0;
    b.v[2] = 10;

    a.v[1] = v->v[1] + (10 * dx);
    b.v[1] = v->v[1] + (10 * dz);

    v->n[0] = (b.v[1] * a.v[2]) - (b.v[2] * a.v[1]);
    v->n[1] = (b.v[2] * a.v[0]) - (b.v[0] * a.v[2]);
    v->n[2] = (b.v[0] * a.v[1]) - (b.v[1] * a.v[0]);

    d = sqrt ((v->n[0] * v->n[0]) + (v->n[1] * v->n[1]) + (v->n[2] * v->n[2]));

    v->n[0] /= d;
    v->n[1] /= d;
    v->n[2] /= d;

}

static void
deformCylinder(CompScreen *s, Water  *w, float progress)
{
    SNOWGLOBE_SCREEN (s);
    CUBE_SCREEN (s);

    int          nVer, nWVer, nRow, nRowS, subdiv;
    Vertex       *v;
    int          i, j, k, l;
    int          br;

    float  ang, r, aStep;
    
    Vertex       *wv;
    
    int bottom = -0.5, size = as->hsize;
 
    Vertex a = {{ 0.0, 0.0, 0.0 }};
    Vertex b = {{ 0.0, 0.0, 0.0 }};
    Vertex c = {{ 0.0, 0.0, 0.0 }};
    Vertex d = {{ 0.0, bottom, 0.0 }};
    Vertex e = {{ 0.0, bottom, 0.0 }};
    
    float    vab[3], vac[3], rb[3], re[3], ri[3];

    if (!w)
	return;
    if (w->sDiv < 0)
	return;
    if (!w->vertices)
	return;
    if (w->size!=size)
	return;

    subdiv = w->sDiv;
    nRow = (subdiv)?(2 << (subdiv - 1)) + 1 : 2;
    nVer = (nRow * (nRow + 1)) / 2;

    nWVer = pow (2, subdiv + 1) + 2;

    ratioRadiusToSideDist = as->radius*as->ratio/as->sideDistance;

    r = cs->distance / cosf (M_PI / size);
    ang = M_PI / size;
    aStep = 2 * M_PI / size;
    
    wv = w->vertices + (size * nVer);


    for (l = 0; l < size; l++)
    {
	v =   w->vertices + (l * nVer);

	d.v[0] = b.v[0] = sin (ang - aStep) * r;
	d.v[2] = b.v[2] = cos (ang - aStep) * r;

	e.v[0] = c.v[0] = sin (ang) * r;
	e.v[2] = c.v[2] = cos (ang) * r;

	for (i = 0; i < 3; i++)
	{
	    vab[i] = b.v[i] - a.v[i];
	    vab[i] /= nRow - 1.0;
	    vac[i] = c.v[i] - a.v[i];
	    vac[i] /= nRow - 1.0;
	}

	//v[0] = a;


	for (i = 1; i < nRow; i++)
	{
	    br = (i * (i + 1)) / 2;
	    for (k = 0; k < 3; k++)
	    {
		rb[k] = a.v[k] + (i * vab[k]);
		re[k] = a.v[k] + (i * vac[k]);
		ri[k] = re[k] - rb[k];
		ri[k] /= i;
	    }
	    
	    for (j = 0; j <= i; j++)
	    {
		v[br + j].v[0] = (rb[0] + (j * ri[0]));
		v[br + j].v[2] = (rb[2] + (j * ri[2]));

		float th = atan2(v[br + j].v[0], v[br + j].v[2]);
		float factor = progress*(ratioRadiusToSideDist-1)*fabsf(cosf(size*th/2))+1;
		
		v[br + j].v[0] *= factor;
		v[br + j].v[2] *= factor;
	    }
	}
	
	Vertex *lVer = wv + (l * nWVer / 2);
	Vertex *hVer = wv + ((l + size) * nWVer / 2);
	
	/*side walls */
	    nRowS = pow (2, subdiv);

	    for (i = 0; i < 3; i++)
	    {
		vab[i] = c.v[i] - b.v[i];
		vab[i] /= nRowS;
	    }

	    for (i = 0; i <= nRowS; i++)
	    {
		for (k = 0; k < 3; k+=2)
		{
		    lVer[i].v[k] = b.v[k] + (i * vab[k]);
		    hVer[i].v[k] = lVer[i].v[k];
		    
		}
		float th = atan2(lVer[i].v[0], lVer[i].v[2]);
		float factor = progress*(ratioRadiusToSideDist-1)*fabsf(cosf(size*th/2))+1;
		    
		for (k = 0; k < 3; k+=2)
		    lVer[i].v[k] *= factor;
		
		for (k = 0; k < 3; k+=2)
		    hVer[i].v[k] *= factor;


		
		/*lVer[i].n[0] = sinf(ang);
		lVer[i].n[1] = 0;
		lVer[i].n[2] = cosf(ang);
		
		hVer[i].n[0] = lVer[i].n[0];
		hVer[i].n[1] = lVer[i].n[1];
		hVer[i].n[2] = lVer[i].n[2];*/
	    }
	
	
	
	ang += aStep;
    }

static void
deformSphere(CompScreen *s, Water  *w, float progress, float waterBottom)
{
    ATLANTIS_SCREEN (s);
    CUBE_SCREEN (s);

    int          nVer, nWVer, nWIdx, nWVer2, nWIdx2, nRow, nRowS, subdiv;
    Vertex       *v;
    int          i, j, k, l;
    int          br;

    float  ang, r, aStep;
    
    Vertex       *wv;
    
    int bottom = -0.5, size = as->hsize;
 
    float ratioRadiusToSideDist, sphereRadiusFactor, sphereRadiusFactor2;
   
    Vertex a = {{ 0.0, 0.0, 0.0 }};
    Vertex b = {{ 0.0, 0.0, 0.0 }};
    Vertex c = {{ 0.0, 0.0, 0.0 }};
    Vertex d = {{ 0.0, bottom, 0.0 }};
    Vertex e = {{ 0.0, bottom, 0.0 }};
    
    float    vab[3], vac[3], rb[3], re[3], ri[3];

    if (!w)
	return;
    if (w->sDiv < 0)
	return;
    if (!w->vertices)
	return;
    if (w->size!=size)
	return;

    subdiv = w->sDiv;
    nRow = (subdiv)?(2 << (subdiv - 1)) + 1 : 2;
    nVer = (nRow * (nRow + 1)) / 2;

    nWIdx = pow (2, subdiv + 1) * 3;
    nWVer = pow (2, subdiv + 1) + 2;

    nWIdx2 = nWIdx * (nRow -1)*2;
    nWVer2 = nWVer * nRow / 2;

    ratioRadiusToSideDist = as->radius*as->ratio/as->sideDistance;

    sphereRadiusFactor  = as->radius/100000;
    sphereRadiusFactor  = progress*(hypotf(sphereRadiusFactor, 0.5f)/sphereRadiusFactor-1);
    //sphereRadiusFactor2 = sphereRadiusFactor*cosf(waterBottom*PI)+1;
    //sphereRadiusFactor  = sphereRadiusFactor*cosf(w->bh*PI)+1;
    sphereRadiusFactor2 = sphereRadiusFactor*cosf(w->bh*PI)+1;
    
    r = cs->distance / cosf (M_PI / size);
    ang = M_PI / size;
    aStep = 2 * M_PI / size;
    
    wv = w->vertices + (size * nVer);

    if (nWVer2 * size != w->nWVer2 && w->vertices2)
    {
	free (w->vertices2);
	w->vertices2 = NULL;
    }
    if (nWIdx2 * size != w->nWIdx2 && w->indices2)
    {
	free (w->indices2);
	w->indices2 = NULL;
    }

    w->nWVer2 = nWVer2 * size;
    w->nWIdx2 = nWIdx2 * size;
    
    if (!w->vertices2)
    {
	w->vertices2 = calloc (1,sizeof (Vertex) * w->nWVer2);
	if (!w->vertices2)
	    return;
    }

    if (!w->indices2)
    {
	w->indices2 = calloc (1,sizeof (int) * w->nWIdx2);
    	if (!w->indices2)
    	    return;
    }
    

    for (l = 0; l < size; l++)
    {
	v =   w->vertices + (l * nVer);

	d.v[0] = b.v[0] = sin (ang - aStep) * r;
	d.v[2] = b.v[2] = cos (ang - aStep) * r;

	e.v[0] = c.v[0] = sin (ang) * r;
	e.v[2] = c.v[2] = cos (ang) * r;

	for (i = 0; i < 3; i++)
	{
	    vab[i] = b.v[i] - a.v[i];
	    vab[i] /= nRow - 1.0;
	    vac[i] = c.v[i] - a.v[i];
	    vac[i] /= nRow - 1.0;
	}

	//v[0] = a;


	for (i = 1; i < nRow; i++)
	{
	    br = (i * (i + 1)) / 2;
	    for (k = 0; k < 3; k++)
	    {
		rb[k] = a.v[k] + (i * vab[k]);
		re[k] = a.v[k] + (i * vac[k]);
		ri[k] = re[k] - rb[k];
		ri[k] /= i;
	    }
	    
	    for (j = 0; j <= i; j++)
	    {
		v[br + j].v[0] = (rb[0] + (j * ri[0]));
		v[br + j].v[2] = (rb[2] + (j * ri[2]));

		float th = atan2(v[br + j].v[0], v[br + j].v[2]);

		float factor = progress*(ratioRadiusToSideDist-1)*
			       (fabsf(cosf(size*th/2)))+1;
		factor *= sphereRadiusFactor2;
		
		v[br + j].v[0] *= factor;
		v[br + j].v[2] *= factor;
	    }
	}
	
	//Vertex *lVer = w->vertices2 + (l * nWVer / 2);
	//Vertex *hVer = w->vertices2 + ((l + size) * nWVer / 2);

	Vertex *lVer = w->vertices2 + (l * nWVer2 / nRow);

	
	/*side walls */
	    nRowS = pow (2, subdiv);

	    for (i = 0; i < 3; i++)
	    {
		vab[i] = c.v[i] - b.v[i];
		vab[i] /= nRowS;
	    }

	    for (i = 0; i <= nRowS; i++)
	    {
		for (k = 0; k < 3; k++)

		    lVer[i].v[k] = b.v[k] + (i * vab[k]);

		float th = atan2(lVer[i].v[0], lVer[i].v[2]);
		float lFactor = progress*(ratioRadiusToSideDist-1)*
			        (fabsf(cosf(size*th/2)))+1;		
		//lFactor *= sphereRadiusFactor * cosf(w->bh*PI)+1;

		lVer[i].n[0] = (1-progress)*sinf(ang) + progress*sinf(th);
		lVer[i].n[1] = 0;
		lVer[i].n[2] = (1-progress)*cosf(ang) + progress*cosf(th);

		for (j=nRow-1; j>=0; j--)
                {
		    Vertex *hVer = lVer + j * (size * nWVer2 / nRow);
		    
		    for (k = 0; k < 3; k++)
		    {
			hVer[i].v[k] = lVer[i].v[k];
			hVer[i].n[k] = lVer[i].n[k];
		    }
    
		    float hFactor = lFactor * (sphereRadiusFactor * cosf((w->bh-j*(w->bh-waterBottom)/(nRow-1))*PI)+1);
		    
		    for (k = 0; k < 3; k+=2)
			hVer[i].v[k] *= hFactor;
		}

	    }
	    
	    unsigned int * indices = w->indices2 + (l * nWIdx);
	    unsigned int idxBaseL = (l * nWVer / 2);
	    
	    for (j=0; j<nRow-1; j++)
	    {
		unsigned int idxBaseH = idxBaseL + size * nWVer / 2;
		
		for (i = 0; i < nRowS; i++)
		{
		    indices[(i * 6)] = idxBaseL + i;
		    indices[(i * 6) + 1] = idxBaseH + i;
		    indices[(i * 6) + 2] = idxBaseH + i + 1;
		    indices[(i * 6) + 3] = idxBaseL + i + 1;
		    indices[(i * 6) + 4] = idxBaseL + i;
		    indices[(i * 6) + 5] = idxBaseH + i + 1;
		}
		idxBaseL = idxBaseH;
		indices += 2*nWIdx2 / (nRow-1);
	    }
	    

	ang += aStep;
    }
}

void
updateHeight (Water  *w, Water *w2, Bool rippleEffect, int currentDeformation)
{
    int offset;

    Bool useOtherWallVertices;
    Vertex * vertices;
    
    int i, j;
    
    if (!w)
	return;

    offset = w->nSVer/2 + 1;
    rippleEffect = (rippleEffect && w->rippleFactor);

    useOtherWallVertices = (currentDeformation == DeformationSphere &&
	    		    w->vertices2);
    vertices = (useOtherWallVertices ? w->vertices2 - w->nSVer : w->vertices);
    
    for (i = 0; i < w->nSVer; i++)
	setAmplitude(&w->vertices[i], w->bh, w->wave1, w->wave2, w->wa,
	             w->swa, w->wf, w->swf,
	             (rippleEffect ? w->rippleFactor[i] : 0),
	             (rippleEffect ? w->rippleFactor[(i+offset)%w->nSVer] :
	             		     0) );
    
    for (i = w->nSVer; i < w->nSVer + (w->nWVer / 2); i++)
        setAmplitude(&vertices[i], w->bh, w->wave1, w->wave2, w->wa,
		     w->swa, w->wf, w->swf, 0, 0);

    if (useOtherWallVertices)
    {
	int nRow = (w->sDiv)?(2 << (w->sDiv - 1)) + 1 : 2;

	Vertex * verticesL = vertices;
	
	for (j=1; j< nRow-1; j++ )
	{
	    vertices += w->nWVer / 2;//(nRow-1)*w->nWVer2/nRow;//   (w->nWVer / 2);
	    
	    for (i=w->nSVer; i < w->nSVer + (w->nWVer / 2); i++)
		vertices[i].v[1] = verticesL[i].v[1]-j*(verticesL[i].v[1]+0.5)/(nRow-1);
	}
	    
	vertices += w->nWVer / 2;
	
	//FIX get rid of the sphere check
	if (w2 && currentDeformation!=DeformationSphere) /* this is okay because ground and water have same grid size */
	    for (i = w2->nSVer; i < w2->nSVer + (w2->nWVer / 2); i++)
	        setAmplitude(&vertices[i], w2->bh, w2->wave1, w2->wave2, w2->wa,
			     w2->swa, w2->wf, w2->swf, 0, 0);
	else
	    for (i = w->nSVer; i < w->nSVer + (w->nWVer / 2); i++)
	        vertices[i].v[1] = -0.5;
    }
}

void
updateDeformation (CompScreen *s, int currentDeformation)
{
    SNOWGLOBE_SCREEN (s);
    CUBE_SCREEN (s);
    
    static const float floatErr = 0.0001;

    Bool deform = FALSE;

    float progress, dummy;
    (*cs->getRotation) (s, &dummy, &dummy, &progress);

    if (currentDeformation == DeformationNone)
    {
	if (as->oldProgress==0.0f)
	    return;
	
	as->oldProgress = 0.0f;
	progress = 0.0f;
	}
    }
    else
    {
	if (fabsf(progress) < floatErr)
	    progress = 0.0f;
	else if (fabsf(1.0f - progress) < floatErr)
	    progress = 1.0f;
	    
	if ((as->oldProgress!=0.0f || progress!=0.0f) &&
		(as->oldProgress!=1.0f || progress!=1.0f))
	{
	    if (progress==0.0f || progress==1.0f)
	    {
		if (as->oldProgress!=progress)
		{
		    deform = TRUE;
		    as->oldProgress = progress;
		}
	    }
	    else if (fabsf(as->oldProgress-progress)>= floatErr)
	    {
		deform = TRUE;
		as->oldProgress = progress;
	    }
	}
    }

    if (deform)
    {
	if (snowglobeGetShowWater (s) || snowglobeGetShowWaterWire (s))
	{
	    switch (currentDeformation)
	    {
	    case DeformationNone :
	    case DeformationCylinder :
		deformCylinder(s, as->water, progress);
		break;
		
	    case DeformationSphere :
		deformSphere(s, as->water, progress,
		             //snowglobeGetShowGround (s) ? as->ground->bh :
		              -0.5);
	    }    
	}

	if (snowglobeGetShowGround (s))
	{
	    switch (currentDeformation)
	    {
	    case DeformationNone :
		progress = 0.0f;
	    case DeformationCylinder :
		deformCylinder (s, as->ground, progress);
		break;
		
	    case DeformationSphere :
		deformSphere (s, as->ground, progress, -0.5);
	    }

	    updateHeight (as->ground, NULL, FALSE, currentDeformation);
	}
    }
}

void
updateWater (CompScreen *s, float time)
{
    SNOWGLOBE_SCREEN (s);
    CUBE_SCREEN (s);

    int sDiv = (snowglobeGetRenderWaves (s))?
	       snowglobeGetGridQuality (s) : 0;
    int size = as->hsize;

    if (!as->water)
	as->water = genWater (size, sDiv, cs->distance, -0.5, snowglobeGetWaveRipple (s));

    if (!as->water)
	return;

    if (as->water->size != size || as->water->sDiv != sDiv ||
	as->water->distance != cs->distance || (snowglobeGetWaveRipple (s) && !as->water->rippleFactor))
    {
	freeWater (as->water);
	as->water = genWater (size, sDiv, cs->distance, -0.5, snowglobeGetWaveRipple (s));

	if (!as->water)
	    return;
    }

    if (snowglobeGetWaveRipple(s))
    {
	as->water->rippleTimer -= (int) (time * 1000);
	if (as->water->rippleTimer <= 0)
	{
	    as->water->rippleTimer += 100;
	    updateRipple(as->water, size);
	}
    }
    
    as->water->wave1 += time * as->speedFactor;
    as->water->wave2 += time * as->speedFactor;

    as->water->wave1 = fmodf (as->water->wave1, 2 * M_PI);
    as->water->wave2 = fmodf (as->water->wave2, 2 * M_PI);

    if (snowglobeGetRenderWaves (s))
    {
	as->water->wa  = snowglobeGetWaveAmplitude (s);
	as->water->swa = snowglobeGetSmallWaveAmplitude (s);
 	as->water->wf  = snowglobeGetWaveFrequency (s);
	as->water->swf = snowglobeGetSmallWaveFrequency (s);
    }
    else
    {
	as->water->wa  = 0.0;
	as->water->swa = 0.0;
 	as->water->wf  = 0.0;
	as->water->swf = 0.0;
    }

    as->water->bh  = -0.5 + snowglobeGetWaterHeight (s);
}

void
updateGround (CompScreen *s, float time)
{
    SNOWGLOBE_SCREEN (s);
    CUBE_SCREEN (s);

    int sDiv = snowglobeGetGridQuality (s);
    int size = as->hsize;

    Bool update = FALSE;

    if (!as->ground)
    {
	as->ground = genWater (size, sDiv, cs->distance, -0.5, FALSE);
	update = TRUE;
    }

    if (!as->ground)
	return;

    if (as->ground->size != size || as->ground->sDiv != sDiv ||
	as->ground->distance != cs->distance)
    {
	freeWater (as->ground);
	as->ground = genWater (size, sDiv, cs->distance, -0.5, FALSE);

	update = TRUE;
	if (!as->ground)
	    return;
    }

    if (!update)
	return;

    as->ground->wave1 = (float)(rand() & 15) / 15.0;
    as->ground->wave2 = (float)(rand() & 15) / 15.0;

    as->ground->bh  = -0.45;
    as->ground->wa  = 0.1;
    as->ground->swa = 0.02;
    as->ground->wf  = 2.0;
    as->ground->swf = 10.0;

    updateHeight (as->ground, NULL, FALSE, DeformationNone);

}

void
drawWater (Water *w, Bool full, Bool wire, int currentDeformation)
{
    static const float mat_shininess[]      = { 50.0 };
    static const float mat_specular[]       = { 0.5, 0.5, 0.5, 1.0 };
    static const float mat_diffuse[]        = { 0.2, 0.2, 0.2, 1.0 };
    static const float mat_ambient[]        = { 0.1, 0.1, 0.1, 1.0 };
    static const float lmodel_ambient[]     = { 0.4, 0.4, 0.4, 1.0 };
    static const float lmodel_localviewer[] = { 0.0 };

    float *v;
    if (!w)
	return;

    glDisable (GL_DEPTH_TEST);

    if (full)
    {
	glMaterialfv (GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glLightModelfv (GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
	glLightModelfv (GL_LIGHT_MODEL_LOCAL_VIEWER, lmodel_localviewer);

	glEnable (GL_COLOR_MATERIAL);
	glEnable (GL_LIGHTING);
	glEnable (GL_LIGHT1);
	glDisable (GL_LIGHT0);

	glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	v = (float *) w->vertices;
	glDisableClientState (GL_TEXTURE_COORD_ARRAY);
	glEnableClientState (GL_NORMAL_ARRAY);

	glVertexPointer (3, GL_FLOAT, 6 * sizeof (float), v);
	glNormalPointer (GL_FLOAT, 6 * sizeof (float), v + 3);
	glDrawElements (GL_TRIANGLES, w->nSIdx, GL_UNSIGNED_INT, w->indices);

	glDisableClientState (GL_NORMAL_ARRAY);
	glDisable (GL_LIGHTING);

	if (currentDeformation == DeformationSphere && w->vertices2 && w->indices2)
	{
	    v = (float *) w->vertices2;
		
	    glVertexPointer (3, GL_FLOAT, 6 * sizeof (float), v);
	    glNormalPointer (GL_FLOAT, 6 * sizeof (float), v + 3);

	    glDrawElements (GL_TRIANGLES, w->nWIdx2, GL_UNSIGNED_INT, w->indices2);
	}
	else
	    glDrawElements (GL_TRIANGLES, w->nWIdx, GL_UNSIGNED_INT, w->indices + w->nSIdx);
	


	glEnableClientState (GL_TEXTURE_COORD_ARRAY);

    glColor4usv (defaultColor);

    }

    if (wire)
    {
	int i, j;

	glDisable (GL_LIGHTING);

	glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	for (i = 0; i < w->nIndices; i+=3)
	{
	    glBegin(GL_LINE_LOOP);

	    for (j = 0; j < 3; j++)
		glVertex3f(w->vertices[w->indices[i + j]].v[0],
			   w->vertices[w->indices[i + j]].v[1],
			   w->vertices[w->indices[i + j]].v[2]);
	    glEnd();
	}
    }

}

void
drawGround (Water *w, Water *g, int currentDeformation)
{
    static const float mat_shininess[]      = { 40.0 };
    static const float mat_specular[]       = { 0.0, 0.0, 0.0, 1.0 };
    static const float mat_diffuse[]        = { -1.0, -1.0, -1.0, 1.0 };
    static const float mat_ambient[]        = { 0.4, 0.4, 0.4, 1.0 };
    static const float lmodel_ambient[]     = { 0.4, 0.4, 0.4, 1.0 };
    static const float lmodel_localviewer[] = { 0.0 };

    float *v;
    float *n;

    if (!g)
	return;

    glEnable (GL_DEPTH_TEST);

    glMaterialfv (GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv (GL_LIGHT_MODEL_LOCAL_VIEWER, lmodel_localviewer);

    glEnable (GL_COLOR_MATERIAL);
    glEnable (GL_LIGHTING);
    glEnable (GL_LIGHT1);
    glDisable (GL_LIGHT0);

    glTexEnvi (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

    v = (float *) g->vertices;
    glDisableClientState (GL_TEXTURE_COORD_ARRAY);

    glVertexPointer (3, GL_FLOAT, 6 * sizeof (float), v);

    if (w)
    {
	 n = (float *) w->vertices;
	glEnableClientState (GL_NORMAL_ARRAY);
	glNormalPointer (GL_FLOAT, 6 * sizeof (float), n + 3);
    }
    else
	glNormal3f (0.0, 1.0, 0.0);

    glDrawElements (GL_TRIANGLES, g->nSIdx, GL_UNSIGNED_INT, g->indices);

    glDisableClientState (GL_NORMAL_ARRAY);
    glDisable (GL_LIGHTING);

    if (currentDeformation == DeformationSphere && g->vertices2 && g->indices2)
    {
	v = (float *) g->vertices2;

	if (w)
	    n = (float *) w->vertices + 6 * g->nSIdx;
	else
	    n = (float *) g->vertices2;
	
	glNormalPointer (GL_FLOAT, 6 * sizeof (float), n + 3);
	glVertexPointer (3, GL_FLOAT, 6 * sizeof (float), v);

	glDrawElements (GL_TRIANGLES, g->nWIdx2, GL_UNSIGNED_INT, g->indices2);
    }
    else
	glDrawElements (GL_TRIANGLES, g->nWIdx, GL_UNSIGNED_INT, g->indices + g->nSIdx);

    glEnableClientState (GL_TEXTURE_COORD_ARRAY);

}

void
drawBottomGround (int size, float distance, float bottom)
{
	glEnable (GL_COLOR_MATERIAL);

    int   i;
    float r = distance / cos (M_PI / size);
    float ang = M_PI / size;
    float aStep = 2 * M_PI / size;

    for (i = 0; i < size; i++)
    {
	glBegin (GL_TRIANGLES);

	glVertex3f (sin (ang - aStep) * r, bottom, cos (ang - aStep) * r);
	glVertex3f (0.0, bottom, 0.0);
	glVertex3f (sin (ang) * r, bottom, cos (ang) * r);
	glEnd ();
	ang += aStep;
    }
}

float
getHeight (Water *w, float x, float z)
{
    if (!w)
	return 0.0;
    return w->bh + (w->wa * sinf (w->wave1 + w->wf * x * z)) +
	   (w->swa * sinf (w->wave2 + w->swf * x * z));
}
