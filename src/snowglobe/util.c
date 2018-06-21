//a collection of utility methods

#include "snowglobe-internal.h"
#include <math.h>
#include <float.h>



float randf (float x) { //return random number in range [0,x)
	return rand()/(((double)RAND_MAX + 1) / x);
}
float minimum (float x, float y) {
	return ((x) < (y) ? (x) : (y));
}
float maximum (float x, float y) {
	return ((x) > (y) ? (x) : (y));
}

int
getCurrentDeformation (CompScreen *s)
{
    CUBE_SCREEN (s);

    CompPlugin *p = NULL;
    const char plugin[] = "cubeaddon";
    p = findActivePlugin (plugin);
    if (p && p->vTable->getObjectOptions)
    {
	CompOption *option;
	int  nOption;
	Bool cylinderManualOnly = FALSE;
	Bool unfoldDeformation = TRUE;

	option = (*p->vTable->getObjectOptions) (p, (CompObject *)s,
		&nOption);
	option = compFindOption (option, nOption, "cylinder_manual_only", 0);

	if (option)
	    if (option->value.b)
		cylinderManualOnly = option->value.b;

	option = (*p->vTable->getObjectOptions) (p, (CompObject *)s,
		&nOption);
	option = compFindOption (option, nOption, "unfold_deformation", 0);

	if (option)
	    if (option->value.b)
		unfoldDeformation = option->value.b;

	if (s->hsize * cs->nOutput > 2 && s->desktopWindowCount &&
	    (cs->rotationState == RotationManual ||
	    (cs->rotationState == RotationChange &&
	    !cylinderManualOnly)) &&
	    (!cs->unfolded || unfoldDeformation))
	{
	    option = (*p->vTable->getObjectOptions) (p, (CompObject *)s,
		      &nOption);
	    option = compFindOption (option, nOption, "deformation", 0);

	    if (option)
		return (option->value.i);
	}
    }
    return DeformationNone;
}

int
getDeformationMode (CompScreen *s)
{
    CompPlugin *p = NULL;
    const char plugin[] = "cubeaddon";
    p = findActivePlugin (plugin);
    if (p && p->vTable->getObjectOptions)
    {
	CompOption *option;
	int  nOption;
	option = (*p->vTable->getObjectOptions) (p, (CompObject *)s,
		  &nOption);
	option = compFindOption (option, nOption, "deformation", 0);

	if (option)
	    return (option->value.i);
    }
    return DeformationNone;
}

float symmDistr() { //returns number in range [-1, 1] with bias towards 0, symmetric about 0.
	float x = 2*randf(1)-1;
	return x*(1-cbrt(1-fabsf(x)));
}
