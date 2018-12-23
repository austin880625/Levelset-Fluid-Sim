#include "common.h"

int min(const float& a, const int& b)
{
	if(a < b) return a;
	else return b;
}
int max(const float& a, const int& b)
{
	if(!(a < b)) return a;
	else return b;
}
