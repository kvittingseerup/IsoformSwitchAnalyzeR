#include "R.h"
#include "R_ext/Rdynload.h"
#include "Rdefines.h"

void findIdenticalExons(int *exon1, int *exon2, int *identical)
{
	if ((exon1[0] == exon2[0]) & (exon1[1] == exon2[1]))
		identical[0] = 1;
	else
		identical[0] = 0;
}

void findOverlappingExons(int *exon1, int *exon2, int *identical)
{

	identical[0] = 0;
	if ((exon1[0] <  exon2[1]) & (exon1[0] >= exon2[0])) identical[0] = 1;
	if ((exon1[1] <= exon2[1]) & (exon1[1] >  exon2[0])) identical[0] = 1;
	if ((exon1[0] <  exon2[0]) & (exon1[1] >  exon2[1])) identical[0] = 1;

}

void determineLeftOverlappingAStype(int *exon1, int *exon2, int *isStrandEqualToPlus, int *ExonIndexesAnalyzed, int *asTypes)
{
	asTypes[0] = 0; asTypes[1] = 0;
	int isExonIndexEqualTo11[2] = {0,0};
	if ( ExonIndexesAnalyzed[0] == 1 ) isExonIndexEqualTo11[0] = 1;
	if ( ExonIndexesAnalyzed[1] == 1 ) isExonIndexEqualTo11[1] = 1;
	if ( isStrandEqualToPlus[0] )
	{
		if( exon1[0] != exon2[0] )
		{
			if ( isExonIndexEqualTo11[0] | isExonIndexEqualTo11[1] )
			{
				return;
			} 
			else {
				asTypes[1] = 1;
			}
		}
	} else
	{
		if( exon1[0] != exon2[0] )
		{
			if ( isExonIndexEqualTo11[0] | isExonIndexEqualTo11[1] )
			{
				return;
			} 
			else {
				asTypes[0] = 1;
			}
		}
	}
}

void determineRightOverlappingAStype(int *exon1, int *exon2, int *isStrandEqualToPlus, int *ExonIndexesAnalyzed, int *numberOfExons, int *asTypes)
{
	asTypes[0] = 0; asTypes[1] = 0;
	int isExonIndexEqualToEnd[2] = {0,0};
	if ( ExonIndexesAnalyzed[0] == numberOfExons[0] ) isExonIndexEqualToEnd[0] = 1;
	if ( ExonIndexesAnalyzed[1] == numberOfExons[1] ) isExonIndexEqualToEnd[1] = 1;
	if ( isStrandEqualToPlus[0] )
	{
		if( exon1[1] != exon2[1] )
		{
			if ( isExonIndexEqualToEnd[0] | isExonIndexEqualToEnd[1] )
			{
				return;
			} 
			else {
				asTypes[0] = 1;
			}
		}
	} else
	{
		if( exon1[1] != exon2[1] )
		{
			if ( isExonIndexEqualToEnd[0] | isExonIndexEqualToEnd[1] )
			{
				return;
			} 
			else {
				asTypes[1] = 1;
			}
		}
	}
}
