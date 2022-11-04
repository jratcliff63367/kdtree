#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <random>
#include <assert.h>

#include "ScopedTime.h"
#include "kdtree.h"

class Point3
{
public:
    float distanceSquared(const Point3 &p) const
    {
        float dx = x-p.x;
        float dy = y-p.y;
        float dz = z-p.z;
        return dx*dx + dy*dy + dz*dz;
    }
    float x,y,z;
};

int main(int argc,const char **argv)
{
    FILE *fph = fopen("kdtree.bin","rb");
    if ( !fph )
    {
        printf("Failed to open 'kdtree.bin'\n");
        exit(1);
    }
    uint32_t pcount=0;
    fread(&pcount,sizeof(pcount),1,fph);
    kdtree::KdTree kdt;
    std::vector< Point3 > points;
    for (uint32_t i=0; i<pcount; i++)
    {
        Point3 p;
        fread(&p,sizeof(p),1,fph);
        points.push_back(p);
        kdtree::KdPoint kp;
        kp.mId = i;
        kp.mPos[0] = p.x;
        kp.mPos[1] = p.y;
        kp.mPos[2] = p.z;
        kdt.addPoint(kp);
    }
    kdt.buildTree();
    Point3 fp;
    fread(&fp,sizeof(fp),1,fph);
    fclose(fph);
    float closest = FLT_MAX;
    uint32_t index = 0;
    for (uint32_t i=0; i<pcount; i++)
    {
        float d2 = fp.distanceSquared(points[i]);
        if ( d2 < closest )
        {
            closest = d2;
            index = i;
        }
    }
    kdtree::KdPoint kp,result;
    kp.mPos[0] = fp.x;
    kp.mPos[1] = fp.y;
    kp.mPos[2] = fp.z;
    kdt.findNearest(kp,result);
    assert( result.mId == index );


    return 0;
}
