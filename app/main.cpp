#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <random>
#include <assert.h>
#include <algorithm>
#include <array>
#include <vector>

#include "ScopedTime.h"
#include "kdtree.h"

//
/**
 * Class for representing a point. coordinate_type must be a numeric type.
 */
template<typename coordinate_type, size_t dimensions>
class point {
public:
    point(coordinate_type x,coordinate_type y,coordinate_type z)
    {
        coords_[0] = x;
        coords_[1] = y;
        coords_[2] = z;
    }
    point(std::array<coordinate_type, dimensions> c) : coords_(c) {}
    point(std::initializer_list<coordinate_type> list) {
        size_t n = std::min(dimensions, list.size());
        std::copy_n(list.begin(), n, coords_.begin());
    }
    /**
     * Returns the coordinate in the given dimension.
     *
     * @param index dimension index (zero based)
     * @return coordinate in the given dimension
     */
    coordinate_type get(size_t index) const {
        return coords_[index];
    }
    /**
     * Returns the distance squared from this point to another
     * point.
     *
     * @param pt another point
     * @return distance squared from this point to the other point
     */
    double distance(const point& pt) const {
        double dist = 0;
        for (size_t i = 0; i < dimensions; ++i) {
            double d = get(i) - pt.get(i);
            dist += d * d;
        }
        return dist;
    }
private:
    std::array<coordinate_type, dimensions> coords_;
};

template<typename coordinate_type, size_t dimensions>
std::ostream& operator<<(std::ostream& out, const point<coordinate_type, dimensions>& pt) {
    out << '(';
    for (size_t i = 0; i < dimensions; ++i) {
        if (i > 0)
            out << ", ";
        out << pt.get(i);
    }
    out << ')';
    return out;
}

/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */
template<typename coordinate_type, size_t dimensions>
class RosettaKdTree {
public:
    typedef point<coordinate_type, dimensions> point_type;
private:
    struct node {
        node(const point_type& pt) : point_(pt), left_(nullptr), right_(nullptr) {}
        coordinate_type get(size_t index) const {
            return point_.get(index);
        }
        double distance(const point_type& pt) const {
            return point_.distance(pt);
        }
        point_type point_;
        node* left_;
        node* right_;
    };
    node* root_ = nullptr;
    node* best_ = nullptr;
    double best_dist_ = 0;
    size_t visited_ = 0;
    std::vector<node> nodes_;

    struct node_cmp {
        node_cmp(size_t index) : index_(index) {}
        bool operator()(const node& n1, const node& n2) const {
            return n1.point_.get(index_) < n2.point_.get(index_);
        }
        size_t index_;
    };

    node* make_tree(size_t begin, size_t end, size_t index) {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin)/2;
        auto i = nodes_.begin();
        std::nth_element(i + begin, i + n, i + end, node_cmp(index));
        index = (index + 1) % dimensions;
        nodes_[n].left_ = make_tree(begin, n, index);
        nodes_[n].right_ = make_tree(n + 1, end, index);
        return &nodes_[n];
    }

    void nearest(node* root, const point_type& point, size_t index) {
        if (root == nullptr)
            return;
        ++visited_;
        double d = root->distance(point);
        if (best_ == nullptr || d < best_dist_) {
            best_dist_ = d;
            best_ = root;
        }
        if (best_dist_ == 0)
            return;
        double dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        nearest(dx > 0 ? root->left_ : root->right_, point, index);
        if (dx * dx >= best_dist_)
            return;
        nearest(dx > 0 ? root->right_ : root->left_, point, index);
    }
public:
    RosettaKdTree(const RosettaKdTree&) = delete;
    RosettaKdTree& operator=(const RosettaKdTree&) = delete;
    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template<typename iterator>
    RosettaKdTree(iterator begin, iterator end) : nodes_(begin, end) {
        root_ = make_tree(0, nodes_.size(), 0);
    }
    
    /**
     * Constructor taking a function object that generates
     * points. The function object will be called n times
     * to populate the tree.
     *
     * @param f function that returns a point
     * @param n number of points to add
     */
    template<typename func>
    RosettaKdTree(func&& f, size_t n) {
        nodes_.reserve(n);
        for (size_t i = 0; i < n; ++i)
            nodes_.push_back(f());
        root_ = make_tree(0, nodes_.size(), 0);
    }

    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const { return nodes_.empty(); }

    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const { return visited_; }

    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double distance() const { return std::sqrt(best_dist_); }

    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& nearest(const point_type& pt) {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0);
        return best_->point_;
    }
};
//

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

typedef point<double, 3> point3d;
typedef RosettaKdTree<double, 3> tree3d;

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
    Point3 nearestMatch;
    for (uint32_t i=0; i<pcount; i++)
    {
        float d2 = fp.distanceSquared(points[i]);
        if ( d2 < closest )
        {
            closest = d2;
            index = i;
            nearestMatch = points[i];
        }
    }

    {
        std::vector< point3d > kpoints;
        for (uint32_t i=0; i<pcount; i++)
        {
            point3d p(points[i].x,points[i].y,points[i].z);
            kpoints.push_back(p);
        }
        tree3d tree(std::begin(kpoints),std::end(kpoints));
        point3d pf(fp.x,fp.y,fp.z);
        point3d n = tree.nearest(pf);
        printf("Result\n");
    }

    kdtree::KdPoint kp,result;
    kp.mPos[0] = fp.x;
    kp.mPos[1] = fp.y;
    kp.mPos[2] = fp.z;
    kdt.findNearest(kp,result);
    assert( result.mId == index );


    return 0;
}
