#pragma once

#include <stdint.h>

#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

namespace kdtree
{

//****
/**
 * Class for representing a point. coordinate_type must be a numeric type.
 */
template<typename coordinate_type, size_t dimensions>
class point 
{
public:
    point(coordinate_type x,coordinate_type y,coordinate_type z)
    {
        coords_[0] = x;
        coords_[1] = y;
        coords_[2] = z;
    }
    point(std::array<coordinate_type, dimensions> c) : coords_(c) 
    {
    }
    point(std::initializer_list<coordinate_type> list) 
    {
        size_t n = std::min(dimensions, list.size());
        std::copy_n(list.begin(), n, coords_.begin());
    }

    /**
     * Returns the coordinate in the given dimension.
     *
     * @param index dimension index (zero based)
     * @return coordinate in the given dimension
     */
    coordinate_type get(size_t index) const 
    {
        return coords_[index];
    }

    /**
     * Returns the distance squared from this point to another
     * point.
     *
     * @param pt another point
     * @return distance squared from this point to the other point
     */
    float distance(const point& pt) const 
    {
        float dist = 0;
        for (size_t i = 0; i < dimensions; ++i) 
        {
            float d = get(i) - pt.get(i);
            dist += d * d;
        }
        return dist;
    }

    void setId(uint32_t id) 
    {
        mId = id;
    }
    uint32_t getId(void) const
    {
        return mId;
    }
private:
    std::array<coordinate_type, dimensions> coords_;
    uint32_t    mId;
};

/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */
template<typename coordinate_type, size_t dimensions>
class RosettaKdTree 
{
public:
    typedef point<coordinate_type, dimensions> point_type;
private:
    struct node 
    {
        node(const point_type& pt) : point_(pt), left_(nullptr), right_(nullptr) 
        {
        }
        coordinate_type get(size_t index) const 
        {
            return point_.get(index);
        }
        float distance(const point_type& pt) const 
        {
            return point_.distance(pt);
        }
        point_type point_;
        node* left_;
        node* right_;
    };
    node* root_ = nullptr;
    node* best_ = nullptr;
    float best_dist_ = 0;
    size_t visited_ = 0;
    std::vector<node> nodes_;

    struct node_cmp 
    {
        node_cmp(size_t index) : index_(index) 
        {
        }

        bool operator()(const node& n1, const node& n2) const 
        {
            return n1.point_.get(index_) < n2.point_.get(index_);
        }

        size_t index_;
    };

    node* make_tree(size_t begin, size_t end, size_t index) 
    {
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

    void nearest(node* root, const point_type& point, size_t index) 
    {
        if (root == nullptr)
            return;
        ++visited_;
        float d = root->distance(point);
        if (best_ == nullptr || d < best_dist_) 
        {
            best_dist_ = d;
            best_ = root;
        }
        if (best_dist_ == 0)
            return;
        float dx = root->get(index) - point.get(index);
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
    RosettaKdTree(iterator begin, iterator end) : nodes_(begin, end) 
    {
        root_ = make_tree(0, nodes_.size(), 0);
    }
    

    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const 
    { 
        return nodes_.empty(); 
    }

    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const 
    { return visited_; 
    }

    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    float distance() const 
    { 
        return std::sqrt(best_dist_); 
    }

    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& nearest(const point_type& pt) 
    {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0);
        return best_->point_;
    }
};

typedef point<float, 3> point3d;
typedef RosettaKdTree<float, 3> tree3d;
//****

class KdPoint
{
public:
    float   x;
    float   y;
    float   z;
    uint32_t mId;
};

using Point3dVector = std::vector< point3d >;

class KdTree
{
public:
    KdTree(void)
    {
    }

    ~KdTree(void)
    {
        delete mTree;
    }

    void reservePoints(uint32_t pcount)
    {
        delete mTree;
        mTree = nullptr;
        mPoints.clear();
        mPoints.reserve(pcount);
    }

    // Add this point...
    void addPoint(const KdPoint &p)
    {
        point3d pp(p.x,p.y,p.z);
        pp.setId(p.mId);
        mPoints.push_back(pp);
    }

    void buildTree(void)
    {
        mTree = new tree3d(std::begin(mPoints),std::end(mPoints));
    }

    float findNearest(const KdPoint &p,KdPoint &result)
    {
        float ret = -1;

        if ( mTree )
        {
            point3d pt(p.x,p.y,p.z);
            point3d n = mTree->nearest(pt);
            result.x = n.get(0);
            result.y = n.get(1);
            result.z = n.get(2);
            result.mId = n.getId();
            ret = n.distance(pt);
        }

        return ret;
    }

private:
    tree3d          *mTree{nullptr};
    Point3dVector   mPoints;
};


}

