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
        mCoordinates[0] = x;
        mCoordinates[1] = y;
        mCoordinates[2] = z;
    }
    point(std::array<coordinate_type, dimensions> c) : mCoordinates(c) 
    {
    }

    /**
     * Returns the coordinate in the given dimension.
     *
     * @param index dimension index (zero based)
     * @return coordinate in the given dimension
     */
    coordinate_type get(size_t index) const 
    {
        return mCoordinates[index];
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
    std::array<coordinate_type, dimensions> mCoordinates;
    uint32_t    mId;
};

// Templated implementation of KdTree
template<typename coordinate_type, size_t dimensions>
class KdTreeTemplate 
{
public:
    typedef point<coordinate_type, dimensions> point_type;
private:
    struct KdNode 
    {
        KdNode(const point_type& pt) : mPoint(pt), mLeft(nullptr), mRight(nullptr) 
        {
        }
        coordinate_type get(size_t index) const 
        {
            return mPoint.get(index);
        }

        float distance(const point_type& pt) const 
        {
            return mPoint.distance(pt);
        }

        point_type mPoint;
        KdNode* mLeft;
        KdNode* mRight;
    };
    KdNode              *mRoot{nullptr};
    KdNode              *mBest ={nullptr};
    float               mBestDistance{0};
    size_t              mVisitCount{0};
    std::vector<KdNode> mNodes;

    class KdNodeCompare 
    {
    public:
        KdNodeCompare(size_t index) : index_(index) 
        {
        }

        bool operator()(const KdNode& n1, const KdNode& n2) const 
        {
            return n1.mPoint.get(index_) < n2.mPoint.get(index_);
        }

        size_t index_;
    };

    KdNode* buildKdTree(size_t begin, size_t end, size_t index) 
    {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin)/2;
        auto i = mNodes.begin();
        std::nth_element(i + begin, i + n, i + end, KdNodeCompare(index));
        index = (index + 1) % dimensions;
        mNodes[n].mLeft  = buildKdTree(begin, n, index);
        mNodes[n].mRight = buildKdTree(n + 1, end, index);
        return &mNodes[n];
    }

    void nearest(KdNode* root, const point_type& point, size_t index) 
    {
        if (root == nullptr)
            return;
        ++mVisitCount;
        float d = root->distance(point);
        if (mBest == nullptr || d < mBestDistance) 
        {
            mBestDistance = d;
            mBest = root;
        }
        if (mBestDistance == 0)
            return;
        float dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        nearest(dx > 0 ? root->mLeft : root->mRight, point, index);
        if (dx * dx >= mBestDistance)
            return;
        nearest(dx > 0 ? root->mRight : root->mLeft, point, index);
    }
public:
    KdTreeTemplate(const KdTreeTemplate&) = delete;
    KdTreeTemplate& operator=(const KdTreeTemplate&) = delete;
    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template<typename iterator>
    KdTreeTemplate(iterator begin, iterator end) : mNodes(begin, end) 
    {
        mRoot = buildKdTree(0, mNodes.size(), 0);
    }
    

    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const 
    { 
        return mNodes.empty(); 
    }

    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const 
    { return mVisitCount; 
    }

    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    float distance() const 
    { 
        return std::sqrt(mBestDistance); 
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
        if (mRoot == nullptr)
            throw std::logic_error("tree is empty");
        mBest = nullptr;
        mVisitCount = 0;
        mBestDistance = 0;
        nearest(mRoot, pt, 0);
        return mBest->mPoint;
    }
};

typedef point<float, 3> point3d;
typedef KdTreeTemplate<float, 3> tree3d;

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

