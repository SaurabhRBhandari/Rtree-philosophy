#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

struct ClusTreeNode;
// Structure to represent a micro-cluster

struct MicroCluster
{
    double weight;
    std::vector<double> sumOfAttributes;
    std::vector<double> sumOfSquares;
    double radius;

    // Constructor
    MicroCluster(int dimension) : weight(0.0), radius(0.0)
    {
        sumOfAttributes.resize(dimension, 0.0);
        sumOfSquares.resize(dimension, 0.0);
    }
};

// Function to calculate Euclidean distance between two vectors
double euclideanDistance(const std::vector<double> &v1, const std::vector<double> &v2)
{
    double distance = 0.0;
    for (size_t i = 0; i < v1.size(); ++i)
    {
        distance += std::pow(v1[i] - v2[i], 2);
    }
    return std::sqrt(distance);
}

// Structure to represent an entry in the ClusTree node
struct ClusTreeEntry
{
    MicroCluster clusterFeature;
    MicroCluster bufferFeature; // Buffer feature of the objects in the node
    ClusTreeNode *child;
    int id;

    // Constructor
    ClusTreeEntry(const MicroCluster &cf, const MicroCluster &bf, ClusTreeNode *ptr) : clusterFeature(cf), bufferFeature(bf), child(ptr) {}
};

bool operator==(const ClusTreeEntry &lhs, const ClusTreeEntry &rhs)
{
    return lhs.clusterFeature.sumOfAttributes == rhs.clusterFeature.sumOfAttributes &&
           lhs.bufferFeature.sumOfAttributes == rhs.bufferFeature.sumOfAttributes &&
           lhs.id == rhs.id;
}

// Structure to represent a ClusTree node
struct ClusTreeNode
{
    bool isLeaf;
    std::vector<ClusTreeEntry> entries;

    // Constructor
    ClusTreeNode(bool leaf) : isLeaf(leaf) {}
};

// Function to initialize centroids using kMeans++ algorithm
std::vector<MicroCluster> initializeCentroids(const std::vector<ClusTreeEntry> &entries, int k)
{
    std::vector<MicroCluster> centroids;
    centroids.reserve(k);

    // Select k random entries as initial centroids
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, entries.size() - 1);
    for (int i = 0; i < k; ++i)
    {
        centroids.push_back(entries[dis(gen)].clusterFeature);
    }

    return centroids;
}

std::vector<std::vector<ClusTreeEntry>> kMeansClustering(const std::vector<ClusTreeEntry> &entries, int k, int maxIterations)
{
    std::vector<MicroCluster> centroids = initializeCentroids(entries, k);
    std::vector<std::vector<ClusTreeEntry>> clusters(k);

    // Perform iterations
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        // Assign each entry to the nearest centroid
        for (const auto &entry : entries)
        {
            double minDist = std::numeric_limits<double>::max();
            int nearestCentroid = -1;
            for (size_t i = 0; i < centroids.size(); ++i)
            {
                double dist = euclideanDistance(entry.clusterFeature.sumOfAttributes, centroids[i].sumOfAttributes);
                if (dist < minDist)
                {
                    minDist = dist;
                    nearestCentroid = i;
                }
            }
            clusters[nearestCentroid].push_back(entry);
        }

        // Update centroids based on the mean of assigned entries
        for (size_t i = 0; i < centroids.size(); ++i)
        {
            if (!clusters[i].empty())
            {
                MicroCluster newCentroid(centroids[i].sumOfAttributes.size());
                for (const auto &entry : clusters[i])
                {
                    for (size_t j = 0; j < newCentroid.sumOfAttributes.size(); ++j)
                    {
                        newCentroid.sumOfAttributes[j] += entry.clusterFeature.sumOfAttributes[j];
                    }
                }
                for (size_t j = 0; j < newCentroid.sumOfAttributes.size(); ++j)
                {
                    newCentroid.sumOfAttributes[j] /= clusters[i].size();
                }
                centroids[i] = newCentroid;
            }
        }
    }

    return clusters;
}

class ClusTree
{
private:
    ClusTreeNode *root;
    int dimension;                              // Dimensionality of data points
    int m, M, l, L;                             // Fanout and leaf node capacity parameters
    std::vector<std::pair<int, int>> mergeList; // Merging list to track merges

    // Function to split an inner node if it exceeds the maximum capacity
    void splitInnerNode(ClusTreeNode *node)
    {
        // Check if the node exceeds the maximum capacity
        if (node->entries.size() > M)
        {
            // Create two new inner nodes
            ClusTreeNode *leftNode = new ClusTreeNode(false);
            ClusTreeNode *rightNode = new ClusTreeNode(false);

            // Perform k-means clustering to split entries into two groups
            std::vector<std::vector<ClusTreeEntry>> clusters = kMeansClustering(node->entries, 2, 10);

            // Assign entries to the left or right node based on the cluster they belong to
            for (const auto &cluster : clusters)
            {
                for (const auto &entry : cluster)
                {
                    if (cluster == clusters[0])
                    {
                        leftNode->entries.push_back(entry);
                    }
                    else
                    {
                        rightNode->entries.push_back(entry);
                    }
                }
            }

            // Clear entries in the original node
            node->entries.clear();

            // Add pointers to child nodes
            node->entries.emplace_back(ClusTreeEntry(MicroCluster(0), MicroCluster(0), leftNode));
            node->entries.emplace_back(ClusTreeEntry(MicroCluster(0), MicroCluster(0), rightNode));
        }
    }

    // Function to calculate pairwise distance between two MicroClusters
    double distance(const MicroCluster &mc1, const MicroCluster &mc2)
    {
        double minDistance = std::numeric_limits<double>::max();
        for (size_t i = 0; i < mc1.sumOfAttributes.size(); ++i)
        {
            double dist = euclideanDistance(mc1.sumOfAttributes, mc2.sumOfAttributes);
            if (dist < minDistance)
            {
                minDistance = dist;
            }
        }
        return minDistance;
    }
    // Function to merge the closest two entries in a leaf node
    void mergeClosestEntries(ClusTreeNode *node)
    {
        // Initialize variables to keep track of the closest entries and their distances
        size_t closestEntry1 = 0, closestEntry2 = 1;
        double minDistance = std::numeric_limits<double>::max();

        // Iterate over all pairs of entries and find the closest pair
        for (size_t i = 0; i < node->entries.size(); ++i)
        {
            for (size_t j = i + 1; j < node->entries.size(); ++j)
            {
                double dist = distance(node->entries[i].clusterFeature, node->entries[j].clusterFeature);
                if (dist < minDistance)
                {
                    minDistance = dist;
                    closestEntry1 = i;
                    closestEntry2 = j;
                }
            }
        }

        // Merge the closest two entries by updating their cluster features and removing one entry
        // For simplicity, here we assume merging involves updating the cluster features only
        for (size_t i = 0; i < dimension; ++i)
        {
            node->entries[closestEntry1].clusterFeature.sumOfAttributes[i] += node->entries[closestEntry2].clusterFeature.sumOfAttributes[i];
            node->entries[closestEntry1].clusterFeature.sumOfSquares[i] += node->entries[closestEntry2].clusterFeature.sumOfSquares[i];
        }
        node->entries[closestEntry1].clusterFeature.weight += node->entries[closestEntry2].clusterFeature.weight;
        // TODO:Recalculate the radius if needed
        // ...

        // Remove the second closest entry
        node->entries.erase(node->entries.begin() + closestEntry2);
    }

    // Function to split a leaf node based on pairwise distances between entries
    void splitBasedOnDistances(ClusTreeNode *node)
    {
        std::sort(node->entries.begin(), node->entries.end(), [&](const ClusTreeEntry &a, const ClusTreeEntry &b)
                  {
            // You may need to implement custom comparison logic based on your specific criterion
            // For example, compare sum of distances to other entries
            return distance(a.clusterFeature, b.clusterFeature) < distance(b.clusterFeature, a.clusterFeature); });

        // Split the node into two groups
        // Determine the splitting point (e.g., halfway through the entries)
        size_t splitPoint = node->entries.size() / 2;

        // Create two new leaf nodes and distribute entries between them
        ClusTreeNode *leftNode = new ClusTreeNode(true);
        ClusTreeNode *rightNode = new ClusTreeNode(true);

        for (size_t i = 0; i < splitPoint; ++i)
        {
            leftNode->entries.push_back(node->entries[i]);
        }

        for (size_t i = splitPoint; i < node->entries.size(); ++i)
        {
            rightNode->entries.push_back(node->entries[i]);
        }

        // Clear entries in the original node and set it to an inner node
        node->entries.clear();
        node->isLeaf = false;

        // Add pointers to child nodes
        node->entries.emplace_back(ClusTreeEntry(MicroCluster(0), MicroCluster(0), leftNode));
        node->entries.emplace_back(ClusTreeEntry(MicroCluster(0), MicroCluster(0), rightNode));
    }

    // Function to split a leaf node if it exceeds the maximum capacity
    void splitLeafNode(ClusTreeNode *node)
    {
        // Check if there is time for a split
        // If no time left, merge closest two entries
        // Else, split based on pairwise distances between entries
        bool hasTimeForSplit = true; // Example: Assume there is always time for a split

        if (!hasTimeForSplit)
        {
            // If no time left, merge closest two entries
            mergeClosestEntries(node);
        }
        else
        {
            // Else, split based on pairwise distances between entries
            splitBasedOnDistances(node);
        }
    }

    // Function to initialize the ClusTree using k-means clustering
    std::vector<MicroCluster> initializeWithKMeansutility(const std::vector<std::vector<double>> &data, int k)
    {
        // Perform k-means clustering on the data to get initial cluster centers
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(0, data.size() - 1);

        std::vector<MicroCluster> initialClusters;
        std::vector<size_t> indices(data.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), gen);

        for (int i = 0; i < k; ++i)
        {
            MicroCluster cluster(dimension);
            for (int j = 0; j < dimension; ++j)
            {
                cluster.sumOfAttributes[j] = data[indices[i]][j];
            }
            cluster.weight = 1.0;
            initialClusters.push_back(cluster);
        }

        return initialClusters;
    }

public:
    // Constructor
    ClusTree(int dim, int minFanout, int maxFanout, int minLeafCapacity, int maxLeafCapacity) : root(nullptr), dimension(dim), m(minFanout), M(maxFanout), l(minLeafCapacity), L(maxLeafCapacity) {}

    // Function to initialize the ClusTree using k-means clustering
    void initializeWithKMeans(const std::vector<std::vector<double>> &data, int k)
    {
        std::vector<MicroCluster> initialClusters = initializeWithKMeansutility(data, k);
        for (const auto &cluster : initialClusters)
        {
            insert(cluster, cluster);
        }
    }

    // Function to insert a micro-cluster into the ClusTree
    void insert(const MicroCluster &mc, const MicroCluster &buffer)
    {
        if (root == nullptr)
        {
            root = new ClusTreeNode(true); // Create root as leaf node
            root->entries.emplace_back(mc, buffer, nullptr);
        }
        else
        {
            insert(root, mc, buffer, 0);
        }
    }

    // Function to insert a micro-cluster into the ClusTree
    void insert(ClusTreeNode *node, const MicroCluster &mc, const MicroCluster &buffer, int level)
    {
        if (node->isLeaf)
        {
            // Insert into leaf node
            if (node->entries.size() < L)
            {
                // Leaf node has capacity to accommodate new entry
                node->entries.emplace_back(mc, buffer, nullptr);
            }
            else
            {
                // Leaf node is full, split leaf node
                splitLeafNode(node);

                // Reinsert the micro-cluster into appropriate child node
                double minDistance = std::numeric_limits<double>::max();
                int minChildIndex = -1;
                for (size_t i = 0; i < node->entries.size(); ++i)
                {
                    double dist = distance(node->entries[i].clusterFeature, mc);
                    if (dist < minDistance)
                    {
                        minDistance = dist;
                        minChildIndex = i;
                    }
                }
                insert(node->entries[minChildIndex].child, mc, buffer, level + 1);
            }
        }
        else
        {
            // Insert into inner node
            if (node->entries.size() < M)
            {
                // Inner node has capacity to accommodate new entry
                // Update cluster features and buffer in the entry
                ClusTreeEntry newEntry(mc, buffer, nullptr); // Create a new entry
                // Update cluster features and buffer in the new entry
                // ...
                node->entries.push_back(newEntry); // Add the new entry to the node
            }
            else
            {
                // Inner node is full, split inner node
                splitInnerNode(node);

                // Reinsert the micro-cluster into appropriate child node
                double minDistance = std::numeric_limits<double>::max();
                int minChildIndex = -1;
                for (size_t i = 0; i < node->entries.size(); ++i)
                {
                    double dist = distance(node->entries[i].clusterFeature, mc);
                    if (dist < minDistance)
                    {
                        minDistance = dist;
                        minChildIndex = i;
                    }
                }
                insert(node->entries[minChildIndex].child, mc, buffer, level + 1);
            }
        }
    }
};

int main()
{
    // Example usage
    int dimension = 3;       // Example dimensionality of data points
    int minFanout = 2;       // Example minimum fanout parameter
    int maxFanout = 4;       // Example maximum fanout parameter
    int minLeafCapacity = 2; // Example minimum leaf node capacity parameter
    int maxLeafCapacity = 4; // Example maximum leaf node capacity parameter
    ClusTree clusTree(dimension, minFanout, maxFanout, minLeafCapacity, maxLeafCapacity);

    // Example: Generating random data for k-means clustering
    std::vector<std::vector<double>> data = {
        {0.1, 0.2, 0.3},
        {0.4, 0.5, 0.6},
        {0.7, 0.8, 0.9},
        // Add more data points as needed
    };

    // Initialize ClusTree using k-means clustering
    int k = 2; // Number of clusters
    clusTree.initializeWithKMeans(data, k);
    printf("ClusTree initialized with k-means clustering\n");
    return 0;
}
