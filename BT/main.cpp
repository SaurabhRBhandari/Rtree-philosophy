#include <bits/stdc++.h>
using namespace std;
#include "utils.h"
#include "RStarTree.h"

// shortcuts
typedef RStarBoundedItem<2> BoundedItem;
typedef typename BoundedItem::BoundingBox BoundingBox;

typedef RStarNode<BoundedItem> Node;
typedef RStarLeaf<BoundedItem, int> Leaf;

// acceptors
typedef RStarAcceptOverlapping<Node, Leaf> AcceptOverlapping;
typedef RStarAcceptEnclosing<Node, Leaf> AcceptEnclosing;
typedef RStarAcceptAny<Node, Leaf> AcceptAny;

// predefined visitors
typedef RStarRemoveLeaf<Leaf> RemoveLeaf;
typedef RStarRemoveSpecificLeaf<Leaf> RemoveSpecificLeaf;

// TODO: Some points are getting inserted twice

typedef RStarTree<int, 2, 2, 4> RTree;
typedef RTree::BoundingBox BoundingBox;

BoundingBox bounds(int x, int y, int w, int h)
{
	BoundingBox bb;

	bb.edges[0].first = x;
	bb.edges[0].second = x + w;

	bb.edges[1].first = y;
	bb.edges[1].second = y + h;

	bb.n = 1;
	bb.sum[0] = (x + x + w) / 2.0;
	bb.sum[1] = (y + y + h) / 2.0;

	bb.sum_sq[0] = (x * x + (x + w) * (x + w)) / 2.0;
	bb.sum_sq[1] = (y * y + (y + h) * (y + h)) / 2.0;

	return bb;
}

// struct Visitor
// {
// 	int count;
// 	bool ContinueVisiting;

// 	Visitor() : count(0), ContinueVisiting(true){};

// 	void operator()(const RTree::Leaf *const leaf)
// 	{
// 		std::cout << "#" << count << ": visited " << leaf->leaf << " with bound " << leaf->bound.ToString() << std::endl;
// 		count++;
// 	}
// };

vector<BoundingBox> getFrontier(int n, RTree tree, vector<Node *> &nodes)
{
	// do a BFS traversal, add levels if no. of nodes in the level is less than n
	queue<Node> q;
	q.push(*tree.m_root);
	vector<BoundingBox> frontier;
	// vector<Node *> nodes;
	int no_of_visited_nodes = 1;
	while (no_of_visited_nodes < n && !q.empty())
	{
		Node node = q.front();
		q.pop();
		for (int i = 0; i < node.items.size(); i++)
		{
			if (node.hasLeaves)
				frontier.push_back(node.items[i]->bound);
			else
				q.push(*static_cast<Node *>(node.items[i]));
			no_of_visited_nodes++;
		}
	}
	while (!q.empty())
	{
		Node node = q.front();
		q.pop();
		frontier.push_back(node.bound);
		// nodes.emplace_back(static_cast<Node *>(node));
	}
	// print all nodes
	// for (auto i : nodes)
	// {
	// 	cout << i->bound.ToString() << endl;
	// }
	// cout << nodes.size() << endl;
	// return nodes;
	return frontier;
}

const double PI = 3.14159265358979323846;

double gaussian(double x[], double mean[], double var[], int n)
{
	double tmp = 1;
	double sum = 0;
	for (int i = 0; i < n; i++)
	{
		tmp *= var[i];
		sum += -0.5 * (x[i] - mean[i]) * (x[i] - mean[i]) / var[i];
	}
	return (1.0 / sqrt(pow(2 * PI, n) * tmp)) * exp(sum);
}
double pdq(vector<BoundingBox> frontier, BoundingBox query)
{
	int N = 0;
	double pdq = 0;
	for (auto i : frontier)
	{
		int n = i.noOfObjects();
		N += n;
		double mean[2] = {i.sum[0] / n, i.sum[1] / n};
		double var[2] = {i.sum_sq[0] / n - mean[0] * mean[0], i.sum_sq[1] / n - mean[1] * mean[1]};
		double x[2] = {query.sum[0], query.sum[1]};
		double g = gaussian(x, mean, var, 2);
		pdq += n * g;
	}
	return pdq / N;
}

void advanceFrontier(vector<Node *> &frontier, BoundingBox &query)
{
	// calculate the distance of each frontier from the query, store the one with the minimum distance
	double max_area = 1e9;
	Node *min_bound;
	// BoundingBox min_bound = frontier[0]->bound;
	for (auto node : frontier)
	{
		// find overlap area
		double area = query.overlap(node->bound);
		if (area < max_area)
		{
			max_area = area;
			min_bound = node;
		}
	}
	// remove the bound from the frontier
	frontier.erase(remove(frontier.begin(), frontier.end(), min_bound), frontier.end());
	// add all the children
	for (auto item : min_bound->items)
	{
		frontier.push_back(static_cast<Node *>(item));
	}
}

// take command line
int main()
{
	RTree tree;
	// Visitor x;
	freopen("test-cases/t2.txt", "r", stdin);
	int n;
	cin >> n;
	pair<int, int> X_train[n];
	// int Y_train[n];
	for (int i = 0; i < n; i++)
	{
		cin >> X_train[i].first >> X_train[i].second;
		tree.Insert(i, bounds(X_train[i].first, X_train[i].second, 0, 0));
	}
	vector<Node *> frontier;
	auto query = bounds(-2, -4, 0, 0);
	// auto query = bounds(2, 4, 0, 0);
	auto f = getFrontier(1, tree, frontier);
	// auto f = getFrontier(3, tree, frontier);
	magenta << pdq(f, query) << endl;

	// for (auto i : f)
	// {
	// 	cout << i.ToString() << endl;
	// }
	// BoundingBox bound = bounds(5, 10, 5, 5);
	// cout << "Searching in " << bound.ToString() << endl;
	// x = tree.Query(RTree::AcceptOverlapping(bound), Visitor());
	// vector<BoundingBox> f = tree.getFrontier(0);
	// tree.preOrderTraversal(tree.m_root);
	// red << tree.m_root->bound.ToString() << endl;
	// cout << "Visited " << x.count << " nodes" << endl;
}