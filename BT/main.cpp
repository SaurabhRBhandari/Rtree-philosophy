#include <bits/stdc++.h>
using namespace std;
#include "utils.h"
#include "RStarTree.h"
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

struct Visitor
{
	int count;
	bool ContinueVisiting;

	Visitor() : count(0), ContinueVisiting(true){};

	void operator()(const RTree::Leaf *const leaf)
	{
		std::cout << "#" << count << ": visited " << leaf->leaf << " with bound " << leaf->bound.ToString() << std::endl;
		count++;
	}
};

int main()
{
	RTree tree;
	Visitor x;
	freopen("test-cases/t2.txt", "r", stdin);
	int n;
	cin >> n;
	pair<int, int> X_train[n];
	int Y_train[n];
	for (int i = 0; i < n; i++)
	{
		cin >> X_train[i].first >> X_train[i].second >> Y_train[i];
		tree.Insert(i, bounds(X_train[i].first, X_train[i].second, 0, 0));
	}
	// BoundingBox bound = bounds(5, 10, 5, 5);
	// cout << "Searching in " << bound.ToString() << endl;
	// x = tree.Query(RTree::AcceptOverlapping(bound), Visitor());
	vector<BoundingBox> f = tree.getFrontier(0);
	// auto f = tree.getFrontier(9);
	// cout << f.size() << endl;
	// tree.preOrderTraversal(tree.m_root);
	// red << tree.m_root->bound.ToString() << endl;
	// cout << "Visited " << x.count << " nodes" << endl;
}