// nomnomsolve.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <bitset>
#include <unordered_map>
#include <list>
#include <map>
#include <deque>
#include <queue>
#include <array>
#include <assert.h>
#include <type_traits>
#include <algorithm>
#include <functional>
//#include <boost/pool/object_pool.hpp>

#ifndef _OPENMP
#define _OPENMP 0
#endif

typedef unsigned __int8 uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;

#define TRACKPARENTS 1 // Whether to keep track of parents and print full game solutions (costs memory)
#define ESTIMBEST 27 // The estimated best move count
#define TESTNUMBER 16 // The test to run
#define MAXPLRS 4 // The most players possible in a single game
#define MAXSTACK 8 // Used for optimizing around maximum stacked items
#define CULLMMAP 0 // Use my custom map
#define CULLBMAP 1 // Use Seen Map
#define THREADBATCH 1000
#define MAXITEM 8 // Maximum amount of stacks of any single item
#define MAXITEMNUM 10 // Maximum item number
#define MAXCOST (MAXITEM * 8 * MAXPLRS * 2) // The highest possible cost returned from the heuristic

typedef uint8_t GRID[4][6];
#define GV(itm, num) (itm|((num)<<5))
#define GI(id,y,x) (id[y][x]&0x1F)
#define GN(id,y,x) (id[y][x]>>5)
#define GE(id,y,x) (id[y][x]==0)
#define GS(id,y,x,itm,num) (id[y][x]=GV(itm,num))
#define GR(id,y,x) (id[y][x]=0)

template<typename T>
class cppalloc {
public:
	T * malloc() { return (T*)new T; }
	void free(T * p) { delete p; }
};

/*
template<typename T>
class boostalloc {
	boost::object_pool<T> ba;
public:
	T * malloc() { return ba.malloc(); }
	void free(T * p) { ba.free(p); }
};
*/

cppalloc<class AStarOpen> nodePool;
cppalloc<GRID> gridPool;

#include "tests.h"

uint32_t lastStatus = 0;
volatile unsigned long powCount = 0;
volatile unsigned long poofCount = 0;
volatile unsigned long earlyCount = 0;
volatile unsigned long testCount = 0;
volatile unsigned long solutionCount = 0;
volatile unsigned long pendCount = 0;

class GridCmp {
public:
	bool operator()(const GRID* a, const GRID* b) {
		return memcmp(a, b, sizeof(GRID)) < 0;
	}
};

class AStarOpen {
public:
	GRID grd;

	volatile unsigned long refcount;
	unsigned char mindist;
	unsigned char moves;

#if TRACKPARENTS
	AStarOpen *parent;
#endif // TRACKPARENTS

	void addref() {
#pragma omp atomic
		refcount++;
	}
	void decref() {
#pragma omp atomic
		refcount--;
		if (refcount == 0) {
#if TRACKPARENTS
			if (parent) {
				parent->decref();
				parent = nullptr;
			}
#endif
			nodePool.free(this);
#pragma omp atomic
			poofCount++;
		}
	}

};

unsigned long hashing_func(GRID* key) {
	const unsigned long *bytes = (unsigned long*)key;
	unsigned long hash = 0;
	for (int i = 0; i < sizeof(GRID) / sizeof(unsigned long); ++i) {
		hash = hash * 37 + bytes[i];
	}
	return hash;
}
bool key_equal_fn(GRID* a, GRID* b) {
	return memcmp(a, b, sizeof(GRID)) == 0;
}

#if 0
typedef std::unordered_map <
	GRID*, uint8_t,
	std::function<unsigned long(GRID*)>,
	std::function<bool(GRID*, GRID*)>
> GRIDBMAP;
GRIDBMAP seenBMap(0xFFFF, &hashing_func, &key_equal_fn);
#else
typedef std::map<GRID*, uint8_t, GridCmp> GRIDBMAP;
GRIDBMAP seenBMap;
#endif
typedef std::deque<AStarOpen*> OPENLISTGRP;
typedef std::vector<AStarOpen*> OPENLIST;
OPENLISTGRP openListX[MAXCOST];
volatile unsigned long lowOpenCost = MAXCOST;
volatile unsigned long fastestFound = ESTIMBEST;

struct ThreadData {
	OPENLIST curOpenList;
	OPENLIST newOpenList;
	GRIDBMAP seenList;
} *tdata;
#pragma omp threadprivate(tdata)

const int pxs[5][4][2][2] = {
	{
		{{0, 0}, { 0, 0 }},
		{ { 0, 0 }, { 0, 0 } },
		{ { 0, 0 }, { 0, 0 } },
		{ { 0, 0 }, { 0, 0 } }
	}, {
			{ { 0, 1 }, { 0, 2 } },
			{ { 0, 0 }, { 0, 0 } },
			{ { 0, 0 }, { 0, 0 } },
			{ { 0, 0 }, { 0, 0 } }
	}, {
			{ { 0, 1 }, { 0, 2 } },
			{ { 5, 1 }, { 5, 2 } },
			{ { 0, 0 }, { 0, 0 } },
			{ { 0, 0 }, { 0, 0 } }
	}, {
			{ { 0, 0 }, { 0, 0 } },
			{ { 0, 0 }, { 0, 0 } },
			{ { 0, 0 }, { 0, 0 } },
			{ { 0, 0 }, { 0, 0 } }
	}, {
			{ { 0, 0 }, { 0, 1 } },
			{ { 0, 2 }, { 0, 3 } },
			{ { 5, 0 }, { 5, 1 } },
			{ { 5, 2 }, { 5, 3 } }
	}
};


struct FoundItem {
	unsigned char x;
	unsigned char y;
	char num;
};
struct FoundItemList {
	FoundItem items[MAXITEM];
	int count;
};
typedef FoundItemList FoundItemPre[MAXITEMNUM];

uint32_t startTime;
int xtime() {
	return hrTime() - startTime;
}

void copyGrid(GRID& dest, const GRID& src) {
	memcpy(&dest, &src, sizeof(GRID));
}

bool canEqualExact(char* items, uint32_t ignore, char cur) {
	const char *item = items;
	uint32_t key = 1;
	for (; *item; ++item, key <<= 1) {
		if (ignore&key) continue;
		if (cur - *item == 0) {
			return true;
		}
		if (canEqualExact(items, ignore | key, cur - *item)) {
			return true;
		}
	}
	return false;
}

bool isPossible(const GRID& in) {
	char items[MAXSTACK + 1];
	for (int j = 0; j < TTLPLAYRS; ++j) {
		for (int n = 0; n < 2; ++n) {
			memset(items, 0, sizeof(items));

			char *item = items;
			for (int ix = 0; ix < 6; ++ix) {
				for (int iy = 0; iy < 4; ++iy) {
					if (GI(in, iy, ix) == i[j][n][0]) {
						*item++ = GN(in, iy, ix);
					}
				}
			}

			if (!canEqualExact(items, 0, i[j][n][1])) {
				return false;
			}
		}
	}
	return true;
}

bool checkWin(GRID& in) {
	for (int iq = 0; iq < TTLPLAYRS; ++iq) {
		const int px1 = pxs[TTLPLAYRS][iq][0][0];
		const int py1 = pxs[TTLPLAYRS][iq][0][1];
		const int px2 = pxs[TTLPLAYRS][iq][1][0];
		const int py2 = pxs[TTLPLAYRS][iq][1][1];
		const uint8_t i1i = i[iq][0][0];
		const uint8_t i1n = i[iq][0][1];
		const uint8_t i2i = i[iq][1][0];
		const uint8_t i2n = i[iq][1][1];

		if (GI(in, py1, px1) == i1i && GI(in, py1, px1) == i2i &&
			GN(in, py1, px1) == i1n && GN(in, py1, px1) == i2n) {
			// Works in given order
		} else if (GI(in, py2, px2) == i2i && GI(in, py2, px2) == i1i &&
			GN(in, py2, px2) == i2n && GN(in, py2, px2) == i1n) {
			// Works in flipped order
		} else {
			return false;
		}
	}
	return true;
}

struct BrettTreeNode {};
struct BrettTreeLeaf : public BrettTreeNode {
	BrettTreeLeaf() {
		memset(value, 0, sizeof(value));
	}

	unsigned char value[0xFF];
};
struct BrettTreeBranch : public BrettTreeNode {
	BrettTreeBranch() {
		memset(nextNode, 0, sizeof(nextNode));
	}

	BrettTreeBranch* nextNode[0xFF];
};

template<typename T, int Size>
struct FastAllocator {
	struct GroupEntry {
		GroupEntry() : nextGroup(nullptr) {}

		T items[Size];
		GroupEntry *nextGroup;
	};

	FastAllocator() {
		currentConsumed = 0;
		currentGroup = new GroupEntry();
	}

	T * alloc() {
		if (currentConsumed >= Size) {
			auto newGroup = new GroupEntry;
			newGroup->nextGroup = currentGroup;
			currentGroup = newGroup;
			currentConsumed = 0;
		}
		return &currentGroup->items[currentConsumed++];
	}

	int currentConsumed;
	GroupEntry *currentGroup;
};

struct BrettTree {
	bool doCheck(const uint8_t in[24], unsigned char value) {
		BrettTreeBranch* thisNode = &baseNode;
		
		for (int i = 0; i < 22; ++i) {
			BrettTreeBranch*& nextNode = thisNode->nextNode[in[i]];
			if (!nextNode) {
				nextNode = branchAlloc.alloc();
			}
			thisNode = nextNode;
		}

		auto& tailBranch = thisNode->nextNode[in[22]];
		if (!tailBranch) {
			tailBranch = (BrettTreeBranch*)leafAlloc.alloc();
		}

		auto leafNode = (BrettTreeLeaf*)tailBranch;
		auto& leafValue = leafNode->value[in[23]];
		if (leafValue != 0 && leafValue <= value) {
			return true;
		}
		leafValue = value;
		return false;
	}
	
	FastAllocator<BrettTreeBranch, 0xFFFF> branchAlloc;
	FastAllocator<BrettTreeLeaf, 0xFFF> leafAlloc;
	BrettTreeBranch baseNode;
};

void moveItem(GRID& in, int fx, int fy, int tx, int ty) {
	GS(in,ty,tx, GI(in,fy,fx), GN(in,fy,fx)+GN(in,ty,tx));
	GR(in,fy,fx);
}

void printBoard(GRID& in, AStarOpen *parent);

struct CacheTableValue {
	unsigned char t2[6][4][MAXSTACK] [6][4][MAXSTACK];
	unsigned char t3[6][4][MAXSTACK] [6][4][MAXSTACK] [6][4][MAXSTACK];
};

struct CacheTables {
	CacheTableValue *t[MAXSTACK][4][2];
} caches;

template<int ListLeft, int ListTotal, typename T>
class FindClosestItemClass {
	static_assert(ListLeft > 0, "ListLeft must be greater than 0");

public:
	static int find(const FoundItemList& list, T ignore, unsigned char dist, unsigned char sx, unsigned char sy, unsigned char num) {
		int bestDist = -1;
		unsigned int key = 1;
		for (int i = 0; i < ListTotal; ++i, key <<= 1) {
			if (ignore & key) continue;
			const FoundItem& itm = list.items[i];

			const int ndist = abs(sx - itm.x) + abs(sy - itm.y);
			if (num < itm.num) {
				continue;
			} else if (num == itm.num) {
				if (ndist < bestDist || bestDist == -1) {
					bestDist = dist + ndist;
				}
				continue;
			}

			const int xdist = FindClosestItemClass<ListLeft - 1, ListTotal, T>::find(list, ignore | key, dist + ndist, itm.x, itm.y, num - itm.num);
			if (xdist >= 0 && (xdist < bestDist || bestDist == -1)) {
				bestDist = xdist;
			}
		}
		return bestDist;
	}
};
template<int ListTotal, typename T>
class FindClosestItemClass < 0, ListTotal, T > {
public:
	__forceinline
	static int find(const FoundItemList& list, T ignore, unsigned char dist, unsigned char sx, unsigned char sy, unsigned char num) {
		return -1;
	}
};

int findClosestItemStk(const FoundItemList& list, unsigned char sx, unsigned char sy, unsigned char num) {
	// Max Recursion MAXITEM

	struct {
		int sx;
		int sy;
		unsigned int key;
		int numLeft;
		int dist;
		unsigned int ignore;
		int i;

		__forceinline void setup(unsigned int _ignore, int _dist, int _sx, int _sy, int _numLeft) {
			sx = _sx;
			sy = _sy;
			key = 1;
			numLeft = _numLeft;
			dist = _dist;
			ignore = _ignore;
			i = 0;
		}
	} lclStack[MAXITEM];
	auto* stackBasePtr = lclStack;

	stackBasePtr->setup(0, 0, sx, sy, num);

	int bestDist = 0xFFFFFF;
	int ndist;
	while (1) {
		auto* stackPtr = stackBasePtr;

		if (stackPtr->i >= list.count) {
			if (stackPtr == lclStack) {
				break;
			} else {
				stackBasePtr--;
				continue;
			}
		}

		if (~stackPtr->ignore & stackPtr->key) {
			const FoundItem& itm = list.items[stackPtr->i];

			ndist = stackPtr->dist + abs(stackPtr->sx - itm.x) + abs(stackPtr->sy - itm.y);
			if (stackPtr->numLeft == itm.num) {
				if (ndist < bestDist) {
					bestDist = ndist;
				}
			} else if (num > itm.num) {
				(stackBasePtr + 1)->setup(stackPtr->ignore | stackPtr->key, ndist, itm.x, itm.y, stackPtr->numLeft - itm.num);
				stackBasePtr++;
			}
		}

		stackPtr->i++;
		stackPtr->key <<= 1;
	}

	return bestDist;
}

template<int ListLeft>
int findClosestItem(const FoundItemList& list, unsigned char sx, unsigned char sy, unsigned char num) {
	return findClosestItemStk(list, sx, sy, num);
	//return FindClosestItemClass<ListLeft, ListLeft, uint32_t>::find(list, 0, 0, sx, sy, num);
}

template<int ListLeft>
unsigned char findClosestItemX(const FoundItemList& list, unsigned char sx, unsigned char sy, unsigned char num) {
	int res = findClosestItem<ListLeft>(list, sx, sy, num);
	if (res < 0 || res > 0xFF) {
		return 0xFF;
	} else {
		return res;
	}
}

void buildTable() {
	FILE *rfh = NULL;
	fopen_s(&rfh, "table.cache", "rb");
	if (rfh != NULL) {
		for (int tn = 0; tn < MAXSTACK; ++tn) {
			for (int sx = 0; sx < 2; ++sx) {
				for (int sy = 0; sy < 4; ++sy) {
					CacheTableValue* v = new CacheTableValue;
					caches.t[tn][sy][sx] = v;
					fwrite(v, sizeof(CacheTableValue), 1, rfh);
				}
			}
		}
		fclose(rfh);
		return;
	}

	int xlkp[] = { 0, 5 };
	FoundItemList itms;
	memset(&itms, 0, sizeof(itms));

	for (int tn = 0; tn < MAXSTACK; ++tn) {
		for (int _sx = 0; _sx < 2; ++_sx) {
			int sx = xlkp[_sx];
			for (int sy = 0; sy < 4; ++sy) {
				CacheTableValue* v = new CacheTableValue();
				caches.t[tn][sy][_sx] = v;

				for (int ix1 = 0; ix1 < 6; ++ix1) {
					itms.items[0].x = ix1;
					for (int iy1 = 0; iy1 < 4; ++iy1) {
						itms.items[0].y = iy1;
						for (int in1 = 0; in1 < MAXSTACK; ++in1) {
							itms.items[0].num = in1;

							for (int ix2 = 0; ix2 < 6; ++ix2) {
								itms.items[1].x = ix2;
								for (int iy2 = 0; iy2 < 4; ++iy2) {
									itms.items[1].y = iy2;
									for (int in2 = 0; in2 < MAXSTACK; ++in2) {
										itms.items[1].num = in2;

										itms.count = 2;
										v->t2[ix1][iy1][in1][ix2][iy2][in2] = 
											findClosestItemX<2>(itms, sx, sy, tn);

										for (int ix3 = 0; ix3 < 6; ++ix3) {
											itms.items[2].x = ix3;
											for (int iy3 = 0; iy3 < 4; ++iy3) {
												itms.items[2].y = iy3;
												for (int in3 = 0; in3 < MAXSTACK; ++in3) {
													itms.items[2].num = in3;

													itms.count = 3;
													v->t3[ix1][iy1][in1][ix2][iy2][in2][ix3][iy3][in3] =
														findClosestItemX<3>(itms, sx, sy, tn);

												}
											}
										}
									}
								}
							}
						}
					}
				}

			}
		}
	}

	FILE *wfh = NULL;
	fopen_s(&wfh, "table.cache", "wb");
	if (wfh != NULL) {
		for (int tn = 0; tn < MAXSTACK; ++tn) {
			for (int sx = 0; sx < 2; ++sx) {
				for (int sy = 0; sy < 4; ++sy) {
					fwrite(caches.t[tn][sy][sx], sizeof(CacheTableValue), 1, wfh);
				}
			}
		}
		fclose(wfh);
	}
}

int findClosestItemFast(const FoundItemList& list, int sx, int sy, int num) {
#if 1
	const int xlkp[] = {0, -1, -1, -1, -1, 5};
	switch (list.count) {
	case 0:
		return 0xFFFFFF;
	case 1:
		if (list.items->num == num) {
			return abs(sx - list.items->x) + abs(sy - list.items->y);
		} else {
			return 0xFFFFFF;
		}
	case 2: {
		unsigned char res = caches.t[num][sy][xlkp[sx]]->t2
			[list.items[0].x][list.items[0].y][list.items[0].num]
			[list.items[1].x][list.items[1].y][list.items[1].num];
		return res != 0xFF ? res : 0xFFFFFF;
	}
	case 3: {
		unsigned char res = caches.t[num][sy][xlkp[sx]]->t3
			[list.items[0].x][list.items[0].y][list.items[0].num]
			[list.items[1].x][list.items[1].y][list.items[1].num]
			[list.items[2].x][list.items[2].y][list.items[2].num];
		return res != 0xFF ? res : 0xFFFFFF;
	}
	default:
		return findClosestItemStk(list, sx, sy, num);
	}
#else
	return findClosestItem(list, numList, bDist, sx, sy, num);
#endif
}

int findClosestPlrItem(const FoundItemPre& items, int iq) {
	const int& sx1 = pxs[TTLPLAYRS][iq][0][0];
	const int& sy1 = pxs[TTLPLAYRS][iq][0][1];
	const int& itm1 = i[iq][0][0];
	const int& num1 = i[iq][0][1];
	const int& sx2 = pxs[TTLPLAYRS][iq][1][0];
	const int& sy2 = pxs[TTLPLAYRS][iq][1][1];
	const int& itm2 = i[iq][1][0];
	const int& num2 = i[iq][1][1];

	const FoundItemList& items1 = items[itm1];
	const FoundItemList& items2 = items[itm2];

	const int cost11 = findClosestItemFast(items1, sx1, sy1, num1);
	const int cost12 = findClosestItemFast(items1, sx2, sy2, num1);

	const int cost21 = findClosestItemFast(items2, sx2, sy2, num2);
	const int cost22 = findClosestItemFast(items2, sx1, sy1, num2);

	int ttlCost = std::min(cost11, cost12) + std::min(cost21, cost22);
	if (ttlCost >= 0xFFFFF) {
		return -1;
	} else {
		return ttlCost;
	}
}

int calcDist(const GRID& in) {
	FoundItemPre precalcList;
	memset(precalcList, 0, sizeof(precalcList));

	for (int iy = 0; iy < 4; ++iy) {
		for (int ix = 0; ix < 6; ++ix) {
			const int ii = GI(in, iy, ix);
			if (ii == 0) continue;
			const int ij = precalcList[ii].count++;
			FoundItem& item = precalcList[ii].items[ij];
			item.x = ix;
			item.y = iy;
			item.num = GN(in, iy, ix);
		}
	}

	int ttlCost = 0;
	for (int iq = 0; iq < TTLPLAYRS; ++iq) {
		int cost = findClosestPlrItem(precalcList, iq);
		if (cost == -1) {
			return -1;
		}
		ttlCost += cost;
	}
	assert(ttlCost < MAXCOST);
	return ttlCost;
}

void printBoard(GRID& in, AStarOpen *parent) {
	for (int iy = 0; iy < 4; ++iy) {
		printf("| ");
		for (int ix = 0; ix < 6; ++ix) {
			printf("%02dx%02d ", GN(in,iy,ix), GI(in,iy,ix));
		}

		if (parent) {
			printf ("| | ");
			GRID& pg = parent->grd;
			for (int ipx = 0; ipx < 6; ipx++) {
				if (GI(in,iy,ipx) != GI(pg,iy,ipx) || GN(in,iy,ipx) != GN(pg,iy,ipx)) {
					if (GE(in,iy,ipx)) {
						printf("----- ");
					} else {
						printf("%02dx%02d ", GN(in,iy,ipx), GI(in,iy,ipx));
					}
				} else {
					printf("      ");
				}
			}
		}

		printf("|\n");
	}
}

#if TRACKPARENTS
void recursePrintBoard(AStarOpen *itm) {
	if (itm->parent) {
		recursePrintBoard(itm->parent);
	}

	printf("\\_____________________________________/\n");
	printBoard(itm->grd, itm->parent);
	printf("\\_____________________________________/ \\_____________________________________/\n");
}
#else
void recursePrintBoard(AStarOpen *itm) {
	printf("\\_____________________________________/\n");
	printBoard(itm->grd, nullptr);
	printf("\\_____________________________________/\n");
}
#endif

#if CULLBMAP && _OPENMP
__forceinline bool wasGridSeen(const GRID& grd, int moves) {
	GRIDBMAP::iterator seenItm = seenBMap.find((GRID*)&grd);
	if (seenItm != seenBMap.end()) {
		if (seenItm->second <= moves) {
			return true;
		}
		seenItm->second = moves;
		return false;
	}
	seenItm = tdata->seenList.find((GRID*)&grd);
	if (seenItm != tdata->seenList.end()) {
		if (seenItm->second <= moves) {
			return true;
		}
		seenItm->second = moves;
		return false;
	}

	GRID *ngrid = gridPool.malloc();
	memcpy(ngrid, &grd, sizeof(GRID));
	tdata->seenList.insert(std::pair<GRID*,unsigned int>(ngrid, moves));
	return false;
}
#elif CULLMMAP
BrettTree myMap;
__forceinline bool wasGridSeen(const GRID& grd, int moves) {
	return myMap.doCheck(&grd[0][0], moves);
}
#elif CULLBMAP
bool wasGridSeen(const GRID& grd, int moves) {
	GRIDBMAP::iterator seenItm = seenBMap.find((GRID*)&grd);
	if (seenItm != seenBMap.end()) {
		if (seenItm->second <= moves) {
			return true;
		}
		seenItm->second = moves;
	} else {
		GRID *ngrid = gridPool.malloc();
		memcpy(ngrid, &grd, sizeof(GRID));
		seenBMap.insert(std::pair<GRID*,unsigned int>(ngrid, moves));
	}
	return false;
}
#else
__forceinline bool wasGridSeen(const GRID& grd, int moves) {
	return false;
}
#endif

void addToOpen(const GRID& grd, unsigned int moves, AStarOpen *parent) {
	/*
	if (!isPossible(grd)) {
		earlyCount++;
		return;
	}
	*/

	int minCost = calcDist(grd);
	if (minCost == -1) {
		earlyCount++;
		return;
	}

	assert(minCost >= 0 && minCost < MAXCOST);
	assert(moves + minCost < MAXCOST);

	if (moves + minCost >= fastestFound) {
		earlyCount++;
		return;
	}

	if (wasGridSeen(grd, moves)) {
		earlyCount++;
		return;
	}

	AStarOpen *itm = nodePool.malloc();

#pragma omp atomic
	powCount++;

	copyGrid(itm->grd, grd);
	itm->mindist = minCost;
	itm->moves = moves;
	itm->refcount = 1;

#if TRACKPARENTS
	if (parent) {
		parent->addref();
		itm->parent = parent;
	} else {
		itm->parent = nullptr;
	}
#endif

#if _OPENMP
	tdata->newOpenList.push_back(itm);
#else
	openListX[minCost].push_back(itm);
	if (minCost < (int)lowOpenCost) {
		lowOpenCost = minCost;
	}
#endif
}

int addMoves(AStarOpen *owner) {
	GRID & in = owner->grd;
	int moves = owner->moves;

	testCount++;
#if !_OPENMP

	static unsigned int lasttick = 0;
	if (testCount % 10000 == 0) {
		unsigned int thistick = hrTime();

		int kmps = -1;
		if (thistick - lasttick > 0) {
			kmps = 10000 / (thistick - lasttick);
			lasttick = thistick;
		}
		printf("[% 7d] tested %dk (%d best) [% 5dkmps %dk,%dk,%dk,%dk]\n", xtime(), testCount/1000, solutionCount?fastestFound:-1, kmps, powCount/1000, poofCount/1000, (powCount-poofCount)/1000, earlyCount/1000);
	}
#endif

	if (owner->moves + owner->mindist >= (int)fastestFound) {
		return 3;
	}

	if (checkWin(in)) {
		recursePrintBoard(owner);
		printf("[% 7d] Found winner using %d moves!\n", xtime(), moves);

		fastestFound = moves;
		solutionCount++;

#if 0
		for (int z = 255; z >= fastestFound; --z) {
			for(OPENLISTGRP::iterator j = openListX[z].begin(); j != openListX[z].end(); ++j) {
				(*j)->decref();
			}
			openListX[z].clear();
		}
#endif

		return 1000;
	}

	GRID tgrd;

	for (int iy = 0; iy < 4; ++iy) {
		for (int ix = 0; ix < 6; ++ix) {
			if (GE(in,iy,ix)) {
				// Cell is empty, can't do anything
				continue;
			}

			// Up
			if (iy >= 1) {
				if (GE(in,iy-1,ix) || GI(in,iy-1,ix) == GI(in,iy,ix)) {
					// Move
					copyGrid(tgrd, in);
					moveItem(tgrd, ix, iy, ix, iy-1);
					addToOpen(tgrd, owner->moves+1, owner);
				}
			}

			// Down
			if (iy <= 2) {
				if (GE(in,iy+1,ix) || GI(in,iy+1,ix) == GI(in,iy,ix)) {
					// Move
					copyGrid(tgrd, in);
					moveItem(tgrd, ix, iy, ix, iy+1);
					addToOpen(tgrd, owner->moves+1, owner);
				}
			}

			// Left
			if (ix >= 1) {
				if (GE(in,iy,ix-1) || GI(in,iy,ix-1) == GI(in,iy,ix)) {
					// Move
					copyGrid(tgrd, in);
					moveItem(tgrd, ix, iy, ix-1, iy);
					addToOpen(tgrd, owner->moves+1, owner);
				}
			}

			// Right
			if (ix <= 4) {
				if (GE(in,iy,ix+1) || GI(in,iy,ix+1) == GI(in,iy,ix)) {
					// Move
					copyGrid(tgrd, in);
					moveItem(tgrd, ix, iy, ix+1, iy);
					addToOpen(tgrd, owner->moves+1, owner);
				}
			}
		}
	}

	return 0;
}

void doSearch() {
#if _OPENMP
#pragma omp parallel
	while(true) {
		uint32_t batchTotal = 0;

#pragma omp critical
		{
			// Merge everything back to global state, we do this to pull
			//   anything added by the master thread in the beginning into
			//   the global state
			for (OPENLIST::iterator i = tdata->newOpenList.begin(); i != tdata->newOpenList.end(); ++i) {
				openListX[(*i)->mindist].push_back(*i);
				if ((*i)->mindist < lowOpenCost) {
					lowOpenCost = (*i)->mindist;
				}
			}
			tdata->newOpenList.clear();

			seenBMap.insert(tdata->seenList.begin(), tdata->seenList.end());
			tdata->seenList.clear();

			// Create a list of things for this thread
			batchTotal = 0;
			while (lowOpenCost < 256 && batchTotal < THREADBATCH) {
				if (openListX[lowOpenCost].empty()) {
					lowOpenCost++;
					continue;
				}

				AStarOpen* itm = openListX[lowOpenCost].front();
				openListX[lowOpenCost].pop_front();

				tdata->curOpenList.push_back(itm);
				batchTotal++;
			}
		}

#pragma omp atomic
		pendCount += batchTotal;

		uint32_t curTime = hrTime();
		if (curTime - lastStatus > 1000) {
			lastStatus = curTime;
			printf("[% 7d] tested %dk (%d best) [%dk %dk %dk %dk %dk]\n", 
				xtime(),
				testCount/1000,
				solutionCount?fastestFound:-1,
				powCount/1000,
				poofCount/1000,
				(powCount-poofCount)/1000,
				earlyCount/1000,
				pendCount/1000);
		}

		if (batchTotal == 0) {
#if 0
			printf("T%d ran out of items\n", omp_get_thread_num());
			break;
#else
			if (pendCount == 0) {
				break;
			}
			sleepWait(10);
			continue;
#endif
		}

		//printf("[% 7d] T%d executing %d items\n", xtime(), omp_get_thread_num(), tdata->curOpenList.size());

		while (!tdata->curOpenList.empty()) {
			AStarOpen *itm = tdata->curOpenList.back();
			tdata->curOpenList.pop_back();

			addMoves(itm);
			itm->decref();
		}

		//printf("[% 7d] T%d finished executing\n", xtime(), omp_get_thread_num());

#pragma omp atomic
		pendCount -= batchTotal;
	}

#else
	while (lowOpenCost < 256) {
		if (openListX[lowOpenCost].empty()) {
			lowOpenCost++;
			continue;
		}

		AStarOpen* itm = openListX[lowOpenCost].front();
		openListX[lowOpenCost].pop_front();

		addMoves(itm);
		itm->decref();
	}
#endif
}

int _tmain(int argc, _TCHAR* argv[])
{
	omp_set_num_threads(omp_get_num_procs() / 2 - 1);

	startTime = hrTime();

	printf("ESTIMBEST: %d\n", ESTIMBEST);
	printf("MAXCOST: %d\n", MAXCOST);
	printf("TESTNUMBER: %d\n", TESTNUMBER);
	printf("MAXPLRS: %d\n", MAXPLRS);
	printf("MAXSTACK: %d\n", MAXSTACK);
	printf("CULLBMAP: %d\n", CULLBMAP);
	printf("_OPENMP: %d\n", _OPENMP);
	printf("TTLPLAYRS: %d\n", TTLPLAYRS);

	printf("sizeof(GRID): %d\n", sizeof(GRID));
	printf("sizeof(AStarOpen): %d\n", sizeof(AStarOpen));
	printf("sizeof(openListX): %d\n", sizeof(openListX));
	printf("\n");

	printf("canEqualExact(2): %d [0]\n", canEqualExact(imptest, 0, 2));
	printf("canEqualExact(5): %d [1]\n", canEqualExact(imptest, 0, 5));
	printf("canEqualExact(7): %d [1]\n", canEqualExact(imptest, 0, 7));
	printf("canEqualExact(11): %d [0]\n", canEqualExact(imptest, 0, 11));
	printf("checkWin(winTest): %d\n", checkWin(winTest));
	printf("checkWin(start): %d\n", checkWin(start));
	printf("\n");

	buildTable();
	printf("[% 7d] Cache Map Ready\n", xtime());

#if _OPENMP
#pragma omp parallel
	tdata = new ThreadData();
#endif

	addToOpen(start, 0, nullptr);
	doSearch();
	printf("[% 7d] Completed! [%d %d %d %d]\n", xtime(), powCount, poofCount, powCount-poofCount, earlyCount);

#if _OPENMP
#pragma omp parallel
	delete tdata;
#endif

	system("PAUSE");
	return 0;
}

