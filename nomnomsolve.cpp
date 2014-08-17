// nomnomsolve.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <unordered_map>
#include <list>
#include <map>
#include <deque>
#include <queue>
#include <array>
#include <assert.h>
#include <boost/pool/object_pool.hpp>

unsigned __int32 crc32_8bytes(const void* data, size_t length, unsigned __int32 previousCrc32 = 0);
unsigned __int32 crc32_4bytes(const void* data, size_t length, unsigned __int32 previousCrc32 = 0);

#define FASTGRID 1
#define TRACKPARENTS 1
#define CULLMAP 0
#define CULLLIST 0
#define CULLBMAP 1
#define RETRACEUP 0

#define MAXCOST 256

#if FASTGRID
typedef unsigned char GRID[4][6];
#define GV(itm, num) ((itm&0x1F)|((num&0x7)<<5))
#define GI(id,y,x) (id[y][x]&0x1F)
#define GN(id,y,x) ((id[y][x]>>5)&0x7)
#define GE(id,y,x) (id[y][x]==0)
#else
typedef unsigned char GRID[4][6][2];
#define GV(itm, num) {itm,num}
#define GI(id) id[0]
#define GN(id) id[1]
#define GE(id) (id[1]==0)
#endif

unsigned int powCount = 0;
unsigned int poofCount = 0;
unsigned int earlyCount = 0;

template<typename T>
class cppalloc {
public:
    T * malloc() {
        return new T;
    }
    void free(T * p) {
        delete p;
    }
};

template<typename T>
class boostalloc {
    boost::object_pool<T> ba;
public:
    T * malloc() {
        return ba.malloc();
    }
    void free(T * p) {
        ba.free(p);
    }
};


cppalloc<class AStarOpen> nodePool;
//boostalloc<class AStarOpen> nodePool;
boost::object_pool<GRID> gridPool;

class AStarOpen {
public:
#if TRACKPARENTS
    AStarOpen() : mindist(0), moves(0), refcount(1), parent(nullptr) {}
#else
    AStarOpen() : mindist(0), moves(0), refcount(1) {}
#endif

    ~AStarOpen() {
    }

    GRID grd;

    unsigned int refcount;
    //unsigned char dist;
    unsigned char mindist;
    unsigned char moves;

#if TRACKPARENTS
    AStarOpen *parent;
#endif // TRACKPARENTS

    void addref() { refcount++; }
    void decref() {
        refcount--;
        if (refcount == 0) {
#if TRACKPARENTS
            if (parent) {
                parent->decref();
                parent = nullptr;
            }
#endif
            nodePool.free(this);
            poofCount++;
        }
    }
    
};
typedef std::list<AStarOpen*> AStarOpenList;


int fastestFound = 400;

#define TESTNUMBER 10
#include "tests.h"


const int pxs[5][4][2][2] = {
    {
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}}
    }, {
        {{0,1}, {0,2}},
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}}
    }, {
        {{0,1}, {0,2}},
        {{5,1}, {5,2}},
        {{0,0}, {0,0}},
        {{0,0}, {0,0}}
    }, {
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}}
    }, {
        {{0,0}, {0,1}},
        {{0,2}, {0,3}},
        {{5,0}, {5,1}},
        {{5,2}, {5,3}}
    }
};


int startTime;
int xtime() {
    return GetTickCount() - startTime;
}

void copyGrid(GRID& in, GRID& out) {
    memcpy(&out, &in, sizeof(GRID));
}

bool checkWin(GRID& in) {
    for (int iq = 0; iq < ttlPlyrs; ++iq) {
        int px1 = pxs[ttlPlyrs][iq][0][0];
        int py1 = pxs[ttlPlyrs][iq][0][1];
        int px2 = pxs[ttlPlyrs][iq][1][0];
        int py2 = pxs[ttlPlyrs][iq][1][1];

        if ((GI(in,py1,px1) == i[iq][0][0] && GN(in,py1,px1) == i[iq][0][1] && GI(in,py2,px2) == i[iq][1][0] && GN(in,py2,px2) == i[iq][1][1]) ||
            (GI(in,py2,px2) == i[iq][0][0] && GN(in,py2,px2) == i[iq][0][1] && GI(in,py1,px1) == i[iq][1][0] && GN(in,py1,px1) == i[iq][1][1])) {
            // Matches so far...
        } else {
            return false;
        }
    }
    return true;
}

typedef std::unordered_map<unsigned int, unsigned char> HASHMAP;

int testCount = 0;
HASHMAP hasSeen;
//AStarOpenList openList;
int solutionCount = 0;

struct SeenGrid {
    SeenGrid(GRID& _grd, unsigned int _moves)
        : moves(_moves) {
            memcpy(&grd, &_grd, sizeof(GRID));
    }

    GRID grd;
    unsigned int moves;
};

class GridCmp {
public:
    typedef unsigned __int64 uint64_t;
    bool operator()(const GRID* a, const GRID* b) {
       return memcmp(a, b, sizeof(GRID)) < 0;
    }
};

std::map<unsigned __int64, std::vector<SeenGrid>> seenGrids;
typedef std::map<GRID*, unsigned char, GridCmp> GRIDBMAP;
GRIDBMAP seenBMap;
typedef std::deque<AStarOpen*> OPENLISTGRP;
OPENLISTGRP openListX[MAXCOST];
int lowOpenCost = 256;

unsigned __int64 gridHashKey(GRID& in) {
    unsigned char *bytes = (unsigned char*)&in;
    unsigned __int64 hash = bytes[0] & 3;
    for (int i = 1; i < 21; ++i) {
        hash = (hash << 3) | (bytes[i] & 7);
    }
    return hash;
}

unsigned int hashGrid(GRID& in) {
    return crc32_4bytes(&in, sizeof(GRID));
}

void moveItem(GRID& in, int fx, int fy, int tx, int ty) {
#if FASTGRID
    in[ty][tx] = GV(GI(in,fy,fx), GN(in,fy,fx)+GN(in,ty,tx));
    in[fy][fx] = 0;
#else
    in[ty][tx][0] = in[fy][fx][0];
    in[ty][tx][1] += in[fy][fx][1];
    in[fy][fx][0] = 0;
    in[fy][fx][1] = 0;
#endif
}

#define XYB(x, y) (1<<(y*6+x))

#if 0
void findClosestItem(GRID& in, int& minCost, int& rCost, int sx, int sy, int itm, int& foundNum, unsigned int &ignore) {
    int bestDist = 10000;
    int foundExtra = 0;
    int bestX = 0;
    int bestY = 0;
    for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 6; ++ix) {
            if (ignore & XYB(ix,iy)) {
                continue;
            }

            if (GI(in,iy,ix) == itm) {
                int fDist = abs(ix-sx) + abs(iy-sy);
                if (fDist < bestDist) {
                    bestDist = fDist;
                    foundExtra = GN(in,iy,ix) - 1;
                    bestX = ix;
                    bestY = iy;
                }
            }
        }
    }
    assert(bestDist < 4 + 6);

    foundNum = foundExtra;
    ignore |= XYB(bestX, bestY);

    

#if 1
    minCost += bestDist;
    rCost += bestDist;
#else
    minCost += bestDist;
    rCost += bestDist * 2;

    if (bestX < sx) {
        std::swap(bestX, sx);
    }
    for (int ix = sx; ix < bestX; ++ix) {
        bool foundHole = false;
        for (int iy = -1; iy <= 1; ++iy) {
            if (bestY+iy < 0 || bestY+iy >= 4) {
                continue;
            }
            if (GE(in, ix, bestY+iy)) {
                foundHole = true;
                break;
            }
        }
        if (!foundHole) {
            rCost++;
        }
    }
#endif
}

unsigned int gnoSearch;
void findClosestItem(GRID& in, int& minCost, int& rCost, int sx, int sy, int itm, int num) {
    gnoSearch = 0;
#if 1
    int totalDist = 0;
    int foundExtra = 0;
    assert(num <= 7);
    for (int i = 0; i < num; ++i) {
        findClosestItem(in, minCost, rCost, sx, sy, itm, foundExtra, gnoSearch);
        i += foundExtra;
    }
#else
    findClosestItem(in, minCost, rCost, sx, sy, itm, gnoSearch);
#endif
}

#else
#endif

void printBoard(GRID& in, AStarOpen *parent);

struct FoundItem {
    unsigned char x;
    unsigned char y;
    char num;
};

#define MAXPLRS 4
#define MAXSTACK 6

int findClosestItem(const FoundItem *list, int numList, unsigned int ignore, int dist, int bDist, int sx, int sy, int num);

unsigned char cacheTable2[MAXSTACK] [2][4] [6][4][MAXSTACK] [6][4][MAXSTACK];
unsigned char cacheTable3[MAXSTACK] [2][4] [6][4][MAXSTACK] [6][4][MAXSTACK] [6][4][MAXSTACK];
void buildTable() {
    int xlkp[] = { 0, 5 };
    FoundItem itms[3];

    for (int tn = 0; tn < MAXSTACK; ++tn) {
        for (int _sx = 0; _sx < 2; ++_sx) {
            int sx = xlkp[_sx];
            for (int sy = 0; sy < 4; ++sy) {

                for (int ix1 = 0; ix1 < 6; ++ix1) {
                    itms[0].x = ix1;
                    for (int iy1 = 0; iy1 < 4; ++iy1) {
                        itms[0].y = iy1;
                        for (int in1 = 0; in1 < MAXSTACK; ++in1) {
                            itms[0].num = in1;
                            for (int ix2 = 0; ix2 < 6; ++ix2) {
                                itms[1].x = ix2;
                                for (int iy2 = 0; iy2 < 4; ++iy2) {
                                    itms[1].y = iy2;
                                    for (int in2 = 0; in2 < MAXSTACK; ++in2) {
                                        itms[1].num = in2;

                                        cacheTable2[tn][_sx][sy][ix1][iy1][in1][ix2][iy2][in2] = 
                                            findClosestItem(itms, 2, 0, 0, -1, sx, sy, tn);

                                        for (int ix3 = 0; ix3 < 6; ++ix3) {
                                            itms[2].x = ix2;
                                            for (int iy3 = 0; iy3 < 4; ++iy3) {
                                                itms[2].y = iy2;
                                                for (int in3 = 0; in3 < MAXSTACK; ++in3) {
                                                    itms[2].num = in2;

                                                    cacheTable3[tn][_sx][sy][ix1][iy1][in1][ix2][iy2][in2][ix3][iy3][in3] = 
                                                        findClosestItem(itms, 3, 0, 0, -1, sx, sy, tn);
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
}

int findClosestItem(const FoundItem *list, int numList, unsigned int ignore, int dist, int bDist, int sx, int sy, int num) {
    int bestDist = -1;
    unsigned int key = 1;
    const FoundItem *mlist = list;
    for (; mlist->num > 0; mlist++, key <<= 1) {
        if (ignore & key) continue;
        const FoundItem& itm = *mlist;

        int nnum = num - itm.num;
        if (nnum < 0) continue;

        int ndist = abs(sx-itm.x) + abs(sy-itm.y);
        if (nnum == 0) return dist + ndist;
        if (bDist > -1 && dist + ndist >= bDist) {
            continue;
        }

        int xdist = findClosestItem(list, numList, ignore | key, dist+ndist, bestDist, itm.x, itm.y, nnum);
        if (xdist >= 0 && (bestDist == -1 || xdist < bestDist)) {
            bestDist = xdist;
        }
    }
    return bestDist;
}

int findClosestItemFast(const FoundItem *list, int numList, unsigned int ignore, int dist, int bDist, int sx, int sy, int num) {
    if (numList == 1) {
        return abs(sx-list->x) + abs(sy-list->y);
    } else if (numList == 2) {
        return cacheTable2[num]
            [sx==0?0:1][sy]
            [list[0].x][list[0].y][list[0].num]
            [list[1].x][list[1].y][list[1].num];
    } else if (numList == 3) {
        return cacheTable3[num]
            [sx==0?0:1][sy]
            [list[0].x][list[0].y][list[0].num]
            [list[1].x][list[1].y][list[1].num]
            [list[2].x][list[2].y][list[2].num];
    } else {
        return findClosestItem(list, numList, ignore, dist, bDist, sx, sy, num);
    }
}

int findClosestPlrItem(const GRID& in, int iq) {
    int sx1 = pxs[ttlPlyrs][iq][0][0];
    int sy1 = pxs[ttlPlyrs][iq][0][1];
    int itm1 = i[iq][0][0];
    int num1 = i[iq][0][1];
    int sx2 = pxs[ttlPlyrs][iq][1][0];
    int sy2 = pxs[ttlPlyrs][iq][1][1];
    int itm2 = i[iq][1][0];
    int num2 = i[iq][1][1];

    int itemCnt1 = 0;
    FoundItem items1[MAXSTACK];
    int itemCnt2 = 0;
    FoundItem items2[MAXSTACK];
    
    memset(items1, 0, sizeof(items1));
    memset(items2, 0, sizeof(items2));

    for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 6; ++ix) {
            if (GI(in,iy,ix) == itm1) {
                items1[itemCnt1].x = ix;
                items1[itemCnt1].y = iy;
                items1[itemCnt1].num = GN(in,iy,ix);
                itemCnt1++;
            }
            if (GI(in,iy,ix) == itm2) {
                items2[itemCnt2].x = ix;
                items2[itemCnt2].y = iy;
                items2[itemCnt2].num = GN(in,iy,ix);
                itemCnt2++;
            }
        }
    }

    int cost11 = findClosestItem(items1, itemCnt1, 0, 0, -1, sx1, sy1, num1);
    int cost12 = findClosestItem(items1, itemCnt1, 0, 0, -1, sx2, sy2, num1);
    if (cost11 == -1 && cost12 == -1) return -1;

    int cost21 = findClosestItem(items2, itemCnt2, 0, 0, -1, sx2, sy2, num2);
    int cost22 = findClosestItem(items2, itemCnt2, 0, 0, -1, sx1, sy1, num2);
    if (cost21 == -1 && cost22 == -1) return -1;

    int ttlCost = 0;
    if (cost11 == -1 || cost12 < cost11) {
        ttlCost += cost12;
    } else {
        ttlCost += cost11;
    }
    if (cost21 == -1 || cost22 < cost21) {
        ttlCost += cost22;
    } else {
        ttlCost += cost21;
    }
    return ttlCost;
}

int calcDist(const GRID& in) {
    int ttlCost = 0;
    for (int iq = 0; iq < ttlPlyrs; ++iq) {
        int cost = findClosestPlrItem(in, iq);
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

void recursePrintBoard(AStarOpen *itm) {
#if TRACKPARENTS
    if (itm->parent) {
        recursePrintBoard(itm->parent);
    }

    printf("\\_____________________________________/\n");
    printBoard(itm->grd, itm->parent);
    printf("\\_____________________________________/ \\_____________________________________/\n");
#else
    printf("\\_____________________________________/\n");
    printBoard(itm->grd, nullptr);
    printf("\\_____________________________________/\n");
#endif
}

void addToOpen(GRID& grd, unsigned int moves, AStarOpen *parent) {
    int minCost = calcDist(grd);
    if (minCost == -1) {
        earlyCount++;
        return;
    }

    assert(minCost >= 0 && minCost < MAXCOST);

    if (moves+minCost < fastestFound) {
#if CULLMAP
        unsigned int hashTag = hashGrid(itm.grd);
        HASHMAP::iterator prevSeen = hasSeen.find(hashTag);
        if (prevSeen != hasSeen.end()) {
            if (itm.moves >= prevSeen->second) {
                itm->decref();
                earlyCount++;
                return;
            }
            prevSeen->second = itm.moves;
        } else {
            hasSeen.insert(std::pair<unsigned int, unsigned char>(hashTag, itm.moves));
        }
#endif
#if CULLLIST
        unsigned __int64 hashKey = gridHashKey(itm.grd);
        std::vector<SeenGrid>& hashGrids = seenGrids[hashKey];
        bool foundMatch = false;
        for (std::vector<SeenGrid>::iterator i = hashGrids.begin(); i != hashGrids.end(); ++i) {
            if (memcmp(&i->grd, &itm.grd, sizeof(GRID)) == 0) {
                if (i->moves < itm.moves) {
                    itm->decref();
                    earlyCount++;
                    return;
                }
                i->moves = itm.moves;
                foundMatch = true;
                break;
            }
        }
        if (!foundMatch) {
            hashGrids.push_back(SeenGrid(itm.grd, itm.moves));
        }
#endif
#if CULLBMAP
        GRIDBMAP::iterator seenItm = seenBMap.find(&grd);
        if (seenItm != seenBMap.end()) {
            //assert(memcmp(seenItm->first, &itm.grd, sizeof(GRID)) == 0);
            if (seenItm->second <= moves) {
                earlyCount++;
                return;
            }
            seenItm->second = moves;
        } else {
            GRID *ngrid = gridPool.malloc();
            memcpy(ngrid, &grd, sizeof(GRID));
            seenBMap.insert(std::pair<GRID*,unsigned int>(ngrid, moves));
        }
#endif
#if RETRACEUP & TRACKPARENTS
        // Skip direct parent, since we had to do something to get here...
        if (parent) {
            AStarOpen *up = parent->parent;
            while (up) {
                if (memcmp(&itm.grd, &up->grd, sizeof(GRID)) == 0) {
                    itm->decref();
                    earlyCount++;
                    return;
                }
                up = up->parent;
            }
        }
#endif

        AStarOpen *itm = nodePool.malloc();
        powCount++;
        copyGrid(grd, itm->grd);
        //itm->dist = rCost;
        itm->mindist = moves + minCost;
        assert(itm->mindist < MAXCOST);
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
        openListX[minCost].push_back(itm);
        if (minCost < lowOpenCost) {
            lowOpenCost = minCost;
        }
    } else {
        earlyCount++;
    }
}

unsigned int lasttick = 0;
int addMoves(AStarOpen *owner) {
    GRID & in = owner->grd;
    int moves = owner->moves;

    testCount++;
    if (testCount % 10000 == 0) {
        unsigned int thistick = GetTickCount();

        int kmps = -1;
        if (thistick - lasttick > 0) {
            kmps = 10000 / (thistick - lasttick);
            lasttick = thistick;
        }
        printf("[% 7d] tested %dk (%d best) [% 5dkmps %dk,%dk,%dk,%dk]\n", xtime(), testCount/1000, solutionCount?fastestFound:-1, kmps, powCount/1000, poofCount/1000, (powCount-poofCount)/1000, earlyCount/1000);
    }

    if (owner->mindist >= fastestFound) {
        return 3;
    }

    if (checkWin(in)) {
        //recursePrintBoard(owner);
        printf("[% 7d] Found winner using %d moves!\n", xtime(), moves);
        fastestFound = moves;
        solutionCount++;

#if 1
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
            if (in[iy][ix] == 0) {
                // Cell is empty, can't do anything
                continue;
            }

            // Up
            if (iy >= 1) {
                if (GE(in,iy-1,ix) || GI(in,iy-1,ix) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(in, tgrd);
                    moveItem(tgrd, ix, iy, ix, iy-1);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }

            // Down
            if (iy <= 2) {
                if (GE(in,iy+1,ix) || GI(in,iy+1,ix) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(in, tgrd);
                    moveItem(tgrd, ix, iy, ix, iy+1);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }

            // Left
            if (ix >= 1) {
                if (GE(in,iy,ix-1) || GI(in,iy,ix-1) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(in, tgrd);
                    moveItem(tgrd, ix, iy, ix-1, iy);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }

            // Right
            if (ix <= 4) {
                if (GE(in,iy,ix+1) || GI(in,iy,ix+1) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(in, tgrd);
                    moveItem(tgrd, ix, iy, ix+1, iy);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }
        }
    }

    return 0;
}

void doSearch() {
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
}

int _tmain(int argc, _TCHAR* argv[])
{
    buildTable();

    startTime = GetTickCount();

    printf("TTL PLRS: %d\n", ttlPlyrs);
    printf("sizeof(GRID): %d\n", sizeof(GRID));
    printf("sizeof(AStarOpen): %d\n", sizeof(AStarOpen));
    printf("\n");

    printf("%d\n", checkWin(winTest));
    printf("%d\n", checkWin(start));
    //printf("%d\n", findClosestItem(distTest, 0, 2, 2, 4));

    addToOpen(start, 0, nullptr);
    doSearch();

    printf("[% 7d] Completed! [%dk %dk %d]\n", xtime(), powCount/1000, poofCount/1000, powCount-poofCount);
    system("PAUSE");
    return 0;
}

