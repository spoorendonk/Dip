#include "solveRelaxed_callback.h"
#include "MCF_Instance.h"

MCF_Instance *instance = NULL;
DecompAlgo *algo = NULL;
DecompApp *app = NULL;

#include <iostream>
#include <climits>
#include <set>
using namespace std;

int minDist(int dist[], bool Set[], int n) //calculate minimum distance
{
   int min = INT_MAX, min_index;
   for (int v = 0; v < n; v++)
      if (Set[v] == false && dist[v] <= min)
         min = dist[v], min_index = v;
   return min_index;
}

void printPath(int parent[], int j)
{

   // Base Case : If j is source
   if (parent[j] == -1)
   {
      printf("%d ", j);
      return;
   }

   printPath(parent, parent[j]);

   printf("%d ", j);
}

void printSol(int dist[], int n, int parent[]) //print the solution
{
   cout << "Vertex Distance from Source" << endl;
   for (int i = 0; i < n; i++)
   {
      cout << " \t\t"
           << i << " \t\t " << dist[i] << endl;
      printPath(parent, i);
      cout << endl;
   }
}

void dijkstra(std::vector<std::vector<int>> const &g, int src, int dist[], int parent[], int n)
{
   bool Set[n];
   for (int i = 0; i < n; i++)
      dist[i] = INT_MAX, Set[i] = false, parent[src] = -1;
   dist[src] = 0;
   for (int c = 0; c < n - 1; c++)
   {
      int u = minDist(dist, Set, n);
      Set[u] = true;
      for (int v = 0; v < n; v++)
         if (!Set[v] && g[u][v] > -1 && dist[u] != INT_MAX && dist[u] + g[u][v] < dist[v])
         {
            dist[v] = dist[u] + g[u][v];
            parent[v] = u;
         }
   }
}

DecompSolverStatus MyRelaxedSolver(const DecompApp *app,
                                   const int whichBlock,
                                   const double *redCostX,
                                   const double target,
                                   DecompVarList &varList)
{
   // hacky way to find out if we are being called from GenerateInitVars
   bool inGenerateInitVarsPhase = target == 9e15;

   int numArcs = instance->m_numArcs;

   // if in generating init vars return if redCosts are negative, we cannot handle negative cycles
   if (inGenerateInitVarsPhase)
   {
      for (int a = 0; a < numArcs; a++)
      {
         DecompConstraintSet *constrSet = app->m_modelRelax.find(whichBlock)->second.getModel();
         int index = constrSet->activeColumns[a];
         if (redCostX[index] < 0)
            return DecompSolverStatus::DecompSolStatOptimal;
      }
   }

   MCF_Instance::commodity *commodities = instance->m_commodities;
   int source = commodities[whichBlock].source;
   int sink = commodities[whichBlock].sink;

   // map arcs to graph
   int n = instance->m_numNodes;
   int dist[n];
   int parent[n];
   std::fill(dist, dist + n, 0);
   std::fill(parent, parent + n, 0);
   std::vector<std::vector<int>> g(n, std::vector<int>(n, -1));

   // arcs to be integers, 1e-5 precision
   double factor = 10000.0;

   double aj_target = target * factor;
   double aj_source = inGenerateInitVarsPhase ? 0 : aj_target;

   map<tuple<int, int>, int> arc_map;
   MCF_Instance::arc *arcs = instance->m_arcs;
   for (int a = 0; a < numArcs; a++)
   {
      int tail = arcs[a].tail;
      int head = arcs[a].head;

      DecompConstraintSet *constrSet = app->m_modelRelax.find(whichBlock)->second.getModel();
      int index = constrSet->activeColumns[a];
      g[tail][head] = int(redCostX[index] * factor);

      if (!inGenerateInitVarsPhase && tail == source)
         g[tail][head] += aj_source;

      arc_map.insert({make_tuple(tail, head), index});
   }

   // solve shortest path problem
   dijkstra(g, source, dist, parent, n);

   // return if no negative path could be found
   if (dist[sink] - aj_source >= aj_target)
      return DecompSolverStatus::DecompSolStatOptimal;

   // get arcs of path
   vector<tuple<int, int>> arc_list;
   int i = sink;
   int j = parent[i];
   while (j != -1)
   {
      arc_list.push_back(make_tuple(j, i));
      i = j;
      j = parent[j];
   }

   // create var
   int len = arc_list.size();
   int *ind = new int[len];
   double *els = new double[len];
   double origCost = 0;
   double redCost = dist[sink] / factor - aj_source;
   i = 0;
   for (auto x : arc_list)
   {
      ind[i] = arc_map[x];
      els[i] = 1.0;
      origCost += arcs[i].weight * commodities[whichBlock].demand;
      i++;
   }

   DecompVar *var = new DecompVar(len, ind, els, redCost, origCost);
   var->setBlockId(whichBlock);

   varList.push_back(var);
   return DecompSolverStatus::DecompSolStatOptimal;
}
