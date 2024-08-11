// Project Identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9

// EECS 281, Project 4 - Zookeeper

#include <iomanip>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <getopt.h>
#include <deque>
#include <queue>   //for PQ
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <cstdlib>
#include "xcode_redirect.hpp"
using namespace std;

enum class VType{
  kNone = 0,
  kWild,
  kBorder,
  kSafe
}; // VType{}

struct Category {
  VType vT = VType::kNone;
};  // Category{}


class Vertex
{
  public:

    Vertex(int xCoord, int yCoord, Category c)
      : xCoordinate(xCoord), yCoordinate(yCoord), category(c){
        visited = false;
        minEWeight = numeric_limits<double>::infinity();
        preVertex = -1;//none before it
      }
    
    Vertex(){
      xCoordinate = 0;
      yCoordinate = 0;
      visited = false;
      minEWeight = numeric_limits<double>::infinity();
      preVertex = -1;//none before it
    }
    
    int getXCoord() const{
      return xCoordinate;
    }

    int getYCoord() const{
      return yCoordinate;
    }

    VType getVType() const{
      return category.vT;
    }

    bool hasVisited(){
      return visited;
    }

    double getMEW() const{
      return minEWeight;
    }

    int getPreVertex() const{
      return preVertex;
    }

    void setDist(double d){
      minEWeight = d;
    }

    void setVisited(){
      visited = true;
    }

    void setPreVex(int p){
      preVertex = p;
    }

  private:
    int xCoordinate;
    int yCoordinate;
    Category category;
    bool visited;
    double minEWeight;
    int preVertex;
};

class Graph
{
    public:
        Graph(size_t num, vector<int> &xCoords, vector<int> &yCoords, int isOnlyMST)
          : numVertices(num), onlyMST(isOnlyMST){
            //myVertices = vertices;
            myVertices.resize(numVertices);
            myPath = {0, 0};
            if(onlyMST == 3)
              myDistances.resize(numVertices, vector<double>(numVertices, numeric_limits<double>::infinity()));
            else 
              myDistances.resize(0);
            xCoordinates= xCoords;
            yCoordinates = yCoords;
            for(size_t i = 0; i < numVertices; ++i){
              Category temp;
              int x = xCoordinates[i];
              int y = yCoordinates[i];
              if(x < 0 && y < 0)
                temp.vT = VType::kWild;
              else if((x == 0 && y <= 0) || (x <= 0 && y == 0))
                temp.vT = VType::kBorder;
              else
                temp.vT = VType::kSafe;
              Vertex v = Vertex(x, y, temp);
              myVertices[i] = v;
            }

            total = 0;
            minIndex = 0;
            bestCost = 0;
            currPath = {0};
            currCost = 0;
            //closingEdge = 0;
            //currEdge = 0;
          }

      void produceMST(size_t startVec){
        //size_t startVec = 0;
        total = 0;
        vector<Vertex> tempVertices = myVertices;
        tempVertices[startVec].setDist(0);
        size_t numTrue = 0;
        for(size_t v = startVec; v < numVertices; ++v){//each vertex in graph
          Vertex* vmin = findMinVertex(tempVertices, startVec);
          vmin->setVisited();
          numTrue++;
          total += vmin->getMEW();
          for(size_t n = startVec; n < numVertices; ++n){//each neighbor of vmin
            Vertex *neighbor = &tempVertices[n];
            double dist = numeric_limits<double>::infinity();
            if(isPossEdge(vmin, neighbor))
              dist = calcDistance(vmin, neighbor);
            if(!neighbor->hasVisited() && dist < neighbor->getMEW()){
              neighbor->setDist(dist);
              neighbor->setPreVex(static_cast<int>(minIndex));
            }
          }  
        }
        if(numTrue == 1 && startVec != (numVertices - 1)){//if one true no edges 
          cerr << "Cannot construct MST" << '\n';
          exit(1);
        }
        cout << total << '\n';
        for(size_t i = 0; i < numVertices; ++i){
          Vertex vTemp = tempVertices[i];
          int prev = vTemp.getPreVertex();
          if(prev != -1){
            cout << min(static_cast<int>(i), prev) << " " << max(static_cast<int>(i), prev) << '\n';
          }
        }
      }

      void prodTempMST(size_t startVec, vector<int> &vec){
        total = 0;
        vector<Vertex> tempVertices;
        size_t numTrue = 0;
        tempVertices.resize(vec.size());
        for(size_t i = 0; i < vec.size(); ++i){
          tempVertices[i] = myVertices[static_cast<size_t>(vec[i])];
        }
        tempVertices[startVec].setDist(0);
        for(size_t v = startVec; v < numVertices; ++v){//each vertex in graph
          Vertex* vmin = findMinVertex(tempVertices, startVec);
          vmin->setVisited();
          numTrue++;
          total += vmin->getMEW();
          for(size_t n = startVec; n < numVertices; ++n){//each neighbor of vmin
            Vertex *neighbor = &tempVertices[n];
            double dist = calcDistance(vmin, neighbor);
            if(!neighbor->hasVisited() && dist < neighbor->getMEW()){
              neighbor->setDist(dist);
              neighbor->setPreVex(static_cast<int>(minIndex));
            }
          }  
        }
        if(numTrue == 1 && startVec != (numVertices - 1)){//if one true no edges 
          cerr << "Cannot construct MST" << '\n';
          exit(1);
        }
      }

      void produceFTSP(){
        size_t edges = 1;
        //vector<Vertex> myVs = myVertices;
        for(size_t n = edges; n < numVertices; ++n){
          size_t pos = findMinEdge(edges, n, myVertices);
          myPath.insert(myPath.begin() + static_cast<int>(pos + 1), static_cast<int>(n));
          ++edges;
        }
        myPath.pop_back();

        //printTSP();
      }

      void produceOTSP(){
        produceFTSP();//(PART B) total is bestCost, myPath is best Path
        bestCost = total;
        currPath = myPath;
        iota(currPath.begin(), currPath.end(), 0);;
        genPerms(1);
        total = bestCost;
        //printTSP();
      }

      bool isPossEdge(Vertex* v1, Vertex* v2){
        if(v1->getVType() == VType::kWild && v2->getVType() == VType::kSafe)
          return false;
        else if(v1->getVType() == VType::kSafe && v2->getVType() == VType::kWild)
            return false;
        return true;
      }
      
      double calcDistance(Vertex* v1, Vertex* v2){
        double x1 = static_cast<double>(v1->getXCoord());
        double y1 = static_cast<double>(v1->getYCoord());
        double x2 = static_cast<double>(v2->getXCoord());
        double y2 = static_cast<double>(v2->getYCoord());
        return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
      }

      Vertex* findMinVertex(vector<Vertex> &myVs, size_t start){
        size_t index = 0;
        double minD = numeric_limits<double>::infinity();
        for(size_t i = start; i < numVertices; i++){
          Vertex *temp = &myVs[i];
          if(!temp->hasVisited() && minD > temp->getMEW()){
            minD = temp->getMEW();
            index = i;
          }
        }
        minIndex = index;
        return &myVs[index];
      }

      size_t findMinEdge(size_t numEdges, size_t k, vector<Vertex> &myVs){
        double minDist = numeric_limits<double>::infinity();
        size_t num1 = 0;
        for(size_t i = 0; i < numEdges; ++i){
          double tempD = calcDistance(&myVs[static_cast<size_t>(myPath[i])], &myVs[k]) + calcDistance(&myVs[k], &myVs[static_cast<size_t>(myPath[i + 1])]) - calcDistance(&myVs[static_cast<size_t>(myPath[i])], &myVs[static_cast<size_t>(myPath[i + 1])]);
          if(minDist > tempD){
            minDist = tempD;
            num1 = i;
          }
        }
        total += minDist;
        return num1;
      }

      void printTSP(){
        cout << total << '\n';
        for(size_t i = 0; i < myPath.size() - 1; ++i){
          cout << myPath[i] << " ";
        }
        cout << myPath[myPath.size() - 1] << '\n';
      }

      void genPerms(size_t permLength) {//path is now currPath
        //double currCost = 0;
        if (permLength == currPath.size()) {
          double closingEdge;
          //if(myDistances[static_cast<size_t>(currPath[0])][static_cast<size_t>(currPath[permLength - 1])] != 0)
          closingEdge = myDistances[static_cast<size_t>(currPath[0])][static_cast<size_t>(currPath[permLength - 1])];
          if(closingEdge == numeric_limits<double>::infinity()){
            closingEdge = calcDistance(&myVertices[static_cast<size_t>(currPath[0])], &myVertices[static_cast<size_t>(currPath[permLength - 1])]);
            myDistances[static_cast<size_t>(currPath[0])][static_cast<size_t>(currPath[permLength - 1])] = closingEdge;
          }
          currCost += closingEdge;
          if(currCost < bestCost){//check if its better
            bestCost = currCost;
            myPath = currPath;
            //cout << "New best cost achieved: " << bestCost << '\n';
          }
          currCost -= closingEdge;//ends up being 0 or really close to 0
          return;
        }  // if ..complete path

        if (!promising(currPath, permLength)) {
          return;
        }  // if ..not promising
        
        for (size_t i = permLength; i < currPath.size(); ++i) {
          double currEdge;
          //if(myDistances[static_cast<size_t>(currPath[permLength])][static_cast<size_t>(currPath[permLength - 1])] != 0)
          
          //double currEdge = calcDistance(&myVertices[static_cast<size_t>(currPath[permLength])], &myVertices[static_cast<size_t>(currPath[permLength - 1])]);
          swap(currPath[permLength], currPath[i]);
          currEdge = myDistances[static_cast<size_t>(currPath[permLength])][static_cast<size_t>(currPath[permLength - 1])];
          if(currEdge == numeric_limits<double>::infinity()){
            currEdge = calcDistance(&myVertices[static_cast<size_t>(currPath[permLength])], &myVertices[static_cast<size_t>(currPath[permLength - 1])]);
            myDistances[static_cast<size_t>(currPath[permLength])][static_cast<size_t>(currPath[permLength - 1])] = currEdge;
          }
          currCost += currEdge;
          genPerms(permLength + 1);
          currCost -= currEdge;//indices of edges will be same think of it w vertices
          swap(currPath[permLength], currPath[i]);
        }  // for ..unpermuted elements
        }  // genPerms()

      bool promising(vector<int> &path, size_t permLength){
        double minADist = numeric_limits<double>::infinity();
        double minCDist = numeric_limits<double>::infinity();
        double tempTotal = 0;
        for(size_t i = permLength; i < myVertices.size(); ++i){
            double tempADist = calcDistance(&myVertices[static_cast<size_t>(path[0])], &myVertices[static_cast<size_t>(path[i])]);
            if(minADist > tempADist)
              minADist = tempADist;
            //if(permLength > 1){
            double tempCDist = calcDistance(&myVertices[static_cast<size_t>(path[permLength - 1])], &myVertices[static_cast<size_t>(path[i])]);
            if(minCDist > tempCDist)
              minCDist = tempCDist;
            //}
            //else 
              //minCDist = 0;
          }
        //tempTotal += minADist;//connecting arms
        //tempTotal += minCDist;
        /*size_t k = path.size() - permLength;
        if(k <= 5)
          return true;*/
          
        prodTempMST(permLength, currPath);//(PART A) gray MST
        double mstCost = total;
        //tempTotal += mstCost;
        
        //tempTotal += currCost;
        tempTotal = minADist + minCDist + mstCost + currCost;
        /*bool promise = false;
        if(tempTotal < bestCost)
          promise = true;
        for (size_t i = 0; i < path.size(); ++i){
          cout << setw(2) << path[i] << ' ';
        }
          cout << setw(4) << permLength << setw(10) << currCost;
          cout << setw(10) << minADist << setw(10) << minCDist;
          cout << setw(10) << mstCost << setw(10) << tempTotal << "  " << promise << '\n';*/
        
        if(tempTotal < bestCost)
          return true;
        return false;
      }

    private:
        size_t numVertices;
        vector<int> xCoordinates;
        vector<int> yCoordinates;
        vector<Vertex> myVertices;
        vector<int> myPath;
        double total;
        size_t minIndex;
        double bestCost;
        vector<int> currPath;
        double currCost;
        int onlyMST;
        vector<vector<double>> myDistances;
        //double closingEdge;
        //double currEdge;
};

enum class Mode {
  kNone = 0,
  kMST,
  kFTSP,
  kOTSP,
};  // Mode{}

struct Options {
  Mode mode = Mode::kNone;
};  // Options{}

void getMode(int argc, char * argv[], Options &options) {
  // These are used with getopt_long()
  opterr = false; // Let us handle all error output for command line options
  int choice;
  int index = 0;
  option long_options[] = {
    { "help",  no_argument,       nullptr, 'h'  },
    { "mode",  required_argument, nullptr, 'm'  },
    { nullptr, 0,                 nullptr, '\0' }
  };  // long_options[]

  while ((choice = getopt_long(argc, argv, "hm:", long_options, &index)) != -1) {
    switch (choice) {
      case 'h':
        cout << "This program simulates zookeeper.\n";
        exit(0);
      case 'm':{
        string arg{optarg};
        if(arg == "MST")
            options.mode = Mode::kMST;
        else if(arg == "FASTTSP")
            options.mode = Mode::kFTSP;
        else if(arg == "OPTTSP")
            options.mode = Mode::kOTSP;
        break;
      }// case 'm'
      default:
        cerr << "Error: invalid option" << endl;
        exit(1);
      }  // switch ..choice
  } // while
}  // getMode()

int main(int argc, char* argv[]) {
    // This should be in all of your projects, speeds up I/O
    ios_base::sync_with_stdio(false);

    xcode_redirect(argc, argv);

    cout << std::setprecision(2); //Always show 2 decimal places
    cout << std::fixed; //Disable scientific notation for large numbers

    // Get the mode from the command line and read in the data
    Options options;
    getMode(argc, argv, options);

    if(options.mode == Mode::kNone){
        cerr << "Mode was not specified!" << '\n';
        exit(1);
    }

    vector<int> xCoordinates;
    vector<int> yCoordinates;
    //vector<Vertex> vertices;
    

    size_t num = 0;
    int xcoord;
    int ycoord;
    size_t index = 0;
    if(!cin.fail()){
        cin >> num;
        xCoordinates.resize(num);
        yCoordinates.resize(num);
        //vertices.resize(num);
    }
    //reads in vertices
    while(index < num){//theres more to read in commands file
        cin >> xcoord;
        cin >> ycoord;
        xCoordinates[index] = xcoord;
        yCoordinates[index] = ycoord;
        /*Category temp;
        if(xcoord < 0 && ycoord < 0)
          temp.vT = VType::kWild;
        else if((xcoord == 0 && ycoord <= 0) || (xcoord <= 0 && ycoord == 0))
          temp.vT = VType::kBorder;
        else
          temp.vT = VType::kSafe;

        Vertex v = Vertex(xcoord, ycoord, temp);
        vertices[index] = v;*/

        index++;
    }
    //Graph myGraph(num, xCoordinates, yCoordinates);
    if(options.mode == Mode::kMST){//part a
      Graph myGraph = Graph(num, xCoordinates, yCoordinates, 1);
      myGraph.produceMST(0);
      //myGraph.printMST();
    }
    else if(options.mode == Mode::kFTSP){
      Graph myGraph(num, xCoordinates, yCoordinates, 2);
      myGraph.produceFTSP();
      myGraph.printTSP();
    }
    else if(options.mode == Mode::kOTSP){
      Graph myGraph(num, xCoordinates, yCoordinates, 3);
      /*cout << "Path                               PL   curCost     arm 1     arm 2       MST     Total  Promising?" << '\n';
      cout << fixed << showpoint << setprecision(2) << boolalpha;*/
      myGraph.produceOTSP();
      myGraph.printTSP();
    }

    return 0;
}
