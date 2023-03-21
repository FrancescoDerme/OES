#include <set>
#include <stack>
#include <float.h>
#include <iostream>
 
#define ROW 15
#define COL 15

#define DO_NOTHING 0
#define SUBSTITUTION 1
#define DELETION 2
#define INSERTION 3

// Creating a shortcut for int, int pair type
typedef std::pair<int, int> Pair;
 
// Creating a shortcut for pair<int, pair<int, int>> type
typedef std::pair<double, std::pair<int, int> > pPair;
 
// A structure to hold the parameters of a position on the map
struct cell {

    // Row and Column index of the cell's parent. The parent of a cell 'A' is not necessarily the cell 'B' which generated 'A'.
    // The parent of a cell 'A' is the cell 'B' which is the second-to-last step in the fastes path going from source to 'A'.
    // In other words, the shortest path to a cell always includes it's parent. This way the problem is recursively solved:
    // when every cell has a parent every shortest path has been found.
    
    // Note that 0 <= y-coordinate <= ROW-1 & 0 <= x-coordinate <= COL-1.
    int parent_i, parent_j;
    
    // This is how many steps - more precisely edits - the current best path takes to get from source to this cell.
    int edit_distance;

    // f = g + h where 'g' is the cost to get to this cell from source, 'h' is the cost, esitimated by a heuristic, to get from this cell to destination.
    // This parameter 'h' is what separates Djkastra's algorithm from A*
    double f;
    
    // edit_from_parent represents which edit action was taken last,
    // current_edit is the string that was created by applying the edit actions from source to this cell.
    std::string edit_from_parent, current_edit;
};

void optimalEditSequence(const std::string& s1, const std::string& s2);

void levenshteinDistance(int grid[][COL], const std::string& s1, const std::string& s2, const std::size_t l1, const std::size_t l2);
bool aStarSearch(int grid[][COL], Pair src, Pair dest, int edit_distance, const std::string& s1, const std::string& s2);

inline bool isValid(int row, int col);
inline bool isUnBlocked(int grid[][COL], int row, int col);
inline bool isDestination(int row, int col, Pair dest);
inline std::string applayEdit(const std::string& original, const std::string& source, const std::string& dest, int edit_action, int index);
inline double calculateHValue(const std::string& s1, const std::string& s2, int row, int col, Pair reach);
inline void tracePath(cell cellDetails[][COL], const std::string& s2, Pair dest);

void printMatrix(int matrix[][15], size_t l1, size_t l2);
void printSet(std::set<pPair> openList);

int main() {
    std::string s1 = "NICOLA";   // starting point
    std::string s2 = "NIBBI"; // end goal
    
    optimalEditSequence(s1, s2);
    return 0;
}

void optimalEditSequence(const std::string& s1, const std::string& s2) {
    const std::size_t l1 = s1.size();
    const std::size_t l2 = s2.size();
    
    // Description of the Grid
    // !(-1) -> The cell is not blocked
    // (-1) -> The cell is blocked
    
    int grid[ROW][COL];
    memset(grid, -1, sizeof grid);
    
    levenshteinDistance(grid, s1, s2, l1, l2);
    
    // Source is the left-most top-most corner
    Pair src = std::make_pair(0, 0);
 
    // Destination is the left-most bottom-most corner
    Pair dest = std::make_pair(l1, l2);
    
    bool edit_sequence_found = aStarSearch(grid, src, dest, grid[l1][l2], s1, s2);
    
    if (edit_sequence_found) printf("Therefore their Levenshtein distance is: %d\n", grid[l1][l2]);
    else printf("Optimal path not found, the Levenshtein distance between %s and %s is: %d\n", s1.c_str(), s2.c_str(), grid[l1][l2]);
}

void levenshteinDistance(int grid[][COL], const std::string& s1, const std::string& s2, const std::size_t l1, const std::size_t l2) {
    int max_size = std::max(ROW, COL);
    
    for (unsigned int i = 0; i < max_size; i++) {
        if (i <= l1) grid[i][0] = i;
        else grid[i][0] = -1;
        
        if (i <= l2) grid[0][i] = i;
        else grid[0][i] = -1;
    }
    
    for (unsigned int j = 1; j <= l2; j++) {
        for (unsigned int i = 1; i <= l1; i++) {
            
            // This slight optimisation takes advantage of the fact that when jumping diagonally costs it's always the best option,
            // otherwise the optimal parent cell for [][] is just the one with the lesser value.
            
            if (s1[i-1] == s2[j-1]) grid[i][j] = grid[i-1][j-1];
            else grid[i][j] = std::min({grid[i-1][j], grid[i][j-1], grid[i-1][j-1]}) + 1;
        }
    }
    
    printMatrix(grid, l1, l2);
}

// Find the shortest path from source to destination using A* Search Algorithm
bool aStarSearch(int grid[][COL], Pair src, Pair dest, int edit_distance, const std::string& s1, const std::string& s2) {

    // Either the source or the destination is out of range
    if (isValid(src.first, src.second) == false ||
        isValid(dest.first, dest.second) == false) {
        printf("Source or destination is invalid\n");
        return false;
    }
 
    // Either the source or the destination is blocked
    if (isUnBlocked(grid, src.first, src.second) == false ||
        isUnBlocked(grid, dest.first, dest.second) == false) {
        printf("Source or destination is blocked\n");
        return false;
    }
 
    // If the destination cell is the same as source cell
    if (isDestination(src.first, src.second, dest) == true) {
        printf("We are already at the destination\n");
        return false;
    }
    
    // This among many other algorithms uses 2 lists during the search: Open List and Closed List.
    // Source cell starts in the Open List meaning that it can visit it's neighbors but hasn't done so quite yet.
    // When the source cell visits it's neighbors in a process called "expansion", they all get transferred to the
    // Open List and are referred to as "generated" while source goes into the Closed List and is referred to as "expanded".
    
    // Therefore the Open List is a collection of all nodes that are neighbors of expanded nodes who generated them: they can visit their neighbors but still have to.
    // It is often implemented as a priority queue but here it's a set of pair of pair like <f, <i, j>> where f = g + h, and i, j are the cell's coordinates.
    // This is convenint as per the C++ standard, iteration over the elements in an std::set proceeds in sorted order as determined by std::less which means
    // that openList.begin() will automatically give back the cell with lowest f in the set that is what we want to expand next.
    
    // The Closed List is a collection of all expanded nodes: nodes that were already "searched" meaning that they have already visited their neighbors.
    // This prevents the search from visiting the same nodes again and again. It is often implemented as a boolean 2D array.
    
    // Note that a cell in the Open List might - as bad as it may sound - potentially change parent if a better route is found.
    // This is due to the fact that Euclidean distance, what we're using as 'h', is an admissable heuristic:
    // it can only underestimate the cost to get to the destination and it doesn't guarantee that the first path explored is also the optimal one.
    // A cell in the closed list cannot change parent beacuse it's neighbors were already searched: no better route can be found.
 
    bool closedList[ROW][COL];
    memset(closedList, false, sizeof(closedList));

    std::set<pPair> openList;
    
    std::stack<Pair> Path;

    // Declare a 2D array of cells to hold their f, g and h values
    cell cellDetails[ROW][COL];
    
    {
        int i, j;
     
        for (i = 0; i < ROW; i++) {
            for (j = 0; j < COL; j++) {

                cellDetails[i][j].f = FLT_MAX;
                cellDetails[i][j].parent_i = -1;
                cellDetails[i][j].parent_j = -1;
                cellDetails[i][j].edit_distance = -1;
                cellDetails[i][j].edit_from_parent = "\0";
                cellDetails[i][j].current_edit = s1;
            }
        }
     
        // Initialising the parameters of the starting node: it is it's own parent
        i = src.first;
        j = src.second;

        cellDetails[i][j].f = 0.0;
        cellDetails[i][j].parent_i = i;
        cellDetails[i][j].parent_j = j;
        cellDetails[i][j].edit_distance = 0;
        
        // Put the starting cell on the open list and set its 'f' to 0
        openList.insert(std::make_pair(0.0, std::make_pair(i, j)));
    }

    // While there exist generated cells which haven't visited their neighbors...
    while (!openList.empty()) {

        pPair expanding_cell = *openList.begin();
        
        // Printing out openList at every iteration is incredibly useful to understand how the search works and to debug
        // printSet(openList);
        
        // Remove this vertex from the open list
        openList.erase(expanding_cell);

        // Add this vertex to the closed list
        int i = expanding_cell.second.first;
        int j = expanding_cell.second.second;

        closedList[i][j] = true;

         // If the expanding cell is (i, j) then
         // South          (i+1, j)
         // East           (i, j+1)
         // South-East     (i+1, j+1)
 
        // To store the 'g', 'h' and 'f' of the 3 successors
        double gNew, hNew, fNew;
        std::string edited;
 
        //----------- 1st Successor (South) ------------
 
        // Only process this cell if this is a valid one
        if (isValid(i + 1, j)) {
            
            // If the south successor is the destination
            if (isDestination(i + 1, j, dest)) {

                // Set the Parent of the destination cell. Note that (i + 1 == dest.first) && (j == dest.second)
                cellDetails[i + 1][j].parent_i = i;
                cellDetails[i + 1][j].parent_j = j;
                cellDetails[i + 1][j].edit_from_parent = "deletion";
                cellDetails[i + 1][j].current_edit = s2;

                if (cellDetails[i][j].edit_distance + 1 == edit_distance) { tracePath(cellDetails, s2, dest); return true; }
                else return false;
            }

            // Only process this cell if it isn't in the Closed list and it isn't blocked
            else if (!closedList[i + 1][j] && isUnBlocked(grid, i + 1, j)) {

                gNew = grid[i + 1][j];
                
                edited = applayEdit(s1, cellDetails[i][j].current_edit, s2, DELETION, i + 1);
                hNew = calculateHValue(edited, s2, i + 1, j, dest);
                
                fNew = gNew + hNew;
 
                // If this cell isn’t on the open list, add it there, make the current square its parent and record its f, g, and h costs. If it is
                // on the open list already, check to see if this path is better based on its 'f'. Note that in the first case the cell's f is FLT_MAX

                if (cellDetails[i + 1][j].f > fNew) {

                    openList.insert(std::make_pair(fNew, std::make_pair(i + 1, j)));

                    // Update the details of this cell
                    cellDetails[i + 1][j].f = fNew;
                    cellDetails[i + 1][j].parent_i = i;
                    cellDetails[i + 1][j].parent_j = j;
                    cellDetails[i + 1][j].edit_distance = cellDetails[i][j].edit_distance + 1;
                    cellDetails[i + 1][j].edit_from_parent = "deletion";
                    cellDetails[i + 1][j].current_edit = edited;
                }
            }
        }


        //----------- 2nd Successor (East) ------------
 
        // Only process this cell if this is a valid one
        if (isValid(i, j + 1)) {

            // If the east successor is the destination
            if (isDestination(i, j + 1, dest)) {

                // Set the Parent of the destination cell. Note that (i 1 == dest.first) && (j + 1 == dest.second)
                cellDetails[i][j + 1].parent_i = i;
                cellDetails[i][j + 1].parent_j = j;
                cellDetails[i][j + 1].edit_from_parent = "insertion";
                cellDetails[i][j + 1].current_edit = s2;

                if (cellDetails[i][j].edit_distance + 1 == edit_distance) { tracePath(cellDetails, s2, dest); return true; }
                else return false;
            }

            // Only process this cell if it isn't in the Closed list and it isn't blocked
            else if (!closedList[i][j + 1] && isUnBlocked(grid, i, j + 1)) {

                gNew = grid[i][j + 1];
                
                edited = applayEdit(s1, cellDetails[i][j].current_edit, s2, INSERTION, i);
                hNew = calculateHValue(edited, s2, i, j + 1, dest);
                
                fNew = gNew + hNew;

                // If this cell isn’t on the open list, add it there, make the current square its parent and record its f, g, and h costs. If it is
                // on the open list already, check to see if this path is better based on its 'f'. Note that in the first case the cell's f is FLT_MAX

                if (cellDetails[i][j + 1].f > fNew) {

                    openList.insert(std::make_pair(fNew, std::make_pair(i, j + 1)));
 
                    // Update the details of this cell
                    cellDetails[i][j + 1].f = fNew;
                    cellDetails[i][j + 1].parent_i = i;
                    cellDetails[i][j + 1].parent_j = j;
                    cellDetails[i][j + 1].edit_distance = cellDetails[i][j].edit_distance + 1;
                    cellDetails[i][j + 1].edit_from_parent = "insertion";
                    cellDetails[i][j + 1].current_edit = edited;
                }
            }
        }

        //----------- 3d Successor (South-East)
 
        // Only process this cell if this is a valid one
        if (isValid(i + 1, j + 1)) {

            // If the south-east successor is the destination
            if (isDestination(i + 1, j + 1, dest)) {

                // Set the Parent of the destination cell. Note that (i + 1 == dest.first) && (j + 1 == dest.second)
                cellDetails[i + 1][j + 1].parent_i = i;
                cellDetails[i + 1][j + 1].parent_j = j;
                cellDetails[i + 1][j + 1].edit_from_parent = (grid[i + 1][j + 1] - grid[i][j] == 0) ? "\0" : "substitution";
                cellDetails[i + 1][j + 1].current_edit = s2;
                
                cellDetails[i + 1][j + 1].edit_distance = cellDetails[i][j].edit_distance + ((grid[i + 1][j + 1] - grid[i][j] == 0) ? 0 : 1);
                
                if (cellDetails[i + 1][j + 1].edit_distance == edit_distance) { tracePath(cellDetails, s2, dest); return true; }
                else return false;
            }

            // Only process this cell if it isn't in the Closed list and it isn't blocked
            else if (!closedList[i + 1][j + 1] && isUnBlocked(grid, i + 1, j + 1)) {

                gNew = grid[i + 1][j + 1];
                
                edited = applayEdit(s1, cellDetails[i][j].current_edit, s2, (grid[i + 1][j + 1] - grid[i][j] == 0) ? DO_NOTHING : SUBSTITUTION, i + 1);
                hNew = calculateHValue(edited, s2, i + 1, j + 1, dest);
                
                fNew = gNew + hNew;

                // If this cell isn’t on the open list, add it there, make the current square its parent and record its f, g, and h costs. If it is
                // on the open list already, check to see if this path is better based on its 'f'. Note that in the first case the cell's f is FLT_MAX

                if (cellDetails[i + 1][j + 1].f > fNew) {

                    openList.insert(std::make_pair(fNew, std::make_pair(i + 1, j + 1)));
 
                    // Update the details of this cell
                    cellDetails[i + 1][j + 1].f = fNew;
                    cellDetails[i + 1][j + 1].parent_i = i;
                    cellDetails[i + 1][j + 1].parent_j = j;
                    cellDetails[i + 1][j + 1].edit_distance = cellDetails[i][j].edit_distance + ((grid[i + 1][j + 1] - grid[i][j] == 0) ? 0 : 1);
                    cellDetails[i + 1][j + 1].edit_from_parent = (grid[i + 1][j + 1] - grid[i][j] == 0) ? "\0" : "substitution";
                    cellDetails[i + 1][j + 1].current_edit = edited;
                }
            }
        }
        
    }
 
    // When the destination cell is not found and the open list is empty, then we conclude that we failed to
    // reach the destination cell. This may happen when the there is no way to destination cell (due to blockages)
    printf("Failed to find the Destination Cell\n");
    return false;
}
 
// Check whether given cell (row, col) is a valid,
// returns true if row number and column number is in range.
inline bool isValid(int row, int col) {
    return (row >= 0) && (row < ROW) && (col >= 0) && (col < COL);
}
 
// Check whether the given cell is blocked,
// returns true if the cell is not blocked.
inline bool isUnBlocked(int grid[][COL], int row, int col) {
    if (grid[row][col] != -1) return (true);
    else return (false);
}
 
// Check whether destination cell has been reached
inline bool isDestination(int row, int col, Pair dest) {
    if (row == dest.first && col == dest.second) return (true);
    else return (false);
}

inline std::string applayEdit(const std::string& original, const std::string& source, const std::string& dest, int edit_action, int index) {
    
    // DO_NOTHING 0
    // SUBSTITUTION 1
    // DELETION 2
    // INSERTION 3
    
    const std::size_t l1 = original.length();
    const std::size_t l2 = source.length();
    
    if (l1 > l2) index -= (l1 - l2);
    else if (l2 > l1) index += (l2 - l1);
    
    std::string edit = source;
    
    switch(edit_action) {
        case INSERTION:
            edit.insert(index, dest, index, 1);
            break;
        case DELETION:
            edit.erase(index - 1, 1);
            break;
        case SUBSTITUTION:
            edit.erase(index - 1, 1);
            edit.insert(index - 1, dest, index - 1, 1);
            break;
        default:
            break;
    }
    
    return edit;
}

// Calculate the 'h' heuristic using the distance formula
inline double calculateHValue(const std::string& s1, const std::string& s2, int row, int col, Pair dest) {
    const int l1 = (int) s1.size();
    const int l2 = (int) s2.size();
    int max_size = std::max(l1, l2);
    
    double number_of_errors = 0;
    for (unsigned int i = 0; i < max_size; i++) if (s1[i] != s2[i]) number_of_errors++;
    double number_of_errors_weighted = number_of_errors / 5;
    
    double euler_distance = ((double)sqrt((row - dest.first) * (row - dest.first) + (col - dest.second) * (col - dest.second)));
    double euler_distance_weighted = euler_distance / 5;
    
    return number_of_errors_weighted + euler_distance_weighted;
}
 
// Trace the path from source to destination. The task is trivial because each cell has a parent: ascending the chain of
// parent-son relationships from destination back to the source gives exactly the reverse of the shortest path from source to destination.
inline void tracePath(cell cellDetails[][COL], const std::string& s2, Pair dest) {
    printf("The optimal edit sequence from %s to %s is:\n", cellDetails[0][0].current_edit.c_str(), s2.c_str());

    int row = dest.first;
    int col = dest.second;
    std::stack<std::string> Path;

    // Loop ends after finding a cell that is it's own parent. Source is the only cell that meets such a requirement
    while (!(cellDetails[row][col].parent_i == row && cellDetails[row][col].parent_j == col)) {

        int parent_row = cellDetails[row][col].parent_i;
        int parent_col = cellDetails[row][col].parent_j;
        
        Path.push(cellDetails[row][col].edit_from_parent);
        Path.push(cellDetails[row][col].current_edit);

        row = parent_row;
        col = parent_col;
    }

    // Print result
    while (!Path.empty()) {
        std::string current_edit = Path.top();
        Path.pop();
        
        std::string edit = Path.top();
        Path.pop();
        
        if (edit != "\0") {
            if (current_edit != s2) printf("- %s to get to %s,\n", edit.c_str(), current_edit.c_str());
            else printf("- %s to finally get to %s.\n", edit.c_str(), current_edit.c_str());
        }
    }
    
    std::cout << std::endl;
    return;
}

void printMatrix(int matrix[][15], size_t l1, size_t l2) {
    std::cout << "Levenshtein matrix" << std::endl;
    
    for (unsigned int i = 0; i <= l1; i++) {
        for (unsigned int j = 0; j <= l2; j++) {
            printf("%2d ", matrix[i][j]);
        } std::cout << std::endl;
    } std::cout << std::endl;
}

// I mean isn't the name self-explanatory?
void printSet(std::set<pPair> openList) {
    for (pPair const& cell : openList) std::cout << cell.first << ' ' << cell.second.first << ' ' << cell.second.second << std::endl;
    std::cout << std::endl;
}
