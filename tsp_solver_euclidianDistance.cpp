#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <cmath>
#include <iomanip>
using namespace std;
using clock_tp = chrono::high_resolution_clock;

// ---------------- compute tour length ----------------
double tourLen(const vector<int>& T, const vector<vector<double>>& dist) {
    double L = 0.0;
    int n = (int)T.size();
    for (int i = 0; i < n - 1; i++)
        L += dist[T[i]][T[i + 1]];
    L += dist[T[n - 1]][T[0]];
    return L;
}

// ---------------- Best-of-N Nearest Neighbor ----------------
vector<int> bestNN(const vector<vector<double>>& dist,
                   mt19937_64& rng,
                   int trials = 40)
{
    int n = (int)dist.size();
    vector<int> bestTour;
    double bestCost = 1e18;
    uniform_int_distribution<int> pickStart(0, n - 1);

    for (int t = 0; t < trials; t++) {
        int start = pickStart(rng);
        vector<bool> used(n, false);
        vector<int> T;
        T.reserve(n);

        int cur = start;
        used[cur] = true;
        T.push_back(cur);

        for (int k = 1; k < n; k++) {
            double best = 1e18;
            int nxt = -1;
            for (int j = 0; j < n; j++) {
                if (!used[j] && dist[cur][j] < best) {
                    best = dist[cur][j];
                    nxt = j;
                }
            }
            cur = nxt;
            used[cur] = true;
            T.push_back(cur);
        }

        double L = tourLen(T, dist);
        if (L < bestCost) {
            bestCost = L;
            bestTour = T;
        }
    }

    return bestTour;
}

// ------------- Build candidate lists (K nearest neighbors) -------------
vector<vector<int>> buildCandidate(const vector<vector<double>>& dist, int K = 50) {
    int n = (int)dist.size();
    vector<vector<int>> cand(n);

    for (int i = 0; i < n; i++) {
        vector<pair<double, int>> tmp;
        tmp.reserve(n - 1);
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            tmp.push_back({ dist[i][j], j });
        }
        sort(tmp.begin(), tmp.end());
        int useK = min(K, (int)tmp.size());
        cand[i].reserve(useK);
        for (int k = 0; k < useK; k++)
            cand[i].push_back(tmp[k].second);
    }
    return cand;
}

// ---------------- FAST 2-opt (Candidate list based) ----------------
bool twoOptFast(vector<int>& T,
                const vector<vector<double>>& dist,
                const vector<vector<int>>& cand)
{
    int n = (int)T.size();
    bool improved = false;

    // position of each node in tour
    vector<int> pos(n);
    for (int i = 0; i < n; i++) pos[T[i]] = i;

    for (int idx = 0; idx < n; idx++) {
        int a = T[idx];
        int idxNext = (idx + 1 == n) ? 0 : (idx + 1);
        int b = T[idxNext];

        for (int nb : cand[a]) {
            int j = pos[nb];
            int jNext = (j + 1 == n) ? 0 : (j + 1);
            int d = T[jNext];

            if (idx == j || idxNext == j) continue;

            double oldCost = dist[a][b] + dist[nb][d];
            double newCost = dist[a][nb] + dist[b][d];

            if (newCost + 1e-12 < oldCost) {
                // Perform 2-opt reversal between (idxNext .. j)
                if (idxNext < j) {
                    reverse(T.begin() + idxNext, T.begin() + j + 1);
                } else {
                    reverse(T.begin() + jNext, T.begin() + idx + 1);
                }

                // Update positions
                for (int k = 0; k < n; k++)
                    pos[T[k]] = k;

                improved = true;
            }
        }
    }

    return improved;
}

// ------------- Or-opt limited (1- and 2-node relocations) -------------
bool orOptLimited(vector<int>& T,
                  const vector<vector<double>>& dist,
                  int MAX_MOVES = 2000)
{
    int n = (int)T.size();
    bool improved = false;
    int moves = 0;

    for (int seg = 1; seg <= 2; seg++) {
        for (int i = 0; i < n && moves < MAX_MOVES; i++) {
            int segEnd = i + seg - 1;
            if (segEnd >= n) break;

            for (int j = 0; j < n && moves < MAX_MOVES; j++) {
                if (j >= i && j <= segEnd) continue;

                int a = (i - 1 + n) % n;
                int b = (segEnd + 1) % n;
                int A = j;
                int C = (j + 1) % n;

                double oldCost = dist[T[a]][T[i]] +
                                 dist[T[segEnd]][T[b]] +
                                 dist[T[A]][T[C]];

                double newCost = dist[T[a]][T[b]] +
                                 dist[T[A]][T[i]] +
                                 dist[T[segEnd]][T[C]];

                if (newCost + 1e-12 < oldCost) {
                    vector<int> newT;
                    newT.reserve(n);

                    for (int k = 0; k < n; k++) {
                        if (k == (A + 1) % n) {
                            for (int t = i; t <= segEnd; t++)
                                newT.push_back(T[t]);
                        }
                        if (k < i || k > segEnd)
                            newT.push_back(T[k]);
                    }

                    T.swap(newT);
                    improved = true;
                    moves++;
                }
            }
        }
    }

    return improved;
}

// ---------------- Strong Double Bridge Perturbation ----------------
void doubleBridgeStrong(vector<int>& T, mt19937_64& rng) {
    int n = (int)T.size();
    if (n < 8) return;

    uniform_int_distribution<int> dist(0, n - 1);
    int a = dist(rng), b = dist(rng), c = dist(rng), d = dist(rng);
    vector<int> idx = { a,b,c,d };
    sort(idx.begin(), idx.end());
    a = idx[0]; b = idx[1]; c = idx[2]; d = idx[3];

    vector<int> p1(T.begin(), T.begin() + a);
    vector<int> p2(T.begin() + a, T.begin() + b);
    vector<int> p3(T.begin() + b, T.begin() + c);
    vector<int> p4(T.begin() + c, T.begin() + d);
    vector<int> p5(T.begin() + d, T.end());

    vector<int> nt;
    nt.reserve(n);

    // reorder segments to make a big "kick"
    nt.insert(nt.end(), p1.begin(), p1.end());
    nt.insert(nt.end(), p3.begin(), p3.end());
    nt.insert(nt.end(), p2.begin(), p2.end());
    nt.insert(nt.end(), p4.begin(), p4.end());
    nt.insert(nt.end(), p5.begin(), p5.end());

    T.swap(nt);
}

// =====================================================================
//                                MAIN
int main(int argc, char* argv[]) {
    auto global_start = clock_tp::now();  // measure full program time
    cout << "Disclamer: Works only for Eculidian distsnces (random distances may give incorrect answer)";
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // RNG
    mt19937_64 rng(
        (uint64_t)chrono::steady_clock::now().time_since_epoch().count()
    );

    // ---------------- read input file name ----------------
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input_file>\n";
        cout << "Example:\n";
        cout << "    " << argv[0] << " TSP_1000_EcludianDistance.txt\n";
        return 1;
    }

    string fname = argv[1];
    cout << "Reading file: " << fname << "\n";

    ifstream f(fname);
    if (!f) {
        cerr << "Cannot open file: " << fname << "\n";
        return 1;
    }

    
    int N;
    f >> N;
    string h1, h2, h3;
    f >> h1 >> h2 >> h3;

    vector<vector<double>> dist(N, vector<double>(N, 0.0));
    int u, v;
    double w;

    while (f >> u >> v >> w) {
        u--; v--;
        dist[u][v] = w;
        dist[v][u] = w;
    }
    f.close();
    
    // -------- Build candidate lists (once) --------
    auto candList = buildCandidate(dist, 50);

    // -------- Initial tour: best-of-40 NN + local search --------
    vector<int> tour = bestNN(dist, rng, 40);
    twoOptFast(tour, dist, candList);
    orOptLimited(tour, dist);
    twoOptFast(tour, dist, candList);

    double bestL = tourLen(tour, dist);
    vector<int> best = tour;

    cout << "Initial tour length: " << bestL << "\n";

    // -------- ILS Loop (hard time limit ~55s) --------
    const double LIMIT = 55.0;  // seconds for search part
    auto search_start = clock_tp::now();
    int iter = 0;

    while (true) {
        double t = chrono::duration<double>(clock_tp::now() - search_start).count();
        if (t > LIMIT) break;

        vector<int> T;

        // 70% of the time: perturb current best
        // 30% of the time: new NN start (to diversify more)    
        uniform_real_distribution<double> prob(0.0, 1.0);
        if (prob(rng) < 0.7) {
            T = best;
            doubleBridgeStrong(T, rng);
        } else {
            T = bestNN(dist, rng, 5);  // smaller trials for speed
        }

        // Local search
        twoOptFast(T, dist, candList);
        orOptLimited(T, dist);
        twoOptFast(T, dist, candList);

        double L = tourLen(T, dist);

        if (L < bestL) {
            bestL = L;
            best = T;
            cout << "Improved: " << L << "  at t=" << t << " sec\n";
        }

        iter++;
    }

    double search_time = chrono::duration<double>(clock_tp::now() - search_start).count();
    double total_time = chrono::duration<double>(clock_tp::now() - global_start).count();

    cout << "\nBest length: " << fixed << setprecision(2) << bestL << "\n";
    cout << "Search time (ILS only): " << search_time << " sec\n";
    cout << "Total runtime (file read + setup + search + output): " << total_time << " sec\n";
    cout << "Iterations";
    std::cout << std::scientific << std::setprecision(2) << static_cast<double>(iter) << std::endl;
    for (int i = 0; i < N; i++)
        cout << best[i] + 1 <<  ',';
    cout<< best[0]+1<<'\n';

    return 0;
}
