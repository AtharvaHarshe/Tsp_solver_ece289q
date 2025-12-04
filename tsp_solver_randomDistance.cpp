#include <bits/stdc++.h>
using namespace std;
using Clock = chrono::high_resolution_clock;

// --------------------- Tour length ---------------------
double tourLength(const vector<int>& tour,
                  const vector<vector<double>>& dist) {
    double L = 0.0;
    int n = (int)tour.size();
    for (int i = 0; i < n - 1; ++i)
        L += dist[tour[i]][tour[i+1]];
    L += dist[tour[n-1]][tour[0]];
    return L;
}

// --------------------- Build full distance matrix ---------------------
// Input format:
// N
// Node1 Node2 Distance
// Node1 Node2 Distance
// ...
bool loadInstance(const string& fname,
                  vector<vector<double>>& dist) {
    ifstream fin(fname);
    if (!fin) {
        cerr << "Error: cannot open file " << fname << "\n";
        return false;
    }

    int N;
    fin >> N;
    if (!fin) {
        cerr << "Error: failed to read N\n";
        return false;
    }

    // Read header line (e.g. "Node1 Node2 Distance"), but don't rely on it.
    string h1, h2, h3;
    fin >> h1 >> h2 >> h3;

    dist.assign(N, vector<double>(N, 0.0));

    int u, v;
    double w;
    while (fin >> u >> v >> w) {
        // convert 1-based to 0-based
        --u; --v;
        if (u < 0 || v < 0 || u >= N || v >= N) continue;
        dist[u][v] = w;
        dist[v][u] = w;
    }
    return true;
}

// --------------------- Candidate lists (by cheap edges) ---------------------
// For random weights, we build for each node a list of its K cheapest neighbors.
vector<vector<int>> buildCandidates(const vector<vector<double>>& dist,
                                    int K = 40) {
    int N = (int)dist.size();
    vector<vector<int>> cand(N);

    for (int i = 0; i < N; ++i) {
        vector<pair<double,int>> tmp;
        tmp.reserve(N - 1);
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            tmp.push_back({dist[i][j], j});
        }
        sort(tmp.begin(), tmp.end(),
             [](const auto& a, const auto& b){ return a.first < b.first; });
        int useK = min(K, (int)tmp.size());
        cand[i].reserve(useK);
        for (int k = 0; k < useK; ++k)
            cand[i].push_back(tmp[k].second);
    }
    return cand;
}

// --------------------- 2-opt local search (candidate based) ---------------------
bool twoOptDescent(vector<int>& tour,
                   const vector<vector<double>>& dist,
                   const vector<vector<int>>& cand,
                   const Clock::time_point& deadline)
{
    int N = (int)tour.size();
    vector<int> pos(N);
    for (int i = 0; i < N; ++i) pos[tour[i]] = i;

    bool improved = false;

    while (true) {
        bool anyImproved = false;

        for (int i = 0; i < N; ++i) {
            // Time check
            if (Clock::now() >= deadline)
                return improved; // stop gently

            int a = tour[i];
            int iNext = (i + 1 == N) ? 0 : (i + 1);
            int b = tour[iNext];

            const auto& candList = cand[a];
            for (int nb : candList) {
                int j = pos[nb];
                int jNext = (j + 1 == N) ? 0 : (j + 1);
                int c = tour[j];
                int d = tour[jNext];

                // Avoid adjacent edges & trivial reversals
                if (i == j || iNext == j) continue;
                if (j == i || jNext == i) continue;

                double oldCost = dist[a][b] + dist[c][d];
                double newCost = dist[a][c] + dist[b][d];
                if (newCost + 1e-12 < oldCost) {
                    // perform 2-opt between (iNext .. j)
                    if (iNext < j) {
                        reverse(tour.begin() + iNext, tour.begin() + j + 1);
                    } else {
                        // wrap-around case
                        reverse(tour.begin() + jNext, tour.begin() + i + 1);
                    }
                    // update positions
                    for (int k = 0; k < N; ++k)
                        pos[tour[k]] = k;
                    anyImproved = true;
                    improved = true;
                    break; // restart scanning i
                }
            }

            if (anyImproved)
                break;
        }

        if (!anyImproved)
            break;
    }

    return improved;
}

// --------------------- Random 3-opt sampling ---------------------
// Try a limited number of random 3-opt moves to escape local minima.
// If an improving move is found, apply it and return true.
bool random3OptKick(vector<int>& tour,
                    const vector<vector<double>>& dist,
                    mt19937_64& rng,
                    int maxTries,
                    const Clock::time_point& deadline)
{
    int N = (int)tour.size();
    if (N < 6) return false;

    uniform_int_distribution<int> pick(0, N - 1);

    // position array for quick index lookups
    vector<int> pos(N);
    for (int i = 0; i < N; ++i) pos[tour[i]] = i;

    for (int attempt = 0; attempt < maxTries; ++attempt) {
        if (Clock::now() >= deadline) return false;

        int i = pick(rng);
        int j = pick(rng);
        int k = pick(rng);
        // ensure distinct and ordered indices
        vector<int> idx = {i,j,k};
        sort(idx.begin(), idx.end());
        i = idx[0]; j = idx[1]; k = idx[2];

        // need at least spacing
        if (i == j || j == k) continue;
        int i1 = i;
        int j1 = j;
        int k1 = k;

        int a = tour[i1];
        int b = tour[(i1 + 1) % N];
        int c = tour[j1];
        int d = tour[(j1 + 1) % N];
        int e = tour[k1];
        int f = tour[(k1 + 1) % N];

        double base = dist[a][b] + dist[c][d] + dist[e][f];

        // We will test a few standard 3-opt reconnections.
        // For simplicity, just implement 3 variants.

        // Option 1: reconnect (a-c), (b-e), (d-f)
        double opt1 = dist[a][c] + dist[b][e] + dist[d][f];

        // Option 2: (a-d), (e-b), (c-f)
        double opt2 = dist[a][d] + dist[e][b] + dist[c][f];

        // Option 3: (a-e), (d-b), (c-f)
        double opt3 = dist[a][e] + dist[d][b] + dist[c][f];

        double bestDelta = 0.0;
        int bestType = 0;

        if (opt1 + 1e-12 < base) {
            bestDelta = opt1 - base;
            bestType = 1;
        }
        if (opt2 + 1e-12 < base && opt2 - base < bestDelta) {
            bestDelta = opt2 - base;
            bestType = 2;
        }
        if (opt3 + 1e-12 < base && opt3 - base < bestDelta) {
            bestDelta = opt3 - base;
            bestType = 3;
        }

        if (bestType == 0) continue; // no improvement

        // Apply chosen reconnection by rearranging segments
        // segments: [i1+1 .. j1], [j1+1 .. k1], etc.
        vector<int> newTour;
        newTour.reserve(N);

        auto addSegment = [&](int from, int to, bool rev) {
            // inclusive indices, possibly wrap-around
            vector<int> tmp;
            if (from <= to) {
                for (int x = from; x <= to; ++x)
                    tmp.push_back(tour[x]);
            } else {
                for (int x = from; x < N; ++x) tmp.push_back(tour[x]);
                for (int x = 0; x <= to; ++x) tmp.push_back(tour[x]);
            }
            if (rev) reverse(tmp.begin(), tmp.end());
            newTour.insert(newTour.end(), tmp.begin(), tmp.end());
        };

        // We'll only handle non-wrap segments assuming 0 <= i1 < j1 < k1 < N
        // which is true by construction. That simplifies implementation.
        //
        // segments in order: [0..i1], [i1+1..j1], [j1+1..k1], [k1+1..N-1]
        // We'll handle each option by choosing how to permute/reverse these segments.

        const int A0 = 0, A1 = i1;
        const int B0 = i1+1, B1 = j1;
        const int C0 = j1+1, C1 = k1;
        const int D0 = k1+1, D1 = N-1;

        newTour.clear();
        if (bestType == 1) {
            // Option 1 pattern:
            // Keep [0..i1], then reverse [B], then [C], then [D]
            addSegment(A0, A1, false);
            addSegment(B0, B1, true);
            addSegment(C0, C1, false);
            addSegment(D0, D1, false);
        } else if (bestType == 2) {
            // Option 2:
            // [0..i1], [C rev], [B], [D]
            addSegment(A0, A1, false);
            addSegment(C0, C1, true);
            addSegment(B0, B1, false);
            addSegment(D0, D1, false);
        } else if (bestType == 3) {
            // Option 3:
            // [0..i1], [C], [B rev], [D]
            addSegment(A0, A1, false);
            addSegment(C0, C1, false);
            addSegment(B0, B1, true);
            addSegment(D0, D1, false);
        }

        if ((int)newTour.size() == N) {
            tour.swap(newTour);
            return true;
        }
    }

    return false;
}

// --------------------- Generate random permutation tour ---------------------
vector<int> randomTour(int N, mt19937_64& rng) {
    vector<int> t(N);
    iota(t.begin(), t.end(), 0);
    shuffle(t.begin(), t.end(), rng);
    return t;
}

// --------------------- Main ---------------------
int main(int argc, char** argv) {
    auto global_start = Clock::now();
    cout << "Disclamer: Works only for random distances (Eculidian distsnces may give incorrect answer)";
    // Hard total limit: 55 seconds
    const double TOTAL_LIMIT = 55.0;
    auto deadline =
        global_start + chrono::duration_cast<Clock::duration>(chrono::duration<double>(TOTAL_LIMIT));

    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <input_file>\n";
        cout << "Example:\n";
        cout << "    " << argv[0] << " TSP_1000_randomDistance.txt\n";
        return 1;
    }

    string fname = argv[1];

    vector<vector<double>> dist;
    if (!loadInstance(fname, dist)) {
        return 1;
    }
    int N = (int)dist.size();
    cerr << "Loaded instance with N = " << N << "\n";

    // RNG
    mt19937_64 rng(
        (uint64_t)chrono::steady_clock::now().time_since_epoch().count()
    );

    // Build candidate lists once
    auto cand = buildCandidates(dist, 40);

    // ---- Initial solution: simple nearest-neighbor from random start ----
    // (We do a few NN starts and keep the best one before local search)
    int NN_TRIALS = 10;
    vector<int> bestTour;
    double bestCost = 1e18;

    for (int t = 0; t < NN_TRIALS; ++t) {
        if (Clock::now() >= deadline) break;

        uniform_int_distribution<int> pickStart(0, N-1);
        int start = pickStart(rng);

        vector<bool> used(N,false);
        vector<int> T; T.reserve(N);
        int cur = start;
        used[cur] = true;
        T.push_back(cur);

        for (int k = 1; k < N; ++k) {
            double best = 1e18;
            int nxt = -1;
            for (int j = 0; j < N; ++j) {
                if (!used[j] && dist[cur][j] < best) {
                    best = dist[cur][j];
                    nxt = j;
                }
            }
            cur = nxt;
            used[cur] = true;
            T.push_back(cur);
        }

        double cost = tourLength(T, dist);
        if (cost < bestCost) {
            bestCost = cost;
            bestTour = T;
        }
    }

    // Local search on this best NN tour
    if (!bestTour.empty()) {
        twoOptDescent(bestTour, dist, cand, deadline);
        // Optional few 3-opt kicks and re-2-opt
        for (int k = 0; k < 3; ++k) {
            if (Clock::now() >= deadline) break;
            if (random3OptKick(bestTour, dist, rng, 200, deadline)) {
                twoOptDescent(bestTour, dist, cand, deadline);
            } else break;
        }
        bestCost = tourLength(bestTour, dist);
    } else {
        // Fallback: random tour
        bestTour = randomTour(N, rng);
        bestCost = tourLength(bestTour, dist);
    }

    cerr << "Initial best cost: " << bestCost << "\n";

    // ---- Multi-start loop ----
    int iterations = 0;
    while (Clock::now() < deadline) {
        vector<int> T = randomTour(N, rng);

        twoOptDescent(T, dist, cand, deadline);
        if (Clock::now() >= deadline) break;

        // a couple of random 3-opt kicks
        for (int k = 0; k < 2; ++k) {
            if (Clock::now() >= deadline) break;
            if (random3OptKick(T, dist, rng, 300, deadline)) {
                twoOptDescent(T, dist, cand, deadline);
            } else {
                break;
            }
        }

        double cost = tourLength(T, dist);
        if (cost < bestCost) {
            bestCost = cost;
            bestTour = T;
            cerr << "Improved: " << bestCost
                 << " at t = "
                 << chrono::duration<double>(Clock::now() - global_start).count()
                 << " sec\n";
        }

        ++iterations;
    }

    

    // ------------- Output -------------
    cout << fixed << setprecision(6);
   cout << "\nBest length: " << fixed << setprecision(2) << bestCost << "\n";

    
    cout << "Iterations: " ;
    std::cout << std::scientific << std::setprecision(2) << static_cast<double>(iterations) << std::endl<<'\n';
    cout << "Tour (1-based indices):\n";
    for (int i = 0; i < N; ++i) {
        cout << (bestTour[i] + 1);
        if (i + 1 < N) cout << ",";
    }
   cout<< "," << bestTour[0]+1<<'\n';

    double total_time =
        chrono::duration<double>(Clock::now() - global_start).count();
    
    cout << "Total runtime (s): " << total_time << "\n";

    return 0;
}
