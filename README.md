# TSP Solvers – Euclidean & Random Distance Graphs

This repository contains two C++ programs for solving the **Traveling Salesman Problem (TSP)**:

1. **Euclidean Solver** – optimized for coordinate-based or metric distance matrices  
2. **Random/General Distance Solver** – works for any fully connected weighted graph

Both solvers read an input file, generate a tour, optimize it using local search, and print the best cost & tour path.

---

## 📁 How to Use This Project

### 1. Place files together

Make sure **the solver `.cpp/.exe` file and the input dataset `.txt` file are in the same folder**:

📂 TSP_Project
├ tsp_euclid.cpp / tsp_euclid.exe
├ tsp_random.cpp / tsp_random.exe
├ TSP_1000_euclid.txt
├ RandomGraph_1000.txt
└ README.md




You may use different names, but keep solver + input file together for ease.

---

## 🔧 Compilation (if using `.cpp`)

Open **Command Prompt (Windows)**

### Step 1 — Navigate to the project folder

cd "C:\path\to\TSP_Project"
Example:

cd "C:\Users\xyz\Desktop\tsp\"
Useful commands:

cd        foldername	enter a folder
cd ..	    go up one level
dir	      list files in current directory
cls	      clear screen

Step 2 — Compile using g++

g++ -std=c++17 -O3 tsp_euclid.cpp -o tsp_euclid.exe
g++ -std=c++17 -O3 tsp_random.cpp -o tsp_random.exe
After compiling, you should see .exe files created.

▶ Running the Solvers
A) Euclidean Solver 
Best for coordinate-based or metric distance datasets.


tsp_euclid.exe TSP_1000_euclid.txt
Output example:


Disclamer: Works only for Eculidian distances
Reading file: TSP_1000_euclid.txt
Initial tour length: 412.84
Improved: 392.10 at t=4.12 sec
Improved: 381.44 at t=29.52 sec

Best length: 381.44
Search time: 55.00 sec
Iterations: 3.09e+02
Tour:
5,19,42,71,....
⚠ This program runs for ~1 minute (time-limited).

B) Random/General Distance Solver
Works when distances don't obey geometry (random edge weights, general matrices).

tsp_random.exe RandomGraph_1000.txt
Output example:


Initial cost: 9150
After local search: 6422
Final best tour: 5981


📝 Input File Format (both solvers follow same)

N

Header1 Header2 Header3   ← ignored

u v weight

u v weight

...

Example:

1000

X Y D

1 2 4.52

1 3 6.71

2 3 3.98
