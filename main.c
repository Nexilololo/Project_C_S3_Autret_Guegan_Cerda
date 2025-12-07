#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"
#include "hasse.h"
#include "matrix.h"

void setInitialState(t_matrix *pi, int state) {
    for(int j=0; j<pi->cols; j++) pi->data[0][j] = 0.0;
    if(state > 0 && state <= pi->cols) {
        pi->data[0][state - 1] = 1.0;
    }
}

void setUniformDistribution(t_matrix *pi, int *states, int count) {
    for(int j=0; j<pi->cols; j++) pi->data[0][j] = 0.0;
    double prob = 1.0 / count;
    for(int i=0; i<count; i++) {
        int idx = states[i] - 1;
        if(idx >= 0 && idx < pi->cols) {
            pi->data[0][idx] = prob;
        }
    }
}

int main() {

    printf("\n=== PART 1: LOADING GRAPH ===\n");

    t_adjlist g = readGraph("../data/proba.txt");
    createMermaidFile(g, "../data/mermaid_graph.txt");
    printf("Graph loaded. Mermaid file created.\n");

    printf("\n=== PART 2: COMMUNICATING CLASSES ===\n");

    t_partition p = tarjanAlgorithm(g);
    printf("Found %d communicating classes.\n", p.classes_number);

    t_link_array links = getHasseLinks(g, p);
    createHasseMermaid(p, links, "../data/mermaid_hasse.txt");
    printf("Hasse diagram generated (mermaid_hasse.txt).\n");

    printf("\n--- Class Properties ---\n");
    analyzeGraphProperties(g, p);


    // ==========================================
    // PART 3: SIMULATION (Q1 - Q7)
    // ==========================================
    printf("\n=== PART 3: SIMULATION (Q1) ===\n");

    t_matrix M = createTransitionMatrix(g);

    // We treat a vector as a 1xN matrix to reuse your multiplyMatrix function
    t_matrix pi = createZeroMatrix(1);
    // Manually adjust to 1x27 (since createZeroMatrix makes NxN usually)
    freeMatrix(&pi); // clear the NxN one
    pi.rows = 1;
    pi.cols = M.cols;
    pi.data = (double**)calloc(1, sizeof(double*));
    pi.data[0] = (double*)calloc(pi.cols, sizeof(double));

    // --- SETUP FOR QUESTION 1: Start at State 2 ---
    setInitialState(&pi, 2);

    // File for plotting (Q1b)
    FILE *fp = fopen("../data/test.csv", "w");
    if(fp) fprintf(fp, "Step,State2,State5,State12,State21,State25\n"); // Log specific states of interest

    printf("Starting simulation from State 2...\n");

    // Print initial state (n=0)
    // printMatrix(pi);

    t_matrix pi_curr = createZeroMatrix(1); /* Temp Copy */
    pi_curr.rows = 1; pi_curr.cols = pi.cols;
    pi_curr.data = (double**)malloc(sizeof(double*));
    pi_curr.data[0] = (double*)malloc(pi.cols * sizeof(double));
    for(int j=0; j<pi.cols; j++) pi_curr.data[0][j] = pi.data[0][j];

    int steps_to_check[] = {1, 2, 10, 50}; // From Q1a
    int check_idx = 0;

    for (int n = 1; n <= 50; n++) {
        t_matrix pi_next = multiplyMatrix(pi_curr, M);

        // Write to CSV for Plotting
        if(fp) fprintf(fp, "%d,%f,%f,%f,%f,%f\n", n,
        pi_next.data[0][1],  // State 2
        pi_next.data[0][4],  // State 5
        pi_next.data[0][11], // State 12
        pi_next.data[0][20], // State 21
        pi_next.data[0][24]  // State 25
        );

        // Check if this is a step requested in Q1a
        if (check_idx < 4 && n == steps_to_check[check_idx]) {
            printf("\n--- Distribution after n=%d steps ---\n", n);
            printMatrix(pi_next); // Prints the vector
            check_idx++;
        }

        // Update current vector
        for(int j=0; j<pi.cols; j++) pi_curr.data[0][j] = pi_next.data[0][j];
        freeMatrix(&pi_next); // Don't forget to free the temp result!
    }

    if(fp) { fclose(fp); printf("Data saved to simulation_Q1.csv for plotting.\n"); }
    freeMatrix(&pi_curr);


    // ==========================================
    // PART 4: LIMITING DISTRIBUTIONS (Q8c, Q10)
    // ==========================================
    printf("\n=== PART 4: STATIONARY DISTRIBUTIONS ===\n");

    for (int c = 0; c < p.classes_number; c++) {
        t_class *cls = &p.partition[c];

        // 1. Check if Transient
        int is_transient = 0;
        for(int v = 0; v < cls->vertex_number; v++) {
            int vertexID = cls->vertex[v].ID; // Assuming 1-based ID
            t_cell *curr = g.list[vertexID - 1].head; // 0-based index
            while(curr != NULL) {
                int dest_class = findClassOfVertex(&p, curr->vertex);
                if(dest_class != c) {
                    is_transient = 1;
                    break;
                }
                curr = curr->next;
            }
            if(is_transient) break;
        }

        if (is_transient) {
            printf("Class %d is TRANSIENT. Limiting distribution is 0 everywhere.\n", c);
        } else {
            printf("Class %d is PERSISTENT (Absorbing).\n", c);

            // Extract submatrix corresponding to this class
            t_matrix sub = subMatrix(M, p, c);

            // Try to find convergence
            t_matrix current_pow = createZeroMatrix(sub.rows);
            copyMatrix(&current_pow, sub);
            t_matrix next_pow;

            double diff = 1.0;
            int step = 0;
            int max_steps = 1000;

            while (diff > 0.00001 && step < max_steps) {
                next_pow = multiplyMatrix(current_pow, sub);
                diff = diffMatrix(current_pow, next_pow);

                // Cleanup old matrix
                freeMatrix(&current_pow);
                current_pow = next_pow;
                step++;
            }

            if (step >= max_steps) {
                // If it didn't converge, it is likely Periodic (Q8c)
                printf("  -> WARNING: Did not converge after %d steps.\n", max_steps);
                printf("  -> This class is likely PERIODIC (d > 1). \n");
                // To solve this, you would calculate the average of P^t over the period
            } else {
                // It converged (Aperiodic)
                printf("  -> Converged after %d steps. (Aperiodic, d=1)\n", step);
                printf("  -> Stationary Distribution for this class:\n");
                printMatrix(current_pow);
            }

            // Clean up class specific matrices
            freeMatrix(&current_pow); // Note: sub was copied into this logic flow
            freeMatrix(&sub);
        }
    }

    // Cleanup Global
    freeMatrix(&pi);
    freeMatrix(&M);
    // freeGraph(&g); // Assuming you have this
    // freePartition(&p); // Assuming you have this

    return 0;
}