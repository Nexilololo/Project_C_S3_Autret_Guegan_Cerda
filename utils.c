#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

static char *getID(int i)
{
    // translate from 1,2,3, .. ,500+ to A,B,C,..,Z,AA,AB,...
    static char buffer[10];
    char temp[10];
    int index = 0;

    i--; // Adjust to 0-based index
    while (i >= 0)
    {
        temp[index++] = 'A' + (i % 26);
        i = (i / 26) - 1;
    }

    // Reverse the string to get the correct order
    for (int j = 0; j < index; j++)
    {
        buffer[j] = temp[index - j - 1];
    }
    buffer[index] = '\0';

    return buffer;
}

t_cell* createCell(int vertex, float weight) {
    t_cell* cell = (t_cell*)malloc(sizeof(t_cell));
    cell->vertex = vertex;
    cell->weight = weight;
    cell->next = NULL;
    return cell;
}

t_list CreateEmptyList(){
    t_list newList;
    newList.head = NULL;
    return newList;
}

void addcell(t_list* l,t_cell* c){
    if (l->head==NULL) {
        l->head = c;
    }
    else {
        c->next=l->head;
        l->head=c;
    }
};

void display_list(t_list l,int num) {
    printf("\nList for vertex %d:[head @]",num);
    t_cell* curr=l.head;
    if (curr!=NULL) {
        printf(" -> (%d, %.2f)",curr->vertex,curr->weight);
        curr=curr->next;
        while (curr!=NULL) {
            printf(" @-> (%d, %.2f)",curr->vertex,curr->weight);
            curr=curr->next;
        }
    }
}

t_adjlist createEmptyAdjlist(int nbvertex) {
    t_adjlist adj;
    adj.list_number = nbvertex;
    adj.list= (t_list*)calloc(nbvertex, sizeof(t_list));
    return adj;
}

void display_adjlist(t_adjlist adj) {
    for (int i = 0; i < adj.list_number; i++) {
        // i + 1 for the graph later on
        display_list(adj.list[i], i + 1);
    }
}

t_adjlist readGraph(const char *filename) {
    FILE *file = fopen(filename, "rt");
    int nbvert, start, end;
    float proba;

    if (file == NULL)
    {
        perror("Could not open file for reading");
        exit(EXIT_FAILURE);
    }

    if (fscanf(file, "%d", &nbvert) != 1)
    {
        perror("Could not read number of vertices");
        exit(EXIT_FAILURE);
    }

    t_adjlist adjlist = createEmptyAdjlist(nbvert);

    while (fscanf(file, "%d %d %f", &start, &end, &proba) == 3)
    {
        int index = start - 1;

        t_cell* new_cell = createCell(end, proba);

        new_cell->next = adjlist.list[index].head;
        adjlist.list[index].head = new_cell;
    }

    fclose(file);
    return adjlist;
}

void checkMarkovGraph(t_adjlist* adjlist){
    int confirmation = 1;
    for(int i = 0; i < adjlist->list_number; i++){
        double suma = 0;
        t_cell* curr = adjlist->list[i].head;
        if(curr == NULL){
            printf("The %d vertex has no transitions.\n", i + 1);
            confirmation = 0;
            continue;
        }
        while(curr != NULL){
            suma += curr->weight;
            curr = curr->next;
        }
        if(suma < 0.99 || suma > 1.01){
            printf("the sum of the probabilities of vertex %d is %.2f.\n", i + 1, suma);
            confirmation = 0;
        }
    }
    if(confirmation){
        printf("The graph is a Markov graph.\n");
    }
    else{
        printf("The graph is not a Markov graph.\n");
    }
}

void createMermaidFile(t_adjlist graph, char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Error creating output file");
        exit(EXIT_FAILURE);
    }

    fprintf(f, "---\nconfig:\n   layout: elk\n   theme: neo\n   look: neo\n---\n\nflowchart LR\n");

    // Here we declare the vertexes with their respective letters
    for (int i = 0; i < graph.list_number; i++) {
        fprintf(f, "%s((%d))\n", getID(i + 1), i + 1);
    }

    fprintf(f, "\n");

    // Loop to write the edges of the mermaid file
    for (int i = 0; i < graph.list_number; i++) {
        t_cell *current = graph.list[i].head;
        while (current != NULL) {
            fprintf(f, "%s -->|%.2f|",getID(i + 1),current->weight);
            fprintf(f,"%s\n",getID(current->vertex));
            current = current->next;
        }
    }

    fclose(f);
    printf("\nMermaid file successfully created to %s\n", filename);
}

t_partition tarjanAlgorithm(t_adjlist graph) {
    t_partition partition;
    partition.classes_number = 0;
    partition.partition = NULL;

    int num = 0;

    t_tarjan_vertex *vertices = createTarjanList(graph);

    t_stack stack = createEmptyStack();

    for (int i = 0; i < graph.list_number; i++) {
        if (vertices[i].number == -1) {
            parcoursTarjan(i, &num, &graph, &stack, vertices, &partition);
        }
    }

    free(vertices);

    return partition;
}

t_tarjan_vertex * createTarjanList(t_adjlist graph) {
    t_tarjan_vertex* t_list = (t_tarjan_vertex*)calloc(graph.list_number, sizeof(t_tarjan_vertex));
    for (int i = 0; i < graph.list_number; i++) {
        t_list[i].ID = i + 1;
        t_list[i].number = -1;
        t_list[i].access_number = -1;
        t_list[i].indicator = 0;
    }
    return t_list;
}

void parcoursTarjan(int curr, int *num, t_adjlist *graph, t_stack *stack, t_tarjan_vertex *vertex, t_partition *partition) {
    vertex[curr].number = *num;
    vertex[curr].access_number = *num;
    (*num)++;

    push(stack, curr);
    vertex[curr].indicator = 1; // 1 means "is in the stack"

    t_cell *edge = graph->list[curr].head;
    while (edge != NULL) {
        int neighbor_index = edge->vertex - 1;

        if (vertex[neighbor_index].number == -1) {
            // If neighbor has not been visited yet, it recurse
            parcoursTarjan(neighbor_index, num, graph, stack, vertex, partition);

            if (vertex[neighbor_index].access_number < vertex[curr].access_number) {
                vertex[curr].access_number = vertex[neighbor_index].access_number;
            }
        }
        else if (vertex[neighbor_index].indicator == 1) {
            if (vertex[neighbor_index].number < vertex[curr].access_number) {
                vertex[curr].access_number = vertex[neighbor_index].number;
            }
        }
        edge = edge->next;
    }

    if (vertex[curr].access_number == vertex[curr].number) {
        int class_idx = partition->classes_number;
        partition->partition = (t_class *)realloc(partition->partition, (class_idx + 1) * sizeof(t_class));

        partition->partition[class_idx].vertex = NULL;
        partition->partition[class_idx].vertex_number = 0;
        sprintf(partition->partition[class_idx].name, "C%d", class_idx + 1);

        int popped_node_idx;
        do {
            popped_node_idx = pop(stack);
            vertex[popped_node_idx].indicator = 0;

            int v_count = partition->partition[class_idx].vertex_number;
            partition->partition[class_idx].vertex = (t_tarjan_vertex *)realloc(
                partition->partition[class_idx].vertex,
                (v_count + 1) * sizeof(t_tarjan_vertex)
            );

            partition->partition[class_idx].vertex[v_count] = vertex[popped_node_idx];
            partition->partition[class_idx].vertex_number++;

        } while (popped_node_idx != curr);

        partition->classes_number++;
    }
}

int findClassOfVertex(t_partition *partition, int vertexID)
{
    for (int i = 0; i < partition->classes_number; i++) {
        t_class *c = &partition->partition[i];

        for (int j = 0; j < c->vertex_number; j++) {
            if (c->vertex[j].ID == vertexID) {
                return i;
            }
        }
    }
    return -1;
}

void analyzeGraphProperties(t_adjlist graph, t_partition partition) {
    int num_classes = partition.classes_number;
    for(int c = 0; c < num_classes; c++) {
        t_class *cls = &partition.partition[c];
        int is_transient = 0;
        for(int v = 0; v < cls->vertex_number; v++) {
            int vertexID = cls->vertex[v].ID;
            t_cell *curr = graph.list[vertexID - 1].head;
            while(curr != NULL) {
                int dest_class = findClassOfVertex(&partition, curr->vertex);
                if(dest_class != c) {
                    is_transient = 1;
                    break;
                }
                curr = curr->next;
            }
            if(is_transient) break;
        }
        if(is_transient) {
            printf("Class %s is transient\n", cls->name);
        } else {
            printf("Class %s is persistent\n", cls->name);
            if(cls->vertex_number == 1) {
                printf(" -> Vertex %d is absorbing\n", cls->vertex[0].ID);
            }
        }
    }
    if(num_classes == 1)
        printf("The Markov graph is irreducible\n");
    else
        printf("The Markov graph is not irreducible\n");
}