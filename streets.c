#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"

struct node {
    double lat;
    double lon;
    int id;
    int num_ways;
    int *way_ids;
};

struct way {
    int id;
    char *name;
    float maxspeed;
    bool oneway;
    int num_nodes;
    int *node_ids;
};

struct ssmap {
    struct node **nodes;
    struct way **ways;
    int nr_nodes;   // Max number of nodes
    int nr_ways;    // Max number of ways
};

struct ssmap * 
ssmap_create(int nr_nodes, int nr_ways)
{
    if (nr_nodes == 0 || nr_ways == 0) {
        return NULL;
    }

    struct ssmap *new_ssmap = (struct ssmap *) malloc(sizeof(struct ssmap));

    if (new_ssmap == NULL) {
        return NULL;
    }

    // Allocate memory for nodes and ways arrays
    new_ssmap->nodes = (struct node **) malloc(sizeof(struct node *) * nr_nodes);
    new_ssmap->ways = (struct way **) malloc(sizeof(struct way *) * nr_ways);

    if (new_ssmap->nodes == NULL || new_ssmap->ways == NULL) {
        // Discard any partially created object/array
        free(new_ssmap->nodes);
        free(new_ssmap->ways);
        free(new_ssmap);
        return NULL;
    }

    for (int i = 0; i < nr_nodes; i++) {
        new_ssmap->nodes[i] = NULL;
    }

    for (int i = 0; i < nr_ways; i++) {
        new_ssmap->ways[i] = NULL;
    }

    // Take number of nodes and ways that are expected to be loaded
    new_ssmap->nr_nodes = nr_nodes;
    new_ssmap->nr_ways = nr_ways;

    return new_ssmap;
}

bool
ssmap_initialize(struct ssmap * m)
{
    return true;
}

void
ssmap_destroy(struct ssmap * m)
{
    for (int i = 0; i < m->nr_nodes; i++) {
        free(m->nodes[i]->way_ids);
        free(m->nodes[i]);
    }
    free(m->nodes);

    for (int i = 0; i < m->nr_ways; i++) {
        free(m->ways[i]->node_ids);
        free(m->ways[i]->name);
        free(m->ways[i]);
    }
    free(m->ways);
    free(m);
}

struct way * 
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed, bool oneway, 
              int num_nodes, const int node_ids[num_nodes])
{
    struct way * new_way = (struct way *) malloc(sizeof(struct way));

    if (new_way == NULL) {
        return NULL;
    }

    new_way->id = id;

    new_way->name = (char *) malloc(strlen(name) + 1); // Space for length + null
    if (new_way->name == NULL) {
        free(new_way);
        return NULL;
    }
    strcpy(new_way->name, name); 

    new_way->maxspeed = maxspeed;
    new_way->oneway = oneway;
    new_way->num_nodes = num_nodes;

    new_way->node_ids = (int *) malloc(num_nodes * sizeof(int));
    if (new_way->node_ids == NULL) {
        free(new_way->name);
        free(new_way);
        return NULL;
    }

    for (int i = 0; i < num_nodes; i++) {
        new_way->node_ids[i] = node_ids[i];
    }

    m->ways[id] = new_way;

    return new_way;
}

struct node * 
ssmap_add_node(struct ssmap * m, int id, double lat, double lon, 
               int num_ways, const int way_ids[num_ways])
{
    struct node * new_node = (struct node *) malloc(sizeof(struct node));
    if (new_node == NULL) {
        return NULL;
    }

    new_node->id = id;
    new_node->lat = lat;
    new_node->lon = lon;
    new_node->num_ways = num_ways;

    new_node->way_ids = (int *) malloc(num_ways * sizeof(int));
    if (new_node->way_ids == NULL) {
        free(new_node);
        return NULL;
    }

    for (int i = 0; i < num_ways; i++) {
        new_node->way_ids[i] = way_ids[i];
    }

    m->nodes[id] = new_node;

    return new_node;
}

void
ssmap_print_way(const struct ssmap *m, int id)
{
    int nr_ways = m->nr_ways;   
    // Iterate over all ways
    for (int i = 0; i < nr_ways; i++) {
        if (m->ways[i]->id == id) {
            printf("Way %d: %s\n", id, m->ways[i]->name);
            return;
        }
    }
    printf("error: way %d does not exist.\n", id);
    return;
}

void
ssmap_print_node(const struct ssmap * m, int id)
{
    int nr_nodes = m->nr_nodes; 
    // Iterate over all nodes
    for (int i = 0; i < nr_nodes; i++) {
        if (m->nodes[i]->id == id) {
            printf("Node %d: (%.7lf, %.7lf)\n", id, m->nodes[i]->lat, m->nodes[i]->lon);
            return;
        }
    }
    printf("error: node %d does not exist.\n", id);
    return;
}

void 
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{
    int nr_ways = m->nr_ways;
    // Iterate over all ways
    for (int i = 0; i < nr_ways; i++) {
        // If there is a match between the name of way and name parameter
        if (strstr(m->ways[i]->name, name) != NULL) {
            printf("%d ", m->ways[i]->id);
        }
    }
    printf("\n");
    return;
}


void 
ssmap_find_node_by_names(const struct ssmap * m, const char * name1, const char * name2)
{
    if (name2 == NULL) {
        // Iterate through all nodes
        for (int i = 0; i < m->nr_nodes; i++) {
            // Iterate through all ways in current node
            for (int j = 0; j < m->nodes[i]->num_ways; j++) {
                int current_way_id = m->nodes[i]->way_ids[j];
                struct way *current_way = m->ways[current_way_id];
                // If there is a match between the name of way and name1 parameter
                if (strstr(current_way->name, name1) != NULL) {
                    printf("%d ", m->nodes[i]->id);
                    break;
                }
            }
        }
    }
    else {
        // Iterate through all nodes
        for (int i = 0; i < m->nr_nodes; i++) {
            bool name1_exists = false;
            bool name2_exists = false;
            // Iterate through all ways in current node
            for (int j = 0; j < m->nodes[i]->num_ways; j++) {
                int current_way_id = m->nodes[i]->way_ids[j];
                struct way *current_way = m->ways[current_way_id];
                // If there is a match between the name of way and name1 parameter
                if (strstr(current_way->name, name1) != NULL && !(name1_exists)) {
                    name1_exists = true;
                }
                // If there is a match between the name of way and name2 parameter
                else if (strstr(current_way->name, name2) != NULL) {
                    name2_exists = true;
                }
                // If there is a match between the name of way and both parameters
                if (name1_exists && name2_exists) {
                    printf("%d ", (m->nodes[i]->id));
                    break;
                }
            }
        }
    }
    printf("\n");
    return;
}

    
/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)

/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) {
    double R = 6371.;       
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1); 
    double dlon = d2r(lon2-lon1); 
    double a = pow(sin(dlat/2), 2) + cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    return R * c; 
}


double 
ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size]) {

    double total_travel_time = 0.0;

    // For keeping track of visited nodes
    bool visited[m->nr_nodes];
 
    for (int i = 0; i < m->nr_nodes; i++) {
        visited[i] = false;
    }
 
    // Check that the nodes in the list have valid node ids or have any duplicates
    for (int i = 0; i < size; i++) {
        if (node_ids[i] < 0 || node_ids[i] >= m->nr_nodes) {
            // Error 1: Each node id specified by the user must be valid
            printf("error: node %d does not exist.\n", node_ids[i]);
            return -1.0;
        }
        if (visited[node_ids[i]]) {
            // Error 5: Should not possible to see the same node twice in a reasonable path
            printf("error: node %d appeared more than once.\n", node_ids[i]);
            return -1.0;
        }
        visited[node_ids[i]] = true;
    }
    
    for (int i = 0; i < size - 1; i++) {
        int current_id = node_ids[i];
        int next_id = node_ids[i + 1];
        double current_maxspeed = 0.0;

        int way_found = false;
        int is_adjacent = false;
        int no_reverse = false;

        // Iterate through all nodes associated with current node
        for (int j = 0; j < m->nodes[current_id]->num_ways; j++) {
            int current_way_id = m->nodes[current_id]->way_ids[j];
            struct way *current_way = m->ways[current_way_id];

            // Iterate through all nodes associated with current way
            for (int k = 0; k < current_way->num_nodes - 1; k++) {
                if (current_way->node_ids[k] == next_id) {
                    way_found = true;
                }

                // If current node id and next node id are adjacent in current way
                if ((current_way->node_ids[k] == current_id && current_way->node_ids[k + 1] == next_id) ||
                    (current_way->node_ids[k] == next_id && current_way->node_ids[k + 1] == current_id)) {
                    is_adjacent = true;

                    if (!current_way->oneway || (current_way->oneway && (current_way->node_ids[k] == current_id && current_way->node_ids[k + 1] == next_id))) {
                        no_reverse = true;
                        current_maxspeed = current_way->maxspeed;
                    }
                }
            }
            
            if (current_way->node_ids[current_way->num_nodes - 1] == next_id) {
                way_found = true;
            }
        }

        if (!way_found) {
            // Error 2: Must be a way object between two adjacent nodes in the array
            printf("error: there are no roads between node %d and node %d.\n", current_id, next_id);
            return -1.0;
        }

        if (!is_adjacent) {
            // Error 3: Two adjacent nodes in a path must also be adjacent in the way object they belong to
            printf("error: cannot go directly from node %d to node %d.\n", current_id, next_id);
            return -1.0;
        }

        if (!no_reverse) {
            // Error 4: If a way object is one-way, the adjacent nodes in the path sharing that way object must also be in the same order
            printf("error: cannot go in reverse from node %d to node %d.\n", current_id, next_id);
            return -1.0;
        }

        double distance = distance_between_nodes(m->nodes[current_id], m->nodes[next_id]);
        double time = distance / current_maxspeed;
        total_travel_time += time;
    }

    return total_travel_time * 60;
}

// Defined structure of nodes in min heap priority queue
struct min_heap_node {
    int node_id;
    int prev_node_id;
    double travel_time;
};

// Defined structure of min heap priority queue
struct min_heap {
    struct min_heap_node **array;
    int current_size;  
};


// Helper function to swap two min heap node struct
void swap_nodes(struct min_heap_node **a, struct min_heap_node **b) {
    struct min_heap_node *temp = *a;
    *a = *b;
    *b = temp;
}

// Helper function that operates min-heapify up which is used when inserting a new min heap node
void min_heapify_up(struct min_heap *queue, int index) {
    while (index != 0 && queue->array[(index - 1) / 2]->travel_time > queue->array[index]->travel_time) {
        swap_nodes(&(queue->array[index]), &(queue->array[(index - 1) / 2]));
        // Index becomes parent
        index = (index - 1) / 2;
    }
}

// Helper function that operates min-heapify down which is used when extracting the root node
void min_heapify_down(struct min_heap* queue, int index) {
    int smallest = index;
    int left_child = 2 * index + 1;
    int right_child = 2 * index + 2;

    if (left_child < queue->current_size && queue->array[left_child]->travel_time < queue->array[smallest]->travel_time)
        smallest = left_child;

    if (right_child < queue->current_size && queue->array[right_child]->travel_time < queue->array[smallest]->travel_time)
        smallest = right_child;

    if (smallest != index) {
        // Maintain min heap property
        swap_nodes(&queue->array[index], &queue->array[smallest]);
        min_heapify_down(queue, smallest);
    }
}

// Helper function to insert a new min heap node
bool insert(struct min_heap *queue, int node_id, int prev_node_id, double travel_time) {
    struct min_heap_node *new_node = malloc(sizeof(struct min_heap_node));
    if (new_node == NULL) {
        return false;
    }

    new_node->node_id = node_id;
    new_node->prev_node_id = prev_node_id;
    new_node->travel_time = travel_time;

    queue->array[queue->current_size] = new_node;
    queue->current_size++;

    // Maintain min heap property
    min_heapify_up(queue, queue->current_size - 1);
    return true;
}


// Helper function to extract and return the root
struct min_heap_node *extract_min (struct min_heap *queue) {
    if (queue->current_size == 0){
        return NULL;
    }
    if (queue->current_size == 1) {
        queue->current_size--;
        return queue->array[0];
    }

    struct min_heap_node *root = queue->array[0];
    queue->array[0] = queue->array[queue->current_size - 1];
    queue->current_size--;

    // Maintain min heap property
    min_heapify_down(queue, 0);

    return root;
}


// Helper function for printing path
void print_sequence(const int sequence[], int end_id, int size) {
    int path[size];
    path[0] = end_id;

    int current_size = 1;
    int current_id = end_id;

    while (sequence[current_id] != current_id) {
        if (current_size >= size) {
            return;
        }
        current_id = sequence[current_id];
        path[current_size++] = current_id;
    }

    for (int i = current_size - 1; i >= 0; i--) {
        printf("%d ", path[i]);
    }
    printf("\n");
    return;
}


void 
ssmap_path_create(const struct ssmap * m, int start_id, int end_id)
{  

    if (start_id < 0 || start_id >= m->nr_nodes) {
        printf("error: node %d does not exist.\n", start_id);
        return;
    }

    if (end_id < 0 || end_id >= m->nr_nodes) {
        printf("error: node %d does not exist.\n", end_id);
        return;
    }

    struct min_heap *queue = malloc(sizeof(struct min_heap));
    if (queue == NULL) {
        printf("error: out of memory.\n");
        return;
    }

    queue->array = malloc(m->nr_nodes * sizeof(struct min_heap_node *));
    if (queue->array == NULL) {
        free(queue);
        printf("error: out of memory.\n");
        return;
    }

    queue->current_size = 0;

    if (!insert(queue, start_id, start_id, 0.0)) {
        printf("error: out of memory.\n");
        return;
    }

    int sequence[m->nr_nodes];
    for (int i = 0; i < m->nr_nodes; i++) {
        sequence[i] = -1;
    }
    

    while (queue->current_size > 0) {
        struct min_heap_node *current_root = extract_min(queue);

        if (sequence[current_root->node_id] == -1) {
            sequence[current_root->node_id] = current_root->prev_node_id;
            struct node *current_node = m->nodes[current_root->node_id];

            // Iterate through all ways associated with current node
            for (int index = 0; index < current_node->num_ways; index++) {
                struct way *current_way = m->ways[current_node->way_ids[index]];

                // Iterate through all nodes associated with current way
                for (int node_index = 0; node_index < current_way->num_nodes; node_index++) {

                    if (current_root->node_id == current_way->node_ids[node_index]) {
                        // Used for checking the neighbours of current node
                        int left_index = node_index - 1;
                        int right_index = node_index + 1;

                        if (node_index != 0 && !(current_way->oneway) && sequence[current_way->node_ids[left_index]] == -1) {
                            double travel_time = distance_between_nodes(current_node, m->nodes[current_way->node_ids[left_index]]) / current_way->maxspeed;
                            double cumulative_travel_time = current_root->travel_time + travel_time;
                            // Insert the node to the priority queue
                            if (!insert(queue, current_way->node_ids[left_index], current_way->node_ids[node_index], cumulative_travel_time)) {
                                printf("error: out of memory.\n");
                                return;
                            }
                        }
                        if (node_index != current_way->num_nodes - 1 && sequence[current_way->node_ids[right_index]] == -1) {
                            double travel_time = distance_between_nodes(current_node, m->nodes[current_way->node_ids[right_index]]) / current_way->maxspeed;
                            double cumulative_travel_time = current_root->travel_time + travel_time;
                            // Insert the node to the priority queue
                            if (!insert(queue, current_way->node_ids[right_index], current_way->node_ids[node_index], cumulative_travel_time)) {
                                printf("error: out of memory.\n");
                                return;
                            }
                        }
                    }
                }
            }

        if (current_root->node_id == end_id) {
                print_sequence(sequence, end_id, m->nr_nodes);
                free(current_root);
                // Free remaining min heap nodes in queue
                while (true) {
                    struct min_heap_node *n = extract_min(queue);
                    if (n == NULL) {
                        break;
                    }
                    free(n);
                }
                free(queue->array);
                free(queue);
                return;
            }
        }

        free(current_root);
    }

    printf("error: could not find a path from node %d to node %d.\n", start_id, end_id);
    free(queue->array);
    free(queue);
}