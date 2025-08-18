#ifndef QUEUE_H
#define QUEUE_H

#include <stdlib.h>
#include <stdint.h>

// Define the initial capacity of the queue
#define INITIAL_CAPACITY 16

// Define the Queue structure
typedef struct _queue {
    uint64_t *data;       // Pointer to the dynamic array
    size_t front;       // Index of the front element
    size_t back;        // Index of the back element
    size_t size;        // Number of elements in the queue
    size_t capacity;    // Capacity of the dynamic array
} Queue;

// Function to initialize the queue
void init_queue(Queue *queue);

// Function to initialize the queue
void clear_queue(Queue *queue);

// Function to check if the queue is empty
int is_empty(Queue *queue);

// Function to get the number of elements in the queue
size_t size(Queue *queue);

// Function to enqueue an element to the back of the queue
void enqueue(Queue *queue, uint64_t value);

// Function to dequeue an element from the front of the queue
uint64_t dequeue(Queue *queue);

// Function to get the element at the front of the queue without dequeuing it
uint64_t peek(Queue *queue);

// Function to free the memory allocated for the queue
void free_queue(Queue *queue);

#endif // QUEUE_H
